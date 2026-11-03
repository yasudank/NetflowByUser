# Pointing Optimization for P0 Targets and Multi-Visit Assignment
(`optimize_pointings_p0.py` and the `prefer_visits` option of `run_netflow.py`)

This document describes the purpose, algorithm and usage of [optimize_pointings_p0.py](../src/netflow_by_user/optimize_pointings_p0.py), and the netflow assignment option `prefer_visits` for observations with several patterns (several visits).

---

## 1. Overview

`optimize_pointings_p0` takes a tiling of pointings (e.g. the 7 fields in `hexagons_xmm_lss_flat_centers.ecsv`) as a reference and optimizes the center and position angle (PA) of each field so that the **number of top-priority (P0) targets that can actually get a fiber** is maximized.

* Targets are not counted inside a hexagon; instead, the **targets accessible by each Cobra** are computed with the same instrument model as netflow (black dots, broken Cobras, fiducials).
* A field that violates the guide star constraints (`min_stars_per_cam`, `min_cams_with_stars`) or has a bright star near a broken fiber observes no targets; the number of such fields is minimized first.
* **Several patterns** (K tilings of the same area, observed one after another) are supported. With the total time split into K, the layout first maximizes the number of P0 targets observed at least once, then makes the numbers of visits as even as possible.

---

## 2. Evaluation

### 2.1. One field
For every candidate field (RA, Dec, PA) the following is computed and cached:

1. The number of guide stars on each guide camera AG0–AG5 and whether the constraints are met
2. Whether a bright star lies near a broken fiber
3. The focal plane positions of the targets (within `--target-radius` of the center, priority up to `--max-priority`) and the Cobras that can reach each target

These are computed in parallel with `--processes` worker processes.

### 2.2. Assignable count (maximum flow)
The target → Cobra edges of all fields form the network below, and the number of assignable P0 targets is its maximum flow.

```
source → target (capacity j) → (target, pattern) (capacity 1) → Cobra (capacity 1) → sink
```

* **One fiber per target and pattern**: netflow is run for each pattern separately, so two fibers on the same target in overlapping fields of one pattern still count as one visit. With this constraint the prediction matches netflow.
* Fiber collisions and calibration targets (sky, fluxstd) are ignored, so strictly the numbers are upper bounds. In the XMM-LSS tests they matched netflow exactly for P0.

### 2.3. Score (lexicographic)
Layouts are compared lexicographically by the tuple below (leftmost first):

1. The number of fields violating the guide star (and bright star) constraints (fewer is better)
2. The numbers of P0 targets observed at least 1, 2, …, K times (K: number of patterns)
3. The number of accessible P1–P`max_priority` targets (tie breaker)
4. The uncovered fraction of the reference tiling (smaller is better)

"The number of P0 targets observed at least j times" is T_j − T_{j−1}, where T_j is the maximum flow with target capacity j. The assignable fiber count vectors form a polymatroid, so a single assignment maximizes Σ min(n_t, j) for all j at once; these values are therefore exact for the given fields.

---

## 3. Search Algorithm

### Step 1: Global search
The whole reference tiling is rotated (0–60° in steps of `--global-rot-step`) and shifted (±`--global-shift` in steps of `--global-shift-step`) as a rigid body. Each field may also be rotated by a multiple of 60° (the FoV is hexagonal, so the tiling stays contiguous); these choices are made by coordinate descent.

### Step 2: Local search
Starting from the best `--n-refine` global solutions, each field is moved individually by small amounts (±`--local-shift`, ±`--local-rot`) by coordinate descent. Layouts whose uncovered fraction exceeds `--max-gap` are rejected. Fields that violate the constraints are moved first.

### Step 3: Several patterns
1. Patterns are added one at a time: the existing ones are fixed and the new one is optimized by the two steps above (a pattern given by `--init` starts from that layout without the search).
2. Then each pattern is re-optimized in turn with the others fixed (up to `--cycles` rounds, stopping when nothing improves). The local search starts from the current layout and the best `--cycle-n-refine` global solutions, and the result replaces the current layout only if it is better, so the score never decreases.

---

## 4. Usage

After `pip install -e .`, the `optimize_pointings_p0` command is available (`python src/netflow_by_user/optimize_pointings_p0.py` works too). The targets, Gaia catalog, guide star constraints, etc. are read from the YAML given by `--config` (the same configuration file as `run_netflow`).

### 4.1. One pattern
```bash
optimize_pointings_p0 -c netflow_config_2026Nov.yaml \
    --grid xmm_lss/hexagons_xmm_lss_flat_centers.ecsv \
    -o xmm_lss/hexagons_xmm_lss_flat_centers_p0opt.ecsv \
    --global-shift 0.6 --global-rot-step 10 --processes 20 \
    --add-columns --report-file xmm_lss/p0opt_report.txt --plot xmm_lss/p0opt.png
```

### 4.2. Adding a pattern to an existing one (two-pattern observation)
Fix the existing pattern with `--fixed` and give the new pattern's codes a suffix with `--code-suffix`.
```bash
optimize_pointings_p0 -c netflow_config_2026Nov.yaml \
    --grid xmm_lss/hexagons_xmm_lss_flat_centers.ecsv \
    --fixed xmm_lss/hexagons_xmm_lss_flat_centers_p0opt.ecsv --code-suffix _B \
    -o xmm_lss/hexagons_xmm_lss_flat_centers_p0opt_B.ecsv \
    --global-shift 0.6 --global-rot-step 10 --n-refine 3 --processes 20 \
    --add-columns --report-file xmm_lss/p0opt_B_report.txt --plot xmm_lss/p0opt_B.png
```
A previous result given by `--init` is included as one of the starting points (the result is never worse than it).

### 4.3. Optimizing K patterns together
```bash
optimize_pointings_p0 -c netflow_config_2026Nov.yaml \
    --grid xmm_lss/hexagons_xmm_lss_flat_centers.ecsv \
    --n-patterns 6 --init xmm_lss/hexagons_xmm_lss_flat_centers_p0opt.ecsv \
    --global-shift 0.6 --global-rot-step 10 --n-refine 3 --cycles 4 --processes 20 \
    -o xmm_lss/hexagons_xmm_lss_flat_centers_p0opt_6pat.ecsv \
    --add-columns --report-file xmm_lss/p0opt_6pat_report.txt --plot xmm_lss/p0opt_6pat.png
```
The output `ppc_code` gets `_A`, `_B`, …, and the `pattern` column holds the pattern name.

### 4.4. Run time
With 20 processes: a few to ~15 minutes to add one pattern, ~20 minutes for 4 patterns, and 30 minutes to 1 hour for 6 patterns. Evaluations are shared between patterns and cycles. The work is dominated by instrument model calls from Python, so a GPU does not help.

---

## 5. Command-Line Arguments

| Argument | Default | Description |
|---|---|---|
| `--config`, `-c` | (required) | Pipeline configuration YAML (targets, Gaia catalog, guide star constraints, ...) |
| `--config-yaml` | `config.yaml` | config.yaml with PFS parameters |
| `--grid`, `-i` | (required) | Reference tiling ECSV (`ppc_code`, `ppc_ra`, `ppc_dec`[, `ppc_pa`]) |
| `--output`, `-o` | (required) | Output pointing ECSV |
| `--targets` | from config | Science target CSV |
| `--gaia` | from config | Gaia catalog ECSV |
| `--center RA DEC` | field closest to the grid mean | Center of the target region |
| `--target-radius` | 1.8 | Only targets within this radius (deg) of the center are used |
| `--max-priority` | 2 | Maximum priority loaded (P0 is maximized, others break ties) |
| `--n-patterns` | 1 | Number of patterns to optimize |
| `--fixed` | none | Pattern ECSVs kept fixed and observed in addition (several allowed) |
| `--init` | none | ECSVs used as the initial layouts of the first patterns (several allowed) |
| `--cycles` | 2 | Max rounds re-optimizing each pattern in turn (several patterns) |
| `--cycle-n-refine` | 1 | Number of best global solutions refined in a cycle |
| `--global-shift` | 0.2 | Max rigid shift of the tiling (deg) |
| `--global-shift-step` | 0.1 | Shift step (deg) |
| `--global-rot-step` | 5.0 | Rigid rotation step (deg, searched over 0–60°) |
| `--n-refine` | 5 | Number of best global solutions refined locally |
| `--local-shift` | 0.08 | Max shift of an individual field (deg) |
| `--local-shift-step` | 0.02 | Local shift step (deg) |
| `--local-rot` | 5.0 | Max extra rotation of an individual field (deg) |
| `--local-rot-step` | 2.5 | Local rotation step (deg) |
| `--max-gap` | 0.03 | Max uncovered fraction of the reference tiling |
| `--max-passes` | 3 | Max coordinate descent passes in the local search |
| `--baseline` | none | Pointing ECSVs evaluated for comparison only (not used in the optimization) |
| `--code-suffix` | none | Suffix appended to the output `ppc_code` (e.g. `_B`) |
| `--processes` | min(16, CPUs) | Number of worker processes |
| `--add-columns` | off | Add `ag0`–`ag5`, `n_p0_access`, `n_p0_assigned` to the output |
| `--report-file` | none | Save the report to a text file |
| `--plot` | none | Save a plot of the optimized layout |

### Outputs
* **ECSV**: `ppc_code`, `ppc_ra`, `ppc_dec`, `ppc_pa`, `ppc_resolution`, `ppc_priority` (plus `pattern` with several patterns). It can be used directly as the `pointing_file` of `run_netflow`.
* **Report**: the score; for each field the offset, guide star counts, accessible and assigned P0 counts; the distribution of P0 visits (several patterns); the list of unassignable P0 targets.
* **Plot**: the field hexagons with `ppc_code`, P0 colored by the number of visits, and P1–P2 in gray. P0 targets outside `--target-radius` are shown as open symbols.

---

## 6. Multi-Visit Assignment in netflow (`prefer_visits`)

There are two ways to assign fibers to the K optimized patterns with netflow.

### 6.1. Running each pattern separately
Run `run_netflow` K times, each with a configuration whose `pointing_file` holds one pattern. The patterns are assigned independently, so the same target is not necessarily chosen again; lower priorities, which are short of fibers, end up with many targets observed only once.

### 6.2. Running all patterns together (recommended)
Put the fields of all patterns into one `pointing_file`, add `prefer_visits` to the configuration and solve with a single `run_netflow`:
```yaml
netflow:
  t_obs: 450.0              # exposure per field (t_obs / K to keep the total time of one pattern)
  prefer_visits: 2          # targets accessible from >= 2 fields are observed twice or not at all
  single_visit_weight: 0.5  # value of one visit for targets accessible from one field only (2 visits = 1)
```

* **Targets accessible from two or more fields** require `prefer_visits × t_obs` and are handled by netflow's native all-or-nothing rule. They are never observed only once because of a shortage of fibers, so the visits become "2 : 1 : 0 ≈ x : 0 : (1−x)".
* **Targets accessible from one field only** (hexagon edges, gaps between Cobras, broken Cobras, ...) may be observed once, with value `single_visit_weight × N` (N: the cost of its priority). The default 0.5 gives the same value per fiber as two visits.
* The priority order is preserved (the class costs differ by a factor of 3), and each target gets at most `prefer_visits` visits.
* The implementation has been tested with `prefer_visits: 2`.

The log shows, for each priority, a table of accessible fields (1 / 2+) against assigned visits (0x / 1x / 2x):
```
Science targets by class, accessible visits and assigned visits:
    sci_P0 accessible 1 : 0x:      1, 1x:     54, 2x:      0
    sci_P0 accessible 2+: 0x:      0, 1x:      0, 2x:    129
    ...
```
Without `prefer_visits`, every target is assigned as a single visit as before.

### 6.3. Example (XMM-LSS, within 1.8° of the center)

| Method | P0 (180) 2x / 1x / 0x | P3 (2052) 2x / 1x / 0x | P9 (23489) 2x / 1x / 0x |
|---|---|---|---|
| 1 pattern | – / 153 / 27 | – / 1556 / 496 | – / 8330 / 15159 |
| 2 patterns, run separately | 129 / 47 / 4 | 1161 / 725 / 166 | 3782 / 9080 / 10627 |
| 2 patterns, `prefer_visits: 2` | 129 / 47 / 4 | 1230 / 560 / 262 | 5834 / 4484 / 13171 |

One visit is the full exposure with one pattern and half of it with two patterns. The 129 P0 targets observed twice equal the number accessible from both patterns found by the optimizer.
