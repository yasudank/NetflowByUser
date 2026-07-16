---
marp: true
theme: gaia
_class: lead
paginate: true
backgroundColor: #f5f5f5
color: #333
style: |
  section {
    font-family: 'Helvetica Neue', Arial, sans-serif;
    padding: 40px;
    font-size: 26px;
  }
  h1 {
    color: #0b3c5d;
    font-size: 1.4em;
  }
  h2 {
    color: #328cc1;
    font-size: 1.1em;
    margin-top: 0px;
    margin-bottom: 10px;
  }
  p, ul, ol {
    margin-top: 8px;
    margin-bottom: 8px;
  }
  li {
    margin-top: 4px;
  }
  footer {
    font-size: 0.5em;
    color: #777;
  }
---

# PFS Netflow Observation Planning Pipeline
### Observer-Led Fiber Assignment & Configuration File Generation

An integrated pipeline for PFS (Prime Focus Spectrograph) observation planning. Automates telescope pointing optimization, fiber assignments, and the generation/verification of control files.

---

# 1. Field of View & Guide Star Optimization
## `optimize_hex_fov_with_guidestars.py`

Determines the optimal telescope pointings $(RA, Dec, PA)$ to maximize target coverage and satisfy tracking constraints.

- **Target Coverage**: Rapid in/out checks for science targets inside the hexagonal FoV.
- **Guide Star Selection (Auto Guide)**:
  - Projects Gaia stars to the 6 AG camera boundaries.
  - Retains stars with $12.0 < G < 21.5$ mag. Discards close pairs ($< 1.0\ \text{arcsec}$).
  - Rejects pointings with saturating stars ($G \le 12.0$) inside camera boundaries.
  - Ensures guide stars are secured in at least `min_cams_with_stars` cameras.
- **Broken Fiber Avoidance**:
  - Rejects pointings where a bright star ($G \le 12.0$) falls within $1.5\ \text{arcmin}$ of a broken fiber.

---

# 1.2. Optimization Search Algorithm
## Coarse-to-Fine Search

1. **Coarse Search**:
   - Generates a grid of candidate centers (spacing $0.03^\circ$) and appends target locations.
   - Position Angles (PAs) are scanned in $[0^\circ, 60^\circ)$ with a $5.0^\circ$ step.
   - Sorts candidates by target coverage and evaluates guide star/broken fiber constraints.
   - **Short-circuits** the search once a candidate satisfies all guide star constraints.
2. **Fine Search (Local Refinement)**:
   - Scans a local grid (spacing $0.002^\circ$) and fine PAs (spacing $0.5^\circ$) around the best coarse pointing.
   - Finalizes the parameters that yield the highest valid target coverage.

---

# 2. Fiber Assignment Optimization
## `netflow_solver.py`

Allocates the ~2,400 fiber positioners to targets (science targets, flux standards, and sky fibers) for each exposure. Solved as an integer linear programming (network flow) problem using Gurobi.

- **Target Priority**: Weighting based on target priority values (Priority 0 is highest).
- **Physical Reach**: Fibers can only be assigned to targets within their individual cobra patrol regions (donut-shaped patrol zones).
- **Collision Avoidance**: Rigid mathematical constraints prevent positioner tips from colliding (maintaining a distance $\ge 2.0\ \text{mm}$).
- **Calib/Sky Allocations**: Guarantees a minimum number of sky and flux standard fibers are assigned.

---

# 3. Dummy Target Avoidance
## `run_netflow.py`

Unassigned fiber positioners parked near bright stars ($G \le 12.0$) collect stray light, causing spectral contamination or detector saturation.

### Avoidance Procedure
1. **Detection**: Identifies unassigned healthy fibers within $1.5\ \text{arcmin}$ of a bright star.
2. **Polar Grid Search**:
   - Searches the cobra's patrol region for a safe configuration in the **opposite direction** of the bright star.
   - Enforces collision margins ($\ge 2.0\ \text{mm}$) against neighboring fibers.
   - Selects the position that maximizes distance to the star.
3. **Registration**: Back-projects this position to sky coordinates and appends a dummy target (`priority=4`) to force the fiber to park safely.

---

# 4. Deliverables & Validation
## `make_pfs_design.py`

Transforms the assignment results into files compatible with the telescope and instrument control systems.

- **`pfsDesign` FITS Files**:
  - Contains fiber mappings, target IDs, coordinates, and expected physical coordinates.
  - Used as the primary input for the instrument control software and the data reduction pipeline.
- **OPE Files**:
  - Text command scripts defining telescope target sequences and exposure parameters.
- **Validation**:
  - Simulates fiber movements to verify collision safety and patrol boundary compliance.
