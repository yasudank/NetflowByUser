#!/usr/bin/env python3
"""Optimize a tiled set of PFS pointings to maximize the number of P0 targets
that can actually be assigned a fiber, under the guide star constraints.

The tiling given by --grid (e.g. hexagons_xmm_lss_flat_centers.ecsv) is used as
the reference layout. The search proceeds in two stages:

1. Global stage: the whole tiling is shifted and rotated rigidly, and every
   pointing may additionally be rotated by a multiple of 60 deg. Since the PFS
   FoV is hexagonal, this keeps the tiling contiguous (no gaps).
2. Local stage: each pointing may be shifted / rotated slightly
   (coordinate descent), as long as the uncovered fraction of the reference
   footprint stays below --max-gap.

For each candidate pointing, the targets accessible by each Cobra are computed
with the same instrument model as netflow (black dots, broken Cobras, fiducials),
and the number of P0 targets that can be assigned is evaluated by maximum
bipartite matching over all pointings. Fiber collisions and calibration targets
are not taken into account, so the numbers are upper bounds.
Pointings that violate the guide star constraints (or have a bright star near a
broken fiber) contribute no targets; the number of such pointings is minimized
first.

Several patterns (complete tilings observed one after another) can be used:
--n-patterns K optimizes K patterns, and --fixed adds patterns that are kept
as they are. With K patterns in total a target can get up to K fibers. Layouts
are compared lexicographically by the numbers of P0 targets observed at least
1, 2, ..., K times, which favours covering every target first and then
spreading the extra visits evenly. For given pointings these numbers are
exact: the fiber count vectors form a polymatroid, so one balanced assignment
maximizes sum_t min(n_t, j) for all j at once, and that maximum is the
maximum matching with j copies of each target.
The patterns are first added one at a time with the previous ones fixed, then
re-optimized in turn (--cycles) until no pattern improves.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import warnings
try:
    from erfa import ErfaWarning
    warnings.filterwarnings("ignore", category=ErfaWarning)
except ImportError:
    pass
import argparse
import itertools
import multiprocessing as mp
import numpy as np
import pandas as pd
import yaml
from astropy.table import Table
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_flow

import ets_fiber_assigner.netflow as nf
from optimize_hex_fov_with_guidestars import (
    R_HEX, get_vertices, evaluate_guidestars_single, check_bright_stars_near_broken_fibers)
from optimize_hex_fov_local_search import is_valid_pointing

# Targets farther than this from the pointing center are ignored (same as netflow_solver)
MAX_TARGET_DIST_DEG = 0.8

# Shared state for worker processes (inherited via fork)
_CTX = {}


def normalize_pa(pa):
    return (pa + 180.0) % 360.0 - 180.0


def load_settings(args):
    """Collect parameters from the pipeline config and config.yaml."""
    with open(args.config, "r") as f:
        cfg = yaml.safe_load(f)
    pfs_cfg = {}
    if args.config_yaml and os.path.exists(args.config_yaml):
        with open(args.config_yaml, "r") as f:
            pfs_cfg = (yaml.safe_load(f) or {}).get("pfs", {})
    pfs_cfg.update(cfg.get("pfs", {}))

    inputs = cfg.get("inputs", {})
    netflow_cfg = cfg.get("netflow", {})
    gs_cfg = netflow_cfg.get("guidestars", {})
    catalog_dir = inputs.get("catalog_dir", ".")
    science = inputs.get("science_targets", "targets_all.csv")
    if not os.path.isabs(science):
        science = os.path.join(catalog_dir, science)

    broken_cobras_margin = pfs_cfg.get("brokenCobrasMargin")
    if broken_cobras_margin is None:
        broken_cobras_margin = netflow_cfg.get("broken_cobras_margin", 1.0)
    avoid_fiducials = pfs_cfg.get("avoidFiducials", True)
    if "fiducialsAvoidDistance" in pfs_cfg:
        avoid_fiducials = bool(pfs_cfg["fiducialsAvoidDistance"] > 0)

    return {
        "obstime": cfg["obstime"],
        "science_targets": args.targets or science,
        "gaia_catalog": args.gaia or inputs.get("gaia_catalog", "cosmos/gaia.ecsv"),
        "min_stars": netflow_cfg.get("min_stars_per_cam", 2),
        "min_cams": netflow_cfg.get("min_cams_with_stars", 6),
        "gs_mag_min": gs_cfg.get("mag_min", 12.0),
        "gs_mag_max": gs_cfg.get("mag_max", 21.5),
        "bright_star_mag_limit": netflow_cfg.get("bright_star_mag_limit", 12.0),
        "bright_star_radius_arcmin": netflow_cfg.get("bright_star_radius_arcmin", 1.5),
        "black_dot_radius_margin": pfs_cfg.get("black_dot_radius_margin", 1.65),
        "broken_cobras_margin": broken_cobras_margin,
        "avoid_fiducials": avoid_fiducials,
        "cobra_safety_margin": pfs_cfg.get("cobraSafetyMargin", pfs_cfg.get("cobra_safety_margin", 0.0)),
    }


def load_targets(path, center, radius_deg, max_priority):
    df = pd.read_csv(path)
    id_col = "ob_code" if "ob_code" in df.columns else ("obj_id" if "obj_id" in df.columns else "ID")
    cos_dec = np.cos(np.radians(center[1]))
    dist = np.hypot((df["ra"] - center[0]) * cos_dec, df["dec"] - center[1])
    df = df[(dist < radius_deg) & (df["priority"] <= max_priority)].reset_index(drop=True)
    tgt = [nf.ScienceTarget(r[id_col], r["ra"], r["dec"], r["exptime"], r["priority"], "sci")
           for _, r in df.iterrows()]
    return df, tgt


def load_gaia(path):
    df_gaia = Table.read(path, format="ascii.ecsv").to_pandas()
    df_gaia["magnitude"] = df_gaia["phot_g_mean_mag"]
    df_gaia["color"] = df_gaia["bp_rp"]
    df_gaia = df_gaia.fillna({"parallax": 1.0e-07, "pmra": 0.0, "pmdec": 0.0})
    for col in ["pmra_error", "pmdec_error", "parallax_over_error"]:
        if col not in df_gaia.columns:
            df_gaia[col] = np.nan
    return df_gaia


def evaluate_pointing(key):
    """Accessible (target, cobra) pairs and guide star status for one pointing."""
    from ics.cobraOps.TargetGroup import TargetGroup
    from ics.cobraOps.RandomTargetSelector import RandomTargetSelector

    ra, dec, pa = key
    st, bench = _CTX["settings"], _CTX["bench"]
    df_t = _CTX["df_targets"]

    # Check the constraints first; accessibility is not computed for invalid pointings
    valid = not check_bright_stars_near_broken_fibers(
        ra, dec, pa, _CTX["df_gaia"], st["obstime"], bench,
        radius_deg=st["bright_star_radius_arcmin"] / 60.0,
        max_mag=st["bright_star_mag_limit"])
    counts = [0] * 6
    if valid:
        counts, _ = evaluate_guidestars_single(ra, dec, pa, _CTX["df_gaia"], st["obstime"],
                                               min_mag=st["gs_mag_min"], max_mag=st["gs_mag_max"],
                                               minsep_arcsec=1.0)
        valid = is_valid_pointing(counts, st["min_stars"], st["min_cams"])

    near = np.where(np.hypot((df_t["ra"].values - ra) * np.cos(np.radians(dec)),
                             df_t["dec"].values - dec) < MAX_TARGET_DIST_DEG)[0]
    rows = np.zeros(0, dtype=int)
    cols = np.zeros(0, dtype=int)
    if valid and len(near) > 0:
        tel = nf.Telescope(ra, dec, pa, st["obstime"])
        pos = tel.get_fp_positions([_CTX["targets"][i] for i in near])
        tsel = RandomTargetSelector(bench, TargetGroup(np.array(pos)))
        tsel.calculateAccessibleTargets(safetyMargin=st["cobra_safety_margin"],
                                        brokenCobrasMargin=st["broken_cobras_margin"],
                                        avoidFiducials=st["avoid_fiducials"])
        acc = tsel.accessibleTargetIndices
        cidx, j = np.where(acc >= 0)
        rows = near[acc[cidx, j]].astype(np.int32)
        cols = cidx.astype(np.int32)
    return {"rows": rows, "cols": cols, "gs_counts": list(counts), "valid": bool(valid)}




class Evaluator:
    """Caches pointing evaluations and scores sets of pointings."""

    def __init__(self, pool, n_cobras, is_p0, capacity=1):
        self.pool = pool
        self.n_cobras = n_cobras
        self.is_p0 = is_p0
        # Maximum number of fibers per target (= number of patterns)
        self.capacity = capacity
        self.cache = {}
        # Pointings of the other patterns, always matched together with the trial pointings,
        # and the pattern each of them belongs to
        self.fixed = []
        self.fixed_groups = np.zeros(0, dtype=int)

    @staticmethod
    def key(ra, dec, pa):
        return (round(float(ra), 7), round(float(dec), 7), round(normalize_pa(pa), 4))

    def evaluate(self, keys):
        todo = list({k for k in keys if k not in self.cache})
        if todo:
            for k, res in zip(todo, self.pool.map(evaluate_pointing, todo, chunksize=4)):
                self.cache[k] = res
        return [self.cache[k] for k in keys]

    def set_fixed(self, keys, groups):
        self.fixed = self.evaluate(list(keys))
        self.fixed_groups = np.asarray(groups, dtype=int)

    def p0_edges(self, results):
        """(target, cobra, pattern) of P0 targets; fixed pointings come first in the cobra index.

        The trial pointings (results) form one more pattern, numbered -1.
        """
        rows, cols, grps = [np.zeros(0, dtype=int)], [np.zeros(0, dtype=int)], [np.zeros(0, dtype=int)]
        groups = np.concatenate([self.fixed_groups, np.full(len(results), -1)])
        for i, r in enumerate(self.fixed + list(results)):
            m = self.is_p0[r["rows"]]
            rows.append(r["rows"][m])
            cols.append(r["cols"][m] + i * self.n_cobras)
            grps.append(np.full(m.sum(), groups[i]))
        return np.concatenate(rows), np.concatenate(cols), np.concatenate(grps)

    @staticmethod
    def flow_network(rows, cols, grps):
        """Edges of source -> target -> (target, pattern) -> cobra -> sink.

        A target gets at most one fiber per pattern, as netflow assigns a target only once
        when its exposure time is reached. The source -> target capacities come first.
        """
        _, t_idx = np.unique(rows, return_inverse=True)
        tp_key = rows * (grps.max() + 2) + grps + 1
        tp_u, tp_idx = np.unique(tp_key, return_inverse=True)
        _, tp_t = np.unique(rows[np.unique(tp_key, return_index=True)[1]], return_inverse=True)
        _, c_idx = np.unique(cols, return_inverse=True)
        n_t, n_tp, n_c = t_idx.max() + 1, len(tp_u), c_idx.max() + 1
        o_t, o_tp, o_c = 2, 2 + n_t, 2 + n_t + n_tp
        src = np.concatenate([np.zeros(n_t, int), o_t + tp_t, o_tp + tp_idx, o_c + np.arange(n_c)])
        dst = np.concatenate([o_t + np.arange(n_t), o_tp + np.arange(n_tp), o_c + c_idx, np.ones(n_c, int)])
        return src, dst, n_t, o_c + n_c

    @staticmethod
    def n_fibers(net, j):
        """Maximum number of fibers when each target may get up to j of them."""
        src, dst, n_t, n_nodes = net
        cap = np.ones(len(src), dtype=np.int32)
        cap[:n_t] = j
        graph = csr_matrix((cap, (src, dst)), shape=(n_nodes, n_nodes))
        return int(maximum_flow(graph, 0, 1).flow_value)

    def score(self, keys, gap=0.0, bound=None):
        """Lexicographic score (-invalid, N>=1, ..., N>=capacity, other accessible, -gap).

        If bound is given and the score (without the gap) is certainly lower, None is returned
        as soon as that is known.
        """
        results = self.evaluate(keys)
        sc = [-sum(not r["valid"] for r in self.fixed + results)]
        rows, cols, grps = self.p0_edges(results)
        net = self.flow_network(rows, cols, grps) if len(rows) > 0 else None
        prev = 0
        for j in range(1, self.capacity + 1):
            if bound is not None and tuple(sc) < tuple(bound[:len(sc)]):
                return None
            # sum_t min(n_t, j) is maximized by one balanced assignment for all j at once,
            # so the increment is the number of targets observed at least j times
            n_j = self.n_fibers(net, j) if net is not None and (j == 1 or sc[-1] > 0) else prev
            sc.append(n_j - prev)
            prev = n_j
        others = set()
        for r in self.fixed + results:
            others.update(r["rows"][~self.is_p0[r["rows"]]].tolist())
        sc.append(len(others))
        if bound is not None and tuple(sc) < tuple(bound[:len(sc)]):
            return None
        return tuple(sc) + (-round(gap, 4),)

    def assignment(self, results):
        """Balanced assignment: number of fibers per P0 target and per pointing (fixed first).

        Maximum flow with a convex cost for the k-th fiber of a target; this gives the
        lexicographically optimal (most even) fiber count vector.
        """
        import networkx as nx
        rows, cols, grps = self.p0_edges(results)
        n_obs = np.zeros(len(self.is_p0), dtype=int)
        n_assigned = np.zeros(len(self.fixed) + len(results), dtype=int)
        if len(rows) == 0:
            return n_obs, n_assigned
        g = nx.DiGraph()
        for t in np.unique(rows):
            for k in range(self.capacity):
                g.add_edge("s", ("u", t, k), capacity=1, weight=k + 1)
                g.add_edge(("u", t, k), ("t", t), capacity=1, weight=0)
        for t, c, p in zip(rows, cols, grps):
            g.add_edge(("t", t), ("tp", t, p), capacity=1, weight=0)
            g.add_edge(("tp", t, p), ("c", c), capacity=1, weight=0)
        for c in np.unique(cols):
            g.add_edge(("c", c), "z", capacity=1, weight=0)
        flow = nx.max_flow_min_cost(g, "s", "z")
        for t in np.unique(rows):
            n_obs[t] = sum(flow["s"][("u", t, k)] for k in range(self.capacity))
        for c in np.unique(cols):
            if flow[("c", c)]["z"] > 0:
                n_assigned[c // self.n_cobras] += 1
        return n_obs, n_assigned


def hex_mask(x, y, cx, cy, pa_deg, r=R_HEX):
    """Inside test for a flat-topped hexagon in the local tangent plane (deg)."""
    a = np.radians(pa_deg)
    dx, dy = x - cx, y - cy
    u = dx * np.cos(a) + dy * np.sin(a)
    v = -dx * np.sin(a) + dy * np.cos(a)
    return (np.abs(v) <= r * np.sqrt(3) / 2.0) & (np.abs(v) + np.sqrt(3) * np.abs(u) <= np.sqrt(3) * r)


class GapMeter:
    """Fraction of the reference footprint not covered by the current pointings."""

    def __init__(self, center, ref_layout, step=0.01):
        self.center = center
        self.cos_dec = np.cos(np.radians(center[1]))
        xy = np.array([self.to_plane(ra, dec) for ra, dec, _ in ref_layout])
        lim = np.abs(xy).max() + R_HEX
        g = np.arange(-lim, lim + step, step)
        x, y = np.meshgrid(g, g)
        x, y = x.ravel(), y.ravel()
        inside = np.zeros(len(x), dtype=bool)
        for (px, py), (_, _, pa) in zip(xy, ref_layout):
            inside |= hex_mask(x, y, px, py, pa)
        self.x, self.y = x[inside], y[inside]

    def to_plane(self, ra, dec):
        return ((ra - self.center[0]) * self.cos_dec, dec - self.center[1])

    def __call__(self, layout):
        covered = np.zeros(len(self.x), dtype=bool)
        for ra, dec, pa in layout:
            px, py = self.to_plane(ra, dec)
            covered |= hex_mask(self.x, self.y, px, py, pa)
        return 1.0 - covered.mean()


def rigid_layout(offsets, base_pa, center, theta, dx, dy):
    """Rotate the tiling by theta (deg, same sense as PA) and shift it by (dx, dy) deg."""
    c, s = np.cos(np.radians(theta)), np.sin(np.radians(theta))
    cos_dec = np.cos(np.radians(center[1]))
    layout = []
    for (ox, oy), pa in zip(offsets, base_pa):
        x = c * ox - s * oy + dx
        y = s * ox + c * oy + dy
        layout.append((center[0] + x / cos_dec, center[1] + y, pa + theta))
    return layout


def choose_rotations(ev, layout, n_passes=2):
    """Pick a 60-deg multiple for each pointing by coordinate descent."""
    keys = [[Evaluator.key(ra, dec, pa + 60 * k) for k in range(6)] for ra, dec, pa in layout]
    # Start from the best rotation of each pointing on its own
    choice = []
    for kl in keys:
        res = ev.evaluate(kl)
        choice.append(max(range(6), key=lambda k: (res[k]["valid"],
                                                   int(ev.is_p0[res[k]["rows"]].sum()))))
    best = ev.score([keys[i][k] for i, k in enumerate(choice)])
    for _ in range(n_passes):
        improved = False
        for i in range(len(layout)):
            for k in range(6):
                trial = list(choice)
                trial[i] = k
                sc = ev.score([keys[j][kk] for j, kk in enumerate(trial)], bound=best)
                if sc is not None and sc > best:
                    best, choice, improved = sc, trial, True
        if not improved:
            break
    return choice, best


def local_refine(ev, layout, gap_meter, args):
    """Coordinate descent over small shifts / rotations of individual pointings."""
    cos_dec = np.cos(np.radians(gap_meter.center[1]))
    shifts = np.arange(-args.local_shift, args.local_shift + 1e-9, args.local_shift_step)
    rots = np.arange(-args.local_rot, args.local_rot + 1e-9, args.local_rot_step)
    layout = list(layout)
    best = ev.score([Evaluator.key(*p) for p in layout], gap_meter(layout))
    for n_pass in range(args.max_passes):
        improved = False
        # Invalid pointings first, before the gap budget is used up by the others
        valid = [r["valid"] for r in ev.evaluate([Evaluator.key(*p) for p in layout])]
        for i in sorted(range(len(layout)), key=lambda i: valid[i]):
            ra0, dec0, pa0 = layout[i]
            cands = [(ra0 + sx / cos_dec, dec0 + sy, pa0 + dr + 60 * k)
                     for sx, sy, dr, k in itertools.product(shifts, shifts, rots, range(6))]
            ev.evaluate([Evaluator.key(*c) for c in cands])
            for c in cands:
                trial = layout[:i] + [c] + layout[i + 1:]
                # The gap only breaks ties, so skip its (slower) computation when possible
                sc = ev.score([Evaluator.key(*p) for p in trial], bound=best[:-1])
                if sc is None:
                    continue
                gap = gap_meter(trial)
                if gap > args.max_gap:
                    continue
                sc = sc[:-1] + (-round(gap, 4),)
                if sc > best:
                    best, layout, improved = sc, trial, True
            print(f"    pass {n_pass + 1}, pointing {i}: score={best}", flush=True)
        if not improved:
            break
    return layout, best


def optimize_pattern(ev, geo, args, n_refine, init=None):
    """Global + local search of one pattern; ev.fixed holds the other patterns.

    init (a pattern dict) is refined together with the best global candidates.
    """
    offsets, grid_pa, center = geo["offsets"], geo["grid_pa"], geo["center"]
    # 1. Global stage: rigid shift / rotation of the whole tiling + 60-deg multiples
    thetas = np.arange(0.0, 60.0, args.global_rot_step)
    shifts = np.arange(-args.global_shift, args.global_shift + 1e-9, args.global_shift_step)
    cands = [(th, dx, dy) for th in thetas for dx in shifts for dy in shifts]
    print(f"\nGlobal stage: {len(cands)} rigid layouts x {len(offsets)} pointings x 6 rotations", flush=True)
    layouts = [rigid_layout(offsets, grid_pa, center, *c) for c in cands]
    ev.evaluate([Evaluator.key(ra, dec, pa + 60 * k) for lay in layouts for ra, dec, pa in lay for k in range(6)])
    global_res = []
    for c, lay in zip(cands, layouts):
        choice, sc = choose_rotations(ev, lay)
        lay = [(ra, dec, pa + 60 * k) for (ra, dec, pa), k in zip(lay, choice)]
        global_res.append((sc, c, lay))
    global_res.sort(key=lambda x: x[0], reverse=True)
    starts = [(c, lay, rigid_layout(offsets, grid_pa, center, *c)) for _, c, lay in global_res[:n_refine]]
    for sc, c, _ in global_res[:n_refine]:
        print(f"  theta={c[0]:5.1f} dx={c[1]:+.2f} dy={c[2]:+.2f} -> score={sc}")
    if init is not None:
        starts.append((init["c"], init["layout"], init["ref"]))

    # 2. Local stage
    best = None
    for c, lay, ref in starts:
        print(f"\nLocal stage from " + (f"theta={c[0]:.1f}, dx={c[1]:+.2f}, dy={c[2]:+.2f}" if c else "current layout"),
              flush=True)
        lay2, sc2 = local_refine(ev, lay, GapMeter(center, ref), args)
        if best is None or sc2 > best["score"]:
            best = {"score": sc2, "c": c, "layout": lay2, "ref": ref}
    return best


def format_score(sc, max_priority, gap=True):
    counts = ", ".join(str(n) for n in sc[1:-2])
    return (f"invalid pointings={-sc[0]}, P0 observed >=1..{len(sc) - 3} times=[{counts}], "
            f"P1-P{max_priority} accessible={sc[-2]}" + (f", gap={-sc[-1]:.4f}" if gap else ""))


def describe(ev, patterns, grid_layout, center):
    """Per-pointing summary of all patterns, with the balanced fiber assignment."""
    ev.set_fixed(pattern_keys(patterns), pattern_groups(patterns))
    results = ev.fixed
    n_obs, n_assigned = ev.assignment([])
    cos_dec = np.cos(np.radians(center[1]))
    rows = []
    i = 0
    for pat in patterns:
        for k, (ra, dec, pa) in enumerate(pat["layout"]):
            r = results[i]
            g_ra, g_dec, g_pa = grid_layout[k]
            rows.append({
                "pattern": pat["name"], "fixed": pat.get("fixed", False),
                "ppc_code": pat["codes"][k], "ra": ra, "dec": dec, "pa": normalize_pa(pa),
                "shift_arcmin": 60 * np.hypot((ra - g_ra) * cos_dec, dec - g_dec),
                "dpa": normalize_pa(pa - g_pa),
                "gs_counts": r["gs_counts"], "valid": r["valid"],
                "n_p0_access": len(np.unique(r["rows"][ev.is_p0[r["rows"]]])),
                "n_p0_assigned": int(n_assigned[i]),
            })
            i += 1
    return rows, n_obs


def observation_summary(n_obs, is_p0, capacity):
    n = n_obs[is_p0]
    hist = ", ".join(f"{j}x: {int((n == j).sum())}" for j in range(capacity, -1, -1))
    return f"P0 observation counts: {hist}; mean={n.mean():.3f}, variance={n.var():.3f}"


def write_report(rows, n_obs, df_targets, ev, args, header_lines):
    lines = list(header_lines)
    lines.append(f"{'Pointing Code':<24} {'RA':>10} {'Dec':>10} {'PA':>8} {'shift[am]':>9} {'dPA':>7} "
                 f"{'AG0..AG5':<26} {'GS ok':>5} {'P0 acc':>6} {'P0 asg':>6}")
    for r in rows:
        ag = ",".join(f"{c:d}" for c in r["gs_counts"])
        code = r["ppc_code"] + (" (fixed)" if r["fixed"] else "")
        lines.append(f"{code:<24} {r['ra']:10.5f} {r['dec']:10.5f} {r['pa']:8.2f} "
                     f"{r['shift_arcmin']:9.2f} {r['dpa']:7.1f} {ag:<26} {str(r['valid']):>5} "
                     f"{r['n_p0_access']:6d} {r['n_p0_assigned']:6d}")
    n_p0 = int(ev.is_p0.sum())
    lines.append(f"P0 targets within {args.target_radius} deg: {n_p0}, assignable (upper bound): "
                 f"{int((n_obs[ev.is_p0] > 0).sum())}")
    if ev.capacity > 1:
        lines.append(observation_summary(n_obs, ev.is_p0, ev.capacity))
    unassigned = df_targets[ev.is_p0 & (n_obs == 0)]
    id_col = "ob_code" if "ob_code" in df_targets.columns else df_targets.columns[0]
    lines.append("Unassignable P0: " + ", ".join(str(x) for x in unassigned[id_col]))
    text = "\n".join(lines)
    print(text)
    if args.report_file:
        with open(args.report_file, "w") as f:
            f.write(text + "\n")
        print(f"Saved report to {args.report_file}")


def plot_layout(path, patterns, df_targets, ev, n_obs, center, title, ref_layout=(), df_outer=None, radius=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(11, 8))
    for ra, dec, pa in ref_layout:
        v_ra, v_dec = get_vertices(ra, dec, R_HEX, pa)
        ax.plot(v_ra, v_dec, color="gray", ls=":", lw=1)
    colors = ["C0", "C4", "C1", "C5", "C6", "C9", "C8", "C7"]
    for n, pat in enumerate(patterns):
        color = colors[n % len(colors)]
        for (ra, dec, pa), code in zip(pat["layout"], pat["codes"]):
            v_ra, v_dec = get_vertices(ra, dec, R_HEX, pa)
            ax.plot(v_ra, v_dec, color=color, lw=1.5 if len(patterns) <= 2 else 1.0)
            ax.text(ra, dec, code, ha="center", va="center", fontsize=7 if len(patterns) <= 2 else 5,
                    color=color, zorder=5, bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.7))
    max_pri = df_targets["priority"].max()
    if df_outer is not None and len(df_outer) > 0:
        out_low = df_outer[df_outer["priority"] != 0]
        out_p0 = df_outer[df_outer["priority"] == 0]
        ax.scatter(out_low["ra"], out_low["dec"], s=4, c="lightgray")
        ax.scatter(out_p0["ra"], out_p0["dec"], s=14, facecolors="none", edgecolors="C7",
                   label=f"P0 outside r={radius} deg, not optimized ({len(out_p0)})")
    if radius is not None:
        t = np.linspace(0, 2 * np.pi, 361)
        ax.plot(center[0] + radius * np.cos(t) / np.cos(np.radians(center[1])), center[1] + radius * np.sin(t),
                color="gray", ls="--", lw=0.8)
    low = df_targets[~ev.is_p0]
    ax.scatter(low["ra"], low["dec"], s=4, c="lightgray", label="P1-P%d" % max_pri)
    p0 = df_targets[ev.is_p0]
    n = n_obs[ev.is_p0]
    if ev.capacity == 1:
        ax.scatter(p0["ra"][n > 0], p0["dec"][n > 0], s=14, c="C2", label=f"P0 assignable ({(n > 0).sum()})")
    else:
        for j in range(ev.capacity, 0, -1):
            sel = n == j
            ax.scatter(p0["ra"][sel], p0["dec"][sel], s=14, color=plt.cm.viridis_r(j / ev.capacity),
                       label=f"P0 assignable {j}x ({sel.sum()})")
    ax.scatter(p0["ra"][n == 0], p0["dec"][n == 0], s=20, c="C3", marker="x",
               label=f"P0 not assignable ({(n == 0).sum()})")
    ax.invert_xaxis()
    ax.set_aspect(1.0 / np.cos(np.radians(center[1])))
    ax.set_xlabel("RA [deg]")
    ax.set_ylabel("Dec [deg]")
    ax.set_title(title)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize="small")
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved plot to {path}")


def read_layout(path):
    t = Table.read(path, format="ascii.ecsv")
    pa = t["ppc_pa"] if "ppc_pa" in t.colnames else np.zeros(len(t))
    codes = [str(c) for c in t["ppc_code"]] if "ppc_code" in t.colnames else \
        [f"POINTING_{i + 1}" for i in range(len(t))]
    return [(float(r), float(d), float(p)) for r, d, p in zip(t["ppc_ra"], t["ppc_dec"], pa)], codes


def pattern_keys(patterns):
    return [Evaluator.key(*p) for pat in patterns for p in pat["layout"]]


def pattern_groups(patterns):
    return [n for n, pat in enumerate(patterns) for _ in pat["layout"]]


def main():
    parser = argparse.ArgumentParser(
        description="Optimize pointings to maximize assignable P0 targets under guide star constraints.")
    parser.add_argument("--config", "-c", required=True, help="Pipeline configuration YAML")
    parser.add_argument("--config-yaml", default="config.yaml", help="config.yaml with PFS parameters")
    parser.add_argument("--grid", "-i", required=True, help="Reference tiling ECSV (ppc_code, ppc_ra, ppc_dec[, ppc_pa])")
    parser.add_argument("--output", "-o", required=True, help="Output pointing ECSV")
    parser.add_argument("--targets", default=None, help="Science target CSV (default: from config)")
    parser.add_argument("--gaia", default=None, help="Gaia catalog ECSV (default: from config)")
    parser.add_argument("--center", type=float, nargs=2, default=None, metavar=("RA", "DEC"),
                        help="Center of the target region (default: the grid pointing closest to the grid mean)")
    parser.add_argument("--target-radius", type=float, default=1.8, help="Only targets within this radius (deg) of the center are used")
    parser.add_argument("--max-priority", type=int, default=2, help="Targets up to this priority are loaded (P0 is maximized, others break ties)")
    parser.add_argument("--n-patterns", type=int, default=1, help="Number of patterns (tilings) to optimize")
    parser.add_argument("--fixed", nargs="*", default=[],
                        help="Pointing ECSV files of patterns kept fixed and observed in addition")
    parser.add_argument("--init", nargs="*", default=[],
                        help="Pointing ECSV files used as the initial layouts of the first patterns "
                             "(skips their initial search; they are still re-optimized in the cycles)")
    parser.add_argument("--cycles", type=int, default=2,
                        help="Max rounds re-optimizing each pattern with the others fixed (multi-pattern only)")
    parser.add_argument("--cycle-n-refine", type=int, default=1,
                        help="Number of best global solutions refined in a cycle (besides the current layout)")
    parser.add_argument("--global-shift", type=float, default=0.2, help="Max rigid shift of the tiling (deg)")
    parser.add_argument("--global-shift-step", type=float, default=0.1, help="Rigid shift step (deg)")
    parser.add_argument("--global-rot-step", type=float, default=5.0, help="Rigid rotation step (deg, searched over [0, 60))")
    parser.add_argument("--n-refine", type=int, default=5, help="Number of best global solutions to refine locally")
    parser.add_argument("--local-shift", type=float, default=0.08, help="Max shift of an individual pointing (deg)")
    parser.add_argument("--local-shift-step", type=float, default=0.02, help="Local shift step (deg)")
    parser.add_argument("--local-rot", type=float, default=5.0, help="Max extra rotation of an individual pointing (deg)")
    parser.add_argument("--local-rot-step", type=float, default=2.5, help="Local rotation step (deg)")
    parser.add_argument("--max-gap", type=float, default=0.03, help="Max uncovered fraction of the reference footprint")
    parser.add_argument("--max-passes", type=int, default=3, help="Max coordinate descent passes in the local stage")
    parser.add_argument("--baseline", nargs="*", default=[],
                        help="Pointing ECSV files evaluated for comparison (each as one more pattern besides --fixed)")
    parser.add_argument("--code-suffix", default="", help="Suffix appended to the output ppc_code (e.g. _B)")
    parser.add_argument("--processes", type=int, default=min(16, os.cpu_count() or 1), help="Number of worker processes")
    parser.add_argument("--add-columns", action="store_true", help="Add ag0..ag5 and P0 count columns to the output")
    parser.add_argument("--report-file", default=None, help="Save the report to a text file")
    parser.add_argument("--plot", default=None, help="Save a plot of the optimized layout")
    args = parser.parse_args()
    if len(args.init) > args.n_patterns:
        parser.error("more --init files than --n-patterns")

    settings = load_settings(args)
    t_grid = Table.read(args.grid, format="ascii.ecsv")
    grid_layout, codes = read_layout(args.grid)
    grid_ra = np.array([p[0] for p in grid_layout])
    grid_dec = np.array([p[1] for p in grid_layout])
    grid_pa = np.array([p[2] for p in grid_layout])
    if args.center is None:
        mean_ra, mean_dec = grid_ra.mean(), grid_dec.mean()
        i_c = np.argmin(np.hypot((grid_ra - mean_ra) * np.cos(np.radians(mean_dec)), grid_dec - mean_dec))
        center = (grid_ra[i_c], grid_dec[i_c])
    else:
        center = tuple(args.center)
    cos_dec = np.cos(np.radians(center[1]))
    offsets = [((ra - center[0]) * cos_dec, dec - center[1]) for ra, dec in zip(grid_ra, grid_dec)]
    geo = {"offsets": offsets, "grid_pa": grid_pa, "center": center}

    # Output ppc_code of each new pattern: grid code + suffix (+ _A, _B, ... with several patterns)
    def pattern_codes(n):
        tag = f"_{chr(ord('A') + n)}" if args.n_patterns > 1 else ""
        return [c + args.code_suffix + tag for c in codes]

    fixed = []
    for path in args.fixed:
        lay, fcodes = read_layout(path)
        fixed.append({"name": path, "layout": lay, "codes": fcodes, "fixed": True})
    all_codes = [c for pat in fixed for c in pat["codes"]] + \
        [c for n in range(args.n_patterns) for c in pattern_codes(n)]
    if len(set(all_codes)) < len(all_codes):
        parser.error("duplicated ppc_code between --fixed and the output; set --code-suffix")
    capacity = len(fixed) + args.n_patterns

    print(f"Center: RA={center[0]:.5f}, Dec={center[1]:.5f}, target radius={args.target_radius} deg")
    df_targets, targets = load_targets(settings["science_targets"], center, args.target_radius, args.max_priority)
    is_p0 = (df_targets["priority"].values == 0)
    print(f"Targets loaded: {len(df_targets)} (P0: {is_p0.sum()}), patterns: {capacity} "
          f"({len(fixed)} fixed, {args.n_patterns} optimized)")
    df_gaia = load_gaia(settings["gaia_catalog"])

    import netflow_instrument
    bench = netflow_instrument.getBench(settings["black_dot_radius_margin"])
    _CTX.update(settings=settings, bench=bench, df_targets=df_targets, targets=targets, df_gaia=df_gaia)

    with mp.get_context("fork").Pool(args.processes) as pool:
        ev = Evaluator(pool, bench.cobras.nCobras, is_p0, capacity)
        header = []
        ev.set_fixed(pattern_keys(fixed), pattern_groups(fixed))
        if fixed:
            header.append(f"Fixed {', '.join(args.fixed)}: " + format_score(ev.score([]), args.max_priority))
            print(header[-1], flush=True)
        for path in [args.grid] + list(args.baseline):
            lay, _ = read_layout(path)
            header.append(f"Baseline {path}: " +
                          format_score(ev.score([Evaluator.key(*p) for p in lay]), args.max_priority))
            print(header[-1], flush=True)

        # Add the patterns one at a time, the previous ones fixed
        patterns = []
        for n in range(args.n_patterns):
            print(f"\n=== Pattern {n + 1}/{args.n_patterns} ===", flush=True)
            ev.set_fixed(pattern_keys(fixed + patterns), pattern_groups(fixed + patterns))
            if n < len(args.init):
                lay, _ = read_layout(args.init[n])
                gap = GapMeter(center, lay)(lay)
                pat = {"score": ev.score([Evaluator.key(*p) for p in lay], gap), "c": None, "layout": lay, "ref": lay}
                print(f"Initial layout from {args.init[n]}")
            else:
                pat = optimize_pattern(ev, geo, args, args.n_refine)
            pat.update(name=chr(ord("A") + n), codes=pattern_codes(n))
            patterns.append(pat)
            print(f"Pattern {pat['name']}: " + format_score(pat["score"], args.max_priority), flush=True)

        # Re-optimize each pattern with all the others fixed
        for cycle in range(args.cycles if len(fixed) + len(patterns) > 1 else 0):
            improved = False
            for n, pat in enumerate(patterns):
                print(f"\n=== Cycle {cycle + 1}, pattern {pat['name']} ===", flush=True)
                others = fixed + patterns[:n] + patterns[n + 1:]
                ev.set_fixed(pattern_keys(others), pattern_groups(others))
                current = ev.score([Evaluator.key(*p) for p in pat["layout"]],
                                   GapMeter(center, pat["ref"])(pat["layout"]))
                new = optimize_pattern(ev, geo, args, args.cycle_n_refine, init=pat)
                if new["score"] > current:
                    pat.update(score=new["score"], c=new["c"], layout=new["layout"], ref=new["ref"])
                    improved = True
                print(f"Pattern {pat['name']}: " + format_score(max(new["score"], current), args.max_priority),
                      flush=True)
            if not improved:
                break

        ev.set_fixed(pattern_keys(fixed + patterns), pattern_groups(fixed + patterns))
        total = ev.score([])
        for pat in patterns:
            c = pat["c"]
            header.append(f"Optimized pattern {pat['name']}: " +
                          (f"rigid theta={c[0]:.1f} deg, shift=({c[1]:+.3f}, {c[2]:+.3f}) deg"
                           if c else "from the initial layout"))
        header.append("Total: " + format_score(total, args.max_priority, gap=False))
        rows, n_obs = describe(ev, fixed + patterns, grid_layout, center)

    print()
    write_report(rows, n_obs, df_targets, ev, args, header)

    rows_out = [r for r in rows if not r["fixed"]]
    t_out = Table()
    t_out["ppc_code"] = [r["ppc_code"] for r in rows_out]
    t_out["ppc_ra"] = [r["ra"] for r in rows_out]
    t_out["ppc_dec"] = [r["dec"] for r in rows_out]
    t_out["ppc_pa"] = [r["pa"] for r in rows_out]
    for col, default in [("ppc_resolution", "L"), ("ppc_priority", 0.0)]:
        values = list(t_grid[col]) if col in t_grid.colnames else [default] * len(t_grid)
        t_out[col] = values * args.n_patterns
    if args.n_patterns > 1:
        t_out["pattern"] = [r["pattern"] for r in rows_out]
    if args.add_columns:
        for cam_idx in range(6):
            t_out[f"ag{cam_idx}"] = [r["gs_counts"][cam_idx] for r in rows_out]
        t_out["n_p0_access"] = [r["n_p0_access"] for r in rows_out]
        t_out["n_p0_assigned"] = [r["n_p0_assigned"] for r in rows_out]
    t_out.meta["comments"] = [f"Optimized by optimize_pointings_p0 from {args.grid}"] + header
    t_out.write(args.output, format="ascii.ecsv", overwrite=True)
    print(f"Saved optimized pointings to {args.output}")

    if args.plot:
        df_all = pd.read_csv(settings["science_targets"])
        dist = np.hypot((df_all["ra"] - center[0]) * cos_dec, df_all["dec"] - center[1])
        df_outer = df_all[(dist >= args.target_radius) & (df_all["priority"] <= args.max_priority)]
        n_p0 = n_obs[is_p0]
        title = f"P0 assignable: {int((n_p0 > 0).sum())}/{int(is_p0.sum())} (within r={args.target_radius} deg)"
        if capacity > 1:
            title += f", {capacity} patterns, mean {n_p0.mean():.2f}"
        plot_layout(args.plot, fixed + patterns, df_targets, ev, n_obs, center, title,
                    ref_layout=patterns[0]["ref"] if len(patterns) == 1 and not fixed else (),
                    df_outer=df_outer, radius=args.target_radius)


if __name__ == "__main__":
    main()
