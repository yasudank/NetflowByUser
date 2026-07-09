#!/usr/bin/env python3
import time
import numpy as np
import pandas as pd
from astropy.table import Table
import gurobipy as gp
from gurobipy import GRB

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src", "netflow_by_user")))

from optimize_hex_fov_with_guidestars import (
    evaluate_candidates_chunk,
    evaluate_guidestars_single,
    score_pointing,
    R_HEX
)

def optimize_fovs_milp(df_targets, df_gaia, obstime, num_fovs=1,
                       min_stars_per_cam=2, min_cams_with_stars=6,
                       grid_step=0.05, pa_step=15.0, time_limit=60):
    """
    Method A: Discretization of candidates + Gurobi MILP Solver
    """
    start_time = time.time()
    
    ra = df_targets['ra'].values
    dec = df_targets['dec'].values
    priorities = df_targets['priority'].values
    N_targets = len(ra)
    
    # 1. Generate candidate centers grid
    ra_min, ra_max = ra.min(), ra.max()
    dec_min, dec_max = dec.min(), dec.max()
    
    grid_ra = np.arange(ra_min - R_HEX * 0.4, ra_max + R_HEX * 0.4, grid_step)
    grid_dec = np.arange(dec_min - R_HEX * 0.4, dec_max + R_HEX * 0.4, grid_step)
    grid_ra_mesh, grid_dec_mesh = np.meshgrid(grid_ra, grid_dec)
    
    # Filter grid centers to only keep those close to at least one target to keep candidates small
    grid_coords = np.column_stack([grid_ra_mesh.flatten(), grid_dec_mesh.flatten()])
    target_coords = np.column_stack([ra, dec])
    
    cos_dec_mid = np.cos(np.radians(dec_min + (dec_max - dec_min)/2.0))
    dx = (grid_coords[:, 0, np.newaxis] - target_coords[np.newaxis, :, 0]) * cos_dec_mid
    dy = grid_coords[:, 1, np.newaxis] - target_coords[np.newaxis, :, 1]
    dist = np.sqrt(dx**2 + dy**2)
    
    close_mask = np.any(dist < R_HEX * 1.0, axis=1)
    filtered_grid = grid_coords[close_mask]
    
    cand_ra = np.concatenate([filtered_grid[:, 0], ra])
    cand_dec = np.concatenate([filtered_grid[:, 1], dec])
    
    pas = np.arange(0.0, 60.0, pa_step)
    
    print(f"[MILP] N targets: {N_targets}")
    print(f"[MILP] Candidates grid: {len(cand_ra)} centers, {len(pas)} PAs")
    
    # 2. Pre-evaluate target coverage for all candidates to filter empty/poor regions
    print("[MILP] Pre-evaluating target coverage for all candidates...")
    candidates_with_coverage = []
    checked_set = set()
    
    for pa in pas:
        counts = evaluate_candidates_chunk(ra, dec, cand_ra, cand_dec, pa)
        for idx, count in enumerate(counts):
            if count >= 3: # Keep candidates covering at least 3 targets
                r_t = cand_ra[idx]
                d_t = cand_dec[idx]
                key = (round(r_t, 4), round(d_t, 4), round(pa, 2))
                if key not in checked_set:
                     checked_set.add(key)
                     candidates_with_coverage.append({
                         'ra': r_t,
                         'dec': d_t,
                         'pa': pa,
                         'count': count
                     })
                     
    # Sort by target count descending and keep top 800 to limit expensive guide star checks
    candidates_with_coverage = sorted(candidates_with_coverage, key=lambda x: x['count'], reverse=True)
    candidates_with_coverage = candidates_with_coverage[:800]
    print(f"[MILP] Candidates to check (top 800 with coverage >= 3): {len(candidates_with_coverage)}")
    
    valid_candidates = []
    print("[MILP] Checking guide star constraints for candidates with coverage...")
    
    for cand in candidates_with_coverage:
        r_t, d_t, pa = cand['ra'], cand['dec'], cand['pa']
        # Guide star check
        star_counts, stars_df = evaluate_guidestars_single(
            r_t, d_t, pa, df_gaia, obstime,
            min_mag=12.0, max_mag=21.5, minsep_arcsec=1.0
        )
        cams_ok = sum(1 for count in star_counts if count >= min_stars_per_cam)
        if cams_ok >= min_cams_with_stars:
            valid_candidates.append({
                'ra': r_t,
                'dec': d_t,
                'pa': pa,
                'star_counts': star_counts,
                'stars_df': stars_df
            })
                
    N_cand = len(valid_candidates)
    print(f"[MILP] Found {N_cand} valid candidates satisfying guide star constraints.")
    if N_cand == 0:
        print("[MILP] Warning: No candidate pointings met the guide star constraints!")
        return [], set()
        
    # 3. Compute coverage matrix (N_cand x N_targets)
    print("[MILP] Computing coverage matrix...")
    valid_ra = np.array([c['ra'] for c in valid_candidates])
    valid_dec = np.array([c['dec'] for c in valid_candidates])
    
    # We will compute coverage for each PA separately
    coverage = np.zeros((N_cand, N_targets), dtype=bool)
    
    for pa in pas:
        pa_indices = [idx for idx, c in enumerate(valid_candidates) if c['pa'] == pa]
        if len(pa_indices) == 0:
            continue
        
        pa_ra = valid_ra[pa_indices]
        pa_dec = valid_dec[pa_indices]
        
        # Evaluate chunks for this PA
        counts = evaluate_candidates_chunk(ra, dec, pa_ra, pa_dec, pa)
        
        # Now fill in the boolean coverage matrix
        # Let's do it efficiently
        cos_dec = np.cos(np.radians(pa_dec))
        rad = np.radians(pa)
        cos_pa = np.cos(rad)
        sin_pa = np.sin(rad)
        
        # Compute exact mask
        for i, idx in enumerate(pa_indices):
            c_ra = pa_ra[i]
            c_dec = pa_dec[i]
            c_cos = cos_dec[i]
            
            dx = (ra - c_ra) * c_cos
            dy = dec - c_dec
            
            dx_rot = dx * cos_pa + dy * sin_pa
            dy_rot = -dx * sin_pa + dy * cos_pa
            
            cond1 = np.abs(dy_rot) <= R_HEX * np.sqrt(3)/2.0
            cond2 = np.abs(dy_rot) + np.sqrt(3) * np.abs(dx_rot) <= np.sqrt(3) * R_HEX
            coverage[idx] = cond1 & cond2

    # 4. Formulate and solve ILP with Gurobi
    print("[MILP] Building Gurobi model...")
    env = gp.Env(empty=True)
    env.setParam("OutputFlag", 0)  # Silence Gurobi output
    env.start()
    
    model = gp.Model("FoV_Maximum_Coverage", env=env)
    model.setParam("TimeLimit", time_limit)
    
    # Variables
    # x[j] = 1 if candidate j is selected
    x = model.addVars(N_cand, vtype=GRB.BINARY, name="x")
    # y[i] = 1 if target i is covered
    y = model.addVars(N_targets, vtype=GRB.BINARY, name="y")
    
    # Constraints
    # 1. Select exactly N pointings
    model.addConstr(gp.quicksum(x[j] for j in range(N_cand)) == num_fovs, "NumFields")
    
    # 2. Target coverage constraints
    # y[i] <= sum_{j covering i} x[j]
    for i in range(N_targets):
        covering_candidates = np.where(coverage[:, i])[0]
        if len(covering_candidates) > 0:
            model.addConstr(y[i] <= gp.quicksum(x[j] for j in covering_candidates), f"Cov_{i}")
        else:
            model.addConstr(y[i] == 0, f"Cov_Zero_{i}")
            
    # Objective: Maximize weighted coverage (weight = 3 - priority, priority is 1 or 2)
    # Higher weight for priority 1
    weights = 3.0 - priorities
    model.setObjective(gp.quicksum(weights[i] * y[i] for i in range(N_targets)), GRB.MAXIMIZE)
    
    print("[MILP] Solving ILP...")
    model.optimize()
    
    # Retrieve selected pointings
    selected_pointings = []
    covered_targets = set()
    
    if model.Status == GRB.OPTIMAL or model.Status == GRB.TIME_LIMIT:
        for j in range(N_cand):
            if x[j].X > 0.5:
                selected_pointings.append(valid_candidates[j])
                # Find covered targets by this pointing
                cov_indices = np.where(coverage[j])[0]
                for idx in cov_indices:
                    covered_targets.add(df_targets.iloc[idx]['obj_id'])
        print(f"[MILP] Optimization complete in {time.time() - start_time:.2f}s. Covered: {len(covered_targets)} targets.")
    else:
        print(f"[MILP] Optimization failed with status: {model.Status}")
        
    return selected_pointings, covered_targets
