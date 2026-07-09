#!/usr/bin/env python3
import time
import numpy as np
import pandas as pd
from astropy.table import Table

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from optimize_hex_fov_with_guidestars import (
    optimize_fovs_with_guidestars,
    evaluate_guidestars_single,
    R_HEX
)

def get_pointing_coverage_mask(pointing, ra, dec):
    """
    Compute target coverage mask for a single pointing.
    """
    c_ra = pointing.get('ra', pointing.get('ppc_ra'))
    c_dec = pointing.get('dec', pointing.get('ppc_dec'))
    pa = pointing.get('pa', pointing.get('ppc_pa'))
    
    cos_dec = np.cos(np.radians(c_dec))
    rad = np.radians(pa)
    cos_pa = np.cos(rad)
    sin_pa = np.sin(rad)
    
    dx = (ra - c_ra) * cos_dec
    dy = dec - c_dec
    
    dx_rot = dx * cos_pa + dy * sin_pa
    dy_rot = -dx * sin_pa + dy * cos_pa
    
    cond1 = np.abs(dy_rot) <= R_HEX * np.sqrt(3)/2.0
    cond2 = np.abs(dy_rot) + np.sqrt(3) * np.abs(dx_rot) <= np.sqrt(3) * R_HEX
    return cond1 & cond2

def evaluate_total_score(pointings, ra, dec, priorities):
    """
    Calculate the total weighted coverage score for a set of pointings.
    """
    if len(pointings) == 0:
        return 0.0, set()
        
    combined_mask = np.zeros(len(ra), dtype=bool)
    for p in pointings:
        mask = get_pointing_coverage_mask(p, ra, dec)
        combined_mask |= mask
        
    weights = 3.0 - priorities
    score = np.sum(weights[combined_mask])
    
    # Return score and the set of covered target indices
    covered_indices = np.where(combined_mask)[0]
    return score, covered_indices

def optimize_fovs_local_search(df_targets, df_gaia, obstime, num_fovs=1,
                              min_stars_per_cam=2, min_cams_with_stars=6,
                              initial_pointings=None):
    """
    Method C: Greedy Initialization + Hill-Climbing Coordinate Descent Local Search
    """
    start_time = time.time()
    
    ra = df_targets['ra'].values
    dec = df_targets['dec'].values
    priorities = df_targets['priority'].values
    obj_ids = df_targets['obj_id'].values
    
    # 1. Initialize pointings using Greedy solver if not provided
    if initial_pointings is None:
        print("[Local Search] Running Greedy solver for initial pointings...")
        initial_pointings, _ = optimize_fovs_with_guidestars(
            df_targets, df_gaia, obstime, num_fovs=num_fovs,
            min_stars_per_cam=min_stars_per_cam, min_cams_with_stars=min_cams_with_stars
        )
        
    current_pointings = [dict(p) for p in initial_pointings]
    current_score, covered_indices = evaluate_total_score(current_pointings, ra, dec, priorities)
    
    print(f"[Local Search] Initial Greedy score: {current_score:.1f} (Covered unique targets: {len(covered_indices)})")
    
    # Define step sizes for coordinate search
    step_sizes = [0.02, 0.005, 0.001]
    pa_step_sizes = [5.0, 1.0, 0.2]
    
    improved = True
    iterations = 0
    
    for step, pa_step in zip(step_sizes, pa_step_sizes):
        print(f"[Local Search] Fine-tuning with step size={step:.4f}°, PA step={pa_step:.1f}°")
        step_improved = True
        
        while step_improved:
            step_improved = False
            
            # Loop over all pointings
            for k in range(len(current_pointings)):
                p = current_pointings[k]
                best_p = dict(p)
                best_score = current_score
                
                # Test local perturbations in RA, Dec, and PA
                # 3^3 = 27 candidates, but let's check RA/Dec and PA perturbations
                ra_shifts = [-step, 0.0, step]
                dec_shifts = [-step, 0.0, step]
                pa_shifts = [-pa_step, 0.0, pa_step]
                
                for dr in ra_shifts:
                    for dd in dec_shifts:
                        for dp in pa_shifts:
                            if dr == 0.0 and dd == 0.0 and dp == 0.0:
                                continue
                                
                            p_ra = p.get('ra', p.get('ppc_ra'))
                            p_dec = p.get('dec', p.get('ppc_dec'))
                            p_pa = p.get('pa', p.get('ppc_pa'))
                            
                            test_ra = p_ra + dr
                            test_dec = p_dec + dd
                            test_pa = (p_pa + dp) % 60.0
                            
                            # Verify guide stars
                            star_counts, stars_df = evaluate_guidestars_single(
                                test_ra, test_dec, test_pa, df_gaia, obstime,
                                min_mag=12.0, max_mag=19.0, minsep_arcsec=1.0
                            )
                            cams_ok = sum(1 for count in star_counts if count >= min_stars_per_cam)
                            if cams_ok < min_cams_with_stars:
                                continue  # Violates guide star constraints
                                
                            # Evaluate total score
                            test_pointings = list(current_pointings)
                            test_pointings[k] = {
                                'ra': test_ra,
                                'dec': test_dec,
                                'pa': test_pa,
                                'ppc_ra': test_ra,
                                'ppc_dec': test_dec,
                                'ppc_pa': test_pa,
                                'star_counts': star_counts,
                                'stars_df': stars_df
                            }
                            
                            score, test_cov = evaluate_total_score(test_pointings, ra, dec, priorities)
                            if score > best_score:
                                best_score = score
                                best_p = test_pointings[k]
                                step_improved = True
                                
                if step_improved:
                    current_pointings[k] = best_p
                    current_score = best_score
                    print(f"  [Local Search] Pointing {k+1} shifted. New score: {current_score:.1f}")
                    
            iterations += 1
            if iterations > 100:  # Safety breakout
                break
                
    # Re-evaluate final covered target IDs
    final_score, final_cov_indices = evaluate_total_score(current_pointings, ra, dec, priorities)
    covered_target_ids = set(obj_ids[final_cov_indices])
    
    print(f"[Local Search] Local Search complete in {time.time() - start_time:.2f}s. Final score: {final_score:.1f}, unique targets covered: {len(covered_target_ids)}.")
    return current_pointings, covered_target_ids
