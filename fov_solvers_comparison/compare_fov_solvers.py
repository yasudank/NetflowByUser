#!/usr/bin/env python3
import time
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from astropy.table import Table

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from optimize_hex_fov_with_guidestars import (
    optimize_fovs_with_guidestars,
    get_vertices,
    R_HEX
)
from optimize_hex_fov_milp import optimize_fovs_milp
from optimize_hex_fov_local_search import optimize_fovs_local_search

def plot_pointings(ax, df_targets, pointings, title, color):
    """
    Helper to plot targets and hexagon pointings on a matplotlib axis.
    """
    # Plot targets
    p1 = df_targets[df_targets['priority'] == 1]
    p2 = df_targets[df_targets['priority'] == 2]
    
    ax.scatter(p2['ra'], p2['dec'], c='gray', s=5, alpha=0.3, label='Priority 2')
    ax.scatter(p1['ra'], p1['dec'], c='blue', s=10, alpha=0.5, label='Priority 1')
    
    # Plot hexagon pointings
    for idx, p in enumerate(pointings):
        p_ra = p.get('ra', p.get('ppc_ra'))
        p_dec = p.get('dec', p.get('ppc_dec'))
        p_pa = p.get('pa', p.get('ppc_pa'))
        
        v_ra, v_dec = get_vertices(p_ra, p_dec, R_HEX, p_pa)
        polygon = patches.Polygon(np.column_stack([v_ra, v_dec]), closed=True,
                                  edgecolor=color, facecolor='none', linewidth=2, linestyle='-')
        ax.add_patch(polygon)
        ax.text(p_ra, p_dec, f"#{idx+1}", color=color, fontsize=12, fontweight='bold',
                ha='center', va='center')
        
    ax.set_xlabel('RA (deg)')
    ax.set_ylabel('Dec (deg)')
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.legend(loc='upper right')
    ax.invert_xaxis()  # Sky coordinates usually have RA increasing to the left

def main():
    parser = argparse.ArgumentParser(description="Compare N-pointing FoV optimization solvers.")
    parser.add_argument("--input", default="cosmos/targets_all_20260514.csv", help="Science targets CSV")
    parser.add_argument("--gaia-catalog", default="cosmos/gaia.ecsv", help="Gaia ecsv catalog path")
    parser.add_argument("--num-fovs", type=int, default=4, help="Number of pointings to optimize")
    parser.add_argument("--max-priority", type=int, default=2, help="Max priority to consider")
    parser.add_argument("--min-stars", type=int, default=2, help="Min stars per guide camera")
    parser.add_argument("--min-cams", type=int, default=6, help="Min cameras with enough stars")
    parser.add_argument("--obstime", default="2026-05-09T06:00:00Z", help="Observation time UTC")
    
    args = parser.parse_args()
    
    # 1. Load targets
    print(f"Reading targets from: {args.input}")
    df_targets = pd.read_csv(args.input)
    df_filtered = df_targets[df_targets['priority'] <= args.max_priority].copy()
    print(f"Loaded {len(df_targets)} targets. Filtered to {len(df_filtered)} with priority <= {args.max_priority}")
    
    # 2. Load Gaia catalog
    print(f"Reading Gaia catalog from: {args.gaia_catalog}")
    t_gaia = Table.read(args.gaia_catalog, format="ascii.ecsv")
    df_gaia = t_gaia.to_pandas()
    df_gaia["magnitude"] = df_gaia["phot_g_mean_mag"]
    df_gaia["color"] = df_gaia["bp_rp"]
    if "catalog" not in df_gaia.columns:
        df_gaia["catalog"] = "gaia_dr3"
    for col in ["pmra_error", "pmdec_error", "parallax_over_error", 
                "astrometric_excess_noise", "astrometric_excess_noise_sig", 
                "ruwe", "r_cmodel_mag", "r_cmodel_magerr", 
                "r_extendedness_value", "phot_g_mean_flux_over_error", 
                "in_galaxy_candidates"]:
        if col not in df_gaia.columns:
            df_gaia[col] = np.nan
    df_gaia = df_gaia.fillna({"parallax": 1.0e-07, "pmra": 0.0, "pmdec": 0.0})
    
    results = {}
    
    # --- Run Greedy Solver ---
    print("\n" + "="*50)
    print("Running Greedy Solver...")
    print("="*50)
    t0 = time.time()
    greedy_pointings, greedy_covered_ids = optimize_fovs_with_guidestars(
        df_filtered, df_gaia, args.obstime, num_fovs=args.num_fovs,
        min_stars_per_cam=args.min_stars, min_cams_with_stars=args.min_cams
    )
    t_greedy = time.time() - t0
    
    # Greedy returns a boolean mask for covered targets
    greedy_covered_ids_set = set(df_filtered.loc[greedy_covered_ids, 'obj_id'])
    
    greedy_p1 = len(df_filtered[df_filtered['obj_id'].isin(greedy_covered_ids_set) & (df_filtered['priority'] == 1)])
    greedy_p2 = len(df_filtered[df_filtered['obj_id'].isin(greedy_covered_ids_set) & (df_filtered['priority'] == 2)])
    results['Greedy'] = {
        'time': t_greedy,
        'pointings': greedy_pointings,
        'covered_all': len(greedy_covered_ids_set),
        'covered_p1': greedy_p1,
        'covered_p2': greedy_p2
    }
    
    # --- Run Method A (MILP) ---
    print("\n" + "="*50)
    print("Running Method A (Discretization + Gurobi MILP)...")
    print("="*50)
    t0 = time.time()
    milp_pointings, milp_covered_ids = optimize_fovs_milp(
        df_filtered, df_gaia, args.obstime, num_fovs=args.num_fovs,
        min_stars_per_cam=args.min_stars, min_cams_with_stars=args.min_cams,
        grid_step=0.08, pa_step=20.0
    )
    t_milp = time.time() - t0
    
    milp_p1 = len(df_filtered[df_filtered['obj_id'].isin(milp_covered_ids) & (df_filtered['priority'] == 1)])
    milp_p2 = len(df_filtered[df_filtered['obj_id'].isin(milp_covered_ids) & (df_filtered['priority'] == 2)])
    results['Method A (MILP)'] = {
        'time': t_milp,
        'pointings': milp_pointings,
        'covered_all': len(milp_covered_ids),
        'covered_p1': milp_p1,
        'covered_p2': milp_p2
    }
    
    # --- Run Method C (Local Search) ---
    print("\n" + "="*50)
    print("Running Method C (Greedy Init + Coordinate Descent Local Search)...")
    print("="*50)
    t0 = time.time()
    ls_pointings, ls_covered_ids = optimize_fovs_local_search(
        df_filtered, df_gaia, args.obstime, num_fovs=args.num_fovs,
        min_stars_per_cam=args.min_stars, min_cams_with_stars=args.min_cams,
        initial_pointings=greedy_pointings
    )
    t_ls = time.time() - t0
    
    ls_p1 = len(df_filtered[df_filtered['obj_id'].isin(ls_covered_ids) & (df_filtered['priority'] == 1)])
    ls_p2 = len(df_filtered[df_filtered['obj_id'].isin(ls_covered_ids) & (df_filtered['priority'] == 2)])
    results['Method C (Local Search)'] = {
        'time': t_ls,
        'pointings': ls_pointings,
        'covered_all': len(ls_covered_ids),
        'covered_p1': ls_p1,
        'covered_p2': ls_p2
    }
    
    # --- Print Benchmark Summary ---
    print("\n" + "="*60)
    print("BENCHMARK COMPARISON SUMMARY")
    print("="*60)
    print(f"{'Solver':<30} | {'Time (s)':<10} | {'Unique Cover':<12} | {'Priority 1':<10} | {'Priority 2':<10}")
    print("-"*78)
    for k, v in results.items():
        print(f"{k:<30} | {v['time']:<10.2f} | {v['covered_all']:<12} | {v['covered_p1']:<10} | {v['covered_p2']:<10}")
    print("="*60)
    
    # --- Plot Comparison ---
    fig, axes = plt.subplots(1, 3, figsize=(20, 6), sharex=True, sharey=True)
    colors = {'Greedy': 'red', 'Method A (MILP)': 'green', 'Method C (Local Search)': 'purple'}
    
    for idx, (name, res) in enumerate(results.items()):
        plot_pointings(axes[idx], df_filtered, res['pointings'], f"{name} (Cover={res['covered_all']})", colors[name])
        
    plt.tight_layout()
    plot_filename = "compare_fov_coverage.png"
    plt.savefig(plot_filename, dpi=150)
    print(f"\nSaved comparison plot to: {plot_filename}")

if __name__ == "__main__":
    main()
