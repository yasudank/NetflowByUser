#!/usr/bin/env python3
import os
import math
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from astropy.table import Table

# Regular hexagon parameters
R_HEX = 0.69  # Outer radius of PFS FoV in degrees

def is_inside_hexagon_single(ra, dec, center_ra, center_dec, r, pa_deg):
    """Check if a single target is inside a rotated hexagon (for verification/plotting)"""
    cos_dec = math.cos(math.radians(center_dec))
    dx = (ra - center_ra) * cos_dec
    dy = dec - center_dec
    
    # Rotate targets by -pa_deg (clockwise)
    rad = np.radians(pa_deg)
    cos_pa = np.cos(rad)
    sin_pa = np.sin(rad)
    
    dx_rot = dx * cos_pa + dy * sin_pa
    dy_rot = -dx * sin_pa + dy * cos_pa
    
    # Check against flat-topped hexagon (PA=0)
    return (abs(dy_rot) <= r * np.sqrt(3)/2.0) and (abs(dy_rot) + np.sqrt(3) * abs(dx_rot) <= np.sqrt(3) * r)

def get_vertices(ra, dec, r, pa_deg):
    """Get the 6 vertices of a rotated hexagon centered at (ra, dec) with radius r and rotation pa_deg"""
    cos_dec = math.cos(math.radians(dec))
    
    # Vertices of standard flat-topped hexagon (PA=0)
    angles = np.radians([60 * i for i in range(7)])
    v_dx = r * np.cos(angles)
    v_dy = r * np.sin(angles)
    
    # Rotate vertices by pa_deg (counter-clockwise)
    rad = np.radians(pa_deg)
    cos_pa = np.cos(rad)
    sin_pa = np.sin(rad)
    
    v_dx_rot = v_dx * cos_pa - v_dy * sin_pa
    v_dy_rot = v_dx * sin_pa + v_dy * cos_pa
    
    v_ra = ra + v_dx_rot / cos_dec
    v_dec = dec + v_dy_rot
    return v_ra, v_dec

def evaluate_candidates_chunk(ra, dec, cand_ra, cand_dec, pa_deg, chunk_size=5000):
    """
    Evaluate coverage for candidates in chunks to prevent memory blowup.
    Supports arbitrary rotation angle pa_deg.
    Returns:
        best_idx: Index of the best candidate in chunk
        best_count: Max count of targets covered
        best_mask: Boolean mask of targets covered by the best candidate
    """
    K = len(cand_ra)
    N = len(ra)
    
    covered_counts = np.zeros(K, dtype=int)
    
    # cos(dec) for candidates
    cos_dec = np.cos(np.radians(cand_dec))
    
    best_count = -1
    best_idx = -1
    best_mask = None
    
    rad = np.radians(pa_deg)
    cos_pa = np.cos(rad)
    sin_pa = np.sin(rad)
    
    for start_k in range(0, K, chunk_size):
        end_k = min(start_k + chunk_size, K)
        c_ra_chunk = cand_ra[start_k:end_k]
        c_dec_chunk = cand_dec[start_k:end_k]
        cos_dec_chunk = cos_dec[start_k:end_k][:, np.newaxis]
        
        # dx, dy shape: (chunk_K, N)
        dx = (ra[np.newaxis, :] - c_ra_chunk[:, np.newaxis]) * cos_dec_chunk
        dy = dec[np.newaxis, :] - c_dec_chunk[:, np.newaxis]
        
        # Rotate targets by -pa_deg (clockwise)
        dx_rot = dx * cos_pa + dy * sin_pa
        dy_rot = -dx * sin_pa + dy * cos_pa
        
        # Check against flat-topped hexagon (PA=0)
        cond1 = np.abs(dy_rot) <= R_HEX * np.sqrt(3)/2.0
        cond2 = np.abs(dy_rot) + np.sqrt(3) * np.abs(dx_rot) <= np.sqrt(3) * R_HEX
        mask_chunk = cond1 & cond2
            
        counts = mask_chunk.sum(axis=1)
        covered_counts[start_k:end_k] = counts
        
        max_idx_in_chunk = np.argmax(counts)
        if counts[max_idx_in_chunk] > best_count:
            best_count = counts[max_idx_in_chunk]
            best_idx = start_k + max_idx_in_chunk
            best_mask = mask_chunk[max_idx_in_chunk]
            
    return best_idx, best_count, best_mask

def optimize_fovs(df_targets, num_fovs=1, coarse_grid_step=0.03, fine_grid_step=0.002, pa_step=5.0):
    """
    Greedy maximum coverage solver with coarse-to-fine search and arbitrary PA optimization
    """
    ra = df_targets['ra'].values
    dec = df_targets['dec'].values
    ids = df_targets['obj_id'].values
    
    N = len(ra)
    covered = np.zeros(N, dtype=bool)
    
    pointings = []
    
    # 1. Generate coarse candidate grid over the target region
    ra_min, ra_max = ra.min(), ra.max()
    dec_min, dec_max = dec.min(), dec.max()
    
    # Add buffer of R_HEX to make sure we cover edge targets
    grid_ra = np.arange(ra_min - R_HEX * 0.5, ra_max + R_HEX * 0.5, coarse_grid_step)
    grid_dec = np.arange(dec_min - R_HEX * 0.5, dec_max + R_HEX * 0.5, coarse_grid_step)
    
    grid_ra_mesh, grid_dec_mesh = np.meshgrid(grid_ra, grid_dec)
    coarse_cand_ra = grid_ra_mesh.flatten()
    coarse_cand_dec = grid_dec_mesh.flatten()
    
    # Also add the coordinates of the targets themselves as candidates
    cand_ra = np.concatenate([coarse_cand_ra, ra])
    cand_dec = np.concatenate([coarse_cand_dec, dec])
    
    # Generate candidate PAs (hexagon is 60-degree symmetric)
    coarse_pas = np.arange(0.0, 60.0, pa_step)
    
    print(f"Number of targets: {N}")
    print(f"Number of coarse candidate centers: {len(cand_ra)}")
    print(f"Candidate PAs to evaluate: {coarse_pas}")
    
    for fov_idx in range(num_fovs):
        print(f"\nOptimizing FoV {fov_idx+1}/{num_fovs}...")
        
        # Only evaluate against currently uncovered targets
        active_ra = ra[~covered]
        active_dec = dec[~covered]
        
        if len(active_ra) == 0:
            print("All targets covered!")
            break
            
        best_count = -1
        best_pa = -1.0
        best_center_ra = -1.0
        best_center_dec = -1.0
        best_mask = None
        
        # Coarse search over grid and PAs
        for pa in coarse_pas:
            idx, count, mask = evaluate_candidates_chunk(
                active_ra, active_dec, cand_ra, cand_dec, pa
            )
            if count > best_count:
                best_count = count
                best_pa = pa
                best_center_ra = cand_ra[idx]
                best_center_dec = cand_dec[idx]
                best_mask = mask
                
        print(f"  Coarse best: center=({best_center_ra:.4f}, {best_center_dec:.4f}), PA={best_pa:.1f}°, covered={best_count} targets")
        
        # Fine search / Local refinement around the best coarse candidate
        if best_count > 0:
            fine_range = coarse_grid_step * 1.2
            fine_grid_ra = np.arange(best_center_ra - fine_range, best_center_ra + fine_range, fine_grid_step)
            fine_grid_dec = np.arange(best_center_dec - fine_range, best_center_dec + fine_range, fine_grid_step)
            
            fg_ra_mesh, fg_dec_mesh = np.meshgrid(fine_grid_ra, fine_grid_dec)
            fine_cand_ra = fg_ra_mesh.flatten()
            fine_cand_dec = fg_dec_mesh.flatten()
            
            # Also refine PA locally in steps of 0.5 degrees
            fine_pas = np.arange(best_pa - pa_step, best_pa + pa_step + 0.1, 0.5)
            
            for f_pa in fine_pas:
                f_pa_norm = f_pa % 60.0
                idx, count, mask = evaluate_candidates_chunk(
                    active_ra, active_dec, fine_cand_ra, fine_cand_dec, f_pa_norm
                )
                if count > best_count:
                    print(f"  Local refinement improved coverage: {best_count} -> {count} (PA: {best_pa:.1f}° -> {f_pa_norm:.1f}°)")
                    best_count = count
                    best_center_ra = fine_cand_ra[idx]
                    best_center_dec = fine_cand_dec[idx]
                    best_pa = f_pa_norm
                    best_mask = mask
            
            # Final active mask evaluation for chosen parameters
            _, _, best_mask = evaluate_candidates_chunk(
                active_ra, active_dec, np.array([best_center_ra]), np.array([best_center_dec]), best_pa
            )
                
        # Map the active mask (relative to uncovered targets) back to the full target list
        full_mask = np.zeros(N, dtype=bool)
        uncovered_indices = np.where(~covered)[0]
        full_mask[uncovered_indices[best_mask]] = True
        
        # Mark as covered
        covered[full_mask] = True
        
        pointings.append({
            'ppc_code': f"OPT_FOV_{fov_idx+1}",
            'ppc_ra': best_center_ra,
            'ppc_dec': best_center_dec,
            'ppc_pa': best_pa,
            'covered_count': best_count,
            'covered_target_ids': list(ids[full_mask])
        })
        
        print(f"  Final selected: center=({best_center_ra:.4f}, {best_center_dec:.4f}), PA={best_pa:.2f}°, newly covered={best_count} targets")
        
    return pointings, covered

def main():
    parser = argparse.ArgumentParser(description="Optimize PFS hexagon FoV pointings to maximize target coverage")
    parser.add_argument("--input", default="cosmos/targets_all_20260514.csv", help="Input targets CSV file")
    parser.add_argument("--max-priority", type=int, default=2, help="Filter targets with priority <= max_priority (smaller values mean higher priority)")
    parser.add_argument("--num-fovs", type=int, default=1, help="Number of FoVs to place")
    parser.add_argument("--pa-step", type=float, default=5.0, help="PA search step size in degrees (hexagon is 60-degree symmetric)")
    parser.add_argument("--output", default="optimized_pointings.ecsv", help="Output ECSV file path")
    parser.add_argument("--plot", default="optimized_coverage.png", help="Output PNG plot path")
    args = parser.parse_args()
    
    # 1. Read input CSV
    print(f"Reading target file: {args.input}...")
    if not os.path.exists(args.input):
        print(f"Error: file not found at {args.input}")
        return
        
    df = pd.read_csv(args.input)
    
    # Check necessary columns
    col_mapping = {}
    for col in ['ra', 'dec', 'priority', 'obj_id']:
        found = False
        for c in df.columns:
            if c.lower() == col or c.lower().replace('.', '').replace(' ', '_') == col:
                col_mapping[col] = c
                found = True
                break
        if not found:
            # Fallbacks
            if col == 'obj_id' and 'ob_code' in df.columns:
                col_mapping['obj_id'] = 'ob_code'
            elif col == 'obj_id' and 'ID' in df.columns:
                col_mapping['obj_id'] = 'ID'
            else:
                raise KeyError(f"Could not find equivalent column for '{col}' in input CSV.")
                
    df_clean = pd.DataFrame({
        'ra': df[col_mapping['ra']],
        'dec': df[col_mapping['dec']],
        'priority': df[col_mapping['priority']],
        'obj_id': df[col_mapping['obj_id']]
    })
    
    # 2. Filter by priority
    print(f"Filtering targets with priority <= {args.max_priority}...")
    df_filtered = df_clean[df_clean['priority'] <= args.max_priority].copy()
    print(f"Filtered target count: {len(df_filtered)}")
    
    if len(df_filtered) == 0:
        print("No targets match the priority filter. Exiting.")
        return
        
    # 3. Run Optimization
    pointings, covered = optimize_fovs(df_filtered, num_fovs=args.num_fovs, pa_step=args.pa_step)
    
    # 4. Save results to ECSV
    opt_table = Table(
        names=('ppc_code', 'ppc_ra', 'ppc_dec', 'ppc_pa', 'covered_count'),
        dtype=('S30', 'f8', 'f8', 'f8', 'i4')
    )
    for p in pointings:
        opt_table.add_row((p['ppc_code'], p['ppc_ra'], p['ppc_dec'], p['ppc_pa'], p['covered_count']))
        
    opt_table.write(args.output, format="ascii.ecsv", overwrite=True)
    print(f"\nSaved optimized pointings to {args.output}")
    
    # Print target allocation statistics
    total_covered = sum(p['covered_count'] for p in pointings)
    print(f"\nTotal targets of priority <= {args.max_priority} covered: {total_covered} / {len(df_filtered)} ({100*total_covered/len(df_filtered):.1f}%)")
    
    # 5. Create a beautiful plot
    print("Generating plot...")
    plt.style.use('default')
    fig, ax = plt.subplots(figsize=(10, 10), dpi=150)
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')
    
    # Plot uncovered targets as light semitransparent points
    uncovered_targets = df_filtered[~covered]
    ax.scatter(uncovered_targets['ra'], uncovered_targets['dec'], 
               color='#94a3b8', alpha=0.3, s=15, label='Uncovered Targets', marker='.')
               
    # Plot covered targets
    for p_idx, p in enumerate(pointings):
        covered_ids = p['covered_target_ids']
        df_cov = df_filtered[df_filtered['obj_id'].isin(covered_ids)]
        ax.scatter(df_cov['ra'], df_cov['dec'], 
                   alpha=0.8, s=25, label=f"Covered by FoV #{p_idx+1} ({len(df_cov)})")
                   
    # Plot Hexagons
    for idx, p in enumerate(pointings):
        ra_c, dec_c, pa = p['ppc_ra'], p['ppc_dec'], p['ppc_pa']
        v_ra, v_dec = get_vertices(ra_c, dec_c, R_HEX, pa)
        
        # Hexagon fill
        poly = patches.Polygon(
            np.column_stack((v_ra, v_dec)),
            closed=True,
            facecolor='#38bdf8',
            alpha=0.1,
            zorder=3
        )
        ax.add_patch(poly)
        
        # Hexagon outline
        ax.plot(
            v_ra, v_dec,
            color='#38bdf8',
            linewidth=2.0,
            linestyle='-',
            zorder=3
        )
        
        # Center marker and text
        ax.scatter(ra_c, dec_c, color='#38bdf8', s=40, marker='o', zorder=4)
        ax.text(
            ra_c, dec_c + 0.03, f"FoV #{idx+1}\nPA={pa:.1f}°",
            color='#0f172a', fontsize=9, fontweight='bold',
            ha='center', va='bottom', zorder=5
        )
        
    ax.set_xlabel('RA (deg)', fontsize=12)
    ax.set_ylabel('Dec (deg)', fontsize=12)
    ax.set_title(f'Optimized PFS Hexagonal FoV Placement\n({args.num_fovs} Fields, Max Priority <= {args.max_priority})', 
                 fontsize=14, fontweight='bold', pad=15)
                 
    ax.grid(True, color='#e2e8f0', linestyle='--', alpha=0.5)
    ax.set_aspect('equal')
    
    # Invert RA axis
    ax.invert_xaxis()
    
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(args.plot, bbox_inches='tight', dpi=150, facecolor=fig.get_facecolor(), edgecolor='none')
    plt.close()
    print(f"Saved plot to {args.plot}")

if __name__ == "__main__":
    main()
