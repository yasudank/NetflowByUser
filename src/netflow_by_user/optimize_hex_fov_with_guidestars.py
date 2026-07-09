#!/usr/bin/env python3
import os
import warnings
try:
    from erfa import ErfaWarning
    warnings.filterwarnings('ignore', category=ErfaWarning)
except ImportError:
    pass
import math
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.path as mppath
from astropy.table import Table
from pfs.utils.coordinates.CoordTransp import CoordinateTransform as ctrans
from ets_shuffle.convenience import flag_close_pairs, guidecam_geometry

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
    Evaluate target coverage for candidates in chunks to prevent memory blowup.
    Supports arbitrary rotation angle pa_deg.
    Returns array of covered counts for each candidate.
    """
    K = len(cand_ra)
    N = len(ra)
    covered_counts = np.zeros(K, dtype=int)
    
    cos_dec = np.cos(np.radians(cand_dec))
    
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
            
    return covered_counts

def evaluate_guidestars_single(ra_tel, dec_tel, pa_tel, df_gaia, obstime,
                               min_mag=12.0, max_mag=19.0, minsep_arcsec=1.0):
    """
    Evaluate guide star counts in each of the 6 guide cameras for a single pointing.
    Uses the guide star selection algorithm logic.
    """
    # 1. Filter Gaia stars by a circular search radius (approx. 1.1 degrees)
    # 260.0 * 10.2 / 3600 * 1.5 = 1.105 degrees
    search_radius = (260.0 * 10.2 / 3600.0) * 1.5
    
    dist2 = ((df_gaia['ra'] - ra_tel) * np.cos(np.radians(dec_tel)))**2 + (df_gaia['dec'] - dec_tel)**2
    df_sel = df_gaia[dist2 < search_radius**2].copy()
    
    if len(df_sel) == 0:
        return [0] * 6, []

    # 2. Run coordinate transform to get PFI coordinates (AG frame)
    epoch = df_sel["ref_epoch"].values
    gaiadb_epoch = epoch[0] if len(epoch) > 0 else 2016.0
    
    tmp = np.array([df_sel["ra"].values, df_sel["dec"].values])
    tmp = ctrans(
        xyin=tmp,
        mode="sky_pfi_ag",
        pa=pa_tel,
        cent=np.array([ra_tel, dec_tel]).reshape((2, 1)),
        pm=np.stack([df_sel["pmra"].values, df_sel["pmdec"].values], axis=0),
        par=df_sel["parallax"].values.copy(),
        time=obstime,
        epoch=gaiadb_epoch,
    )
    xypos = np.array([tmp[0, :], tmp[1, :]]).T
    
    # 3. Check guide camera polygons
    agcoord = guidecam_geometry()
    cam_star_counts = []
    selected_stars_all = []
    
    for i in range(agcoord.shape[0]):
        p = mppath.Path(agcoord[i])
        
        # Check which stars are in the slightly enlarged camera footprint (for neighbor checks)
        in_enlarged = p.contains_points(xypos, radius=1.0) # 1mm larger
        df_cam = df_sel[in_enlarged].copy()
        cam_xypos = xypos[in_enlarged]
        
        if len(df_cam) == 0:
            cam_star_counts.append(0)
            continue
            
        # Check for bright stars in the exact boundary causing saturation
        in_exact_all = p.contains_points(cam_xypos)
        if (df_cam[in_exact_all]["magnitude"] <= min_mag).any():
            cam_star_counts.append(-999)
            continue
            
        # Eliminate close neighbors
        flags_close = flag_close_pairs(df_cam["ra"].values, df_cam["dec"].values, minsep_arcsec / 3600.0)
        df_cam = df_cam[~flags_close]
        cam_xypos = cam_xypos[~flags_close]
        
        if len(df_cam) == 0:
            cam_star_counts.append(0)
            continue
            
        # Magnitude range filter
        in_mag_range = (df_cam["magnitude"] > min_mag) & (df_cam["magnitude"] < max_mag)
        df_cam = df_cam[in_mag_range]
        cam_xypos = cam_xypos[in_mag_range]
        
        if len(df_cam) == 0:
            cam_star_counts.append(0)
            continue
            
        # Exact boundary check
        in_exact = p.contains_points(cam_xypos)
        df_cam = df_cam[in_exact].copy()
        df_cam["agid"] = i
        
        cam_star_counts.append(len(df_cam))
        selected_stars_all.append(df_cam)
        
    if len(selected_stars_all) > 0:
        df_selected_stars = pd.concat(selected_stars_all, ignore_index=True)
    else:
        df_selected_stars = pd.DataFrame()
        
    return cam_star_counts, df_selected_stars

def score_pointing(target_count, cam_star_counts, min_stars_per_cam, min_cams_with_stars):
    """
    Compute candidate pointing score based on target coverage and guide star count.
    Score = target_count - 1,000,000 * max(0, min_cams_with_stars - cams_ok)
    """
    if any(count < 0 for count in cam_star_counts):
        return -10000000, 0
        
    cams_ok = sum(1 for count in cam_star_counts if count >= min_stars_per_cam)
    penalty = 1000000 * max(0, min_cams_with_stars - cams_ok)
    return target_count - penalty, cams_ok

def optimize_fovs_with_guidestars(df_targets, df_gaia, obstime, num_fovs=1,
                                 min_stars_per_cam=2, min_cams_with_stars=6,
                                 coarse_grid_step=0.03, fine_grid_step=0.002, pa_step=5.0,
                                 max_gs_checks=500):
    """
    Greedy maximum coverage solver with coarse-to-fine search and guide star constraint optimization
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
            
        # Coarse target coverage evaluation
        candidates_list = []
        for pa in coarse_pas:
            counts = evaluate_candidates_chunk(
                active_ra, active_dec, cand_ra, cand_dec, pa
            )
            for idx, count in enumerate(counts):
                if count > 0:
                    candidates_list.append({
                        'ra': cand_ra[idx],
                        'dec': cand_dec[idx],
                        'pa': pa,
                        'count': count
                    })
                    
        # Sort candidates by covered count descending
        candidates_list = sorted(candidates_list, key=lambda x: x['count'], reverse=True)
        print(f"Found {len(candidates_list)} coarse candidate pointings with coverage > 0.")
        
        # Evaluate guide star constraints for the top candidates
        best_score = -999999999
        best_pa = -1.0
        best_center_ra = -1.0
        best_center_dec = -1.0
        best_target_count = -1
        best_cams_ok = 0
        best_star_counts = [0]*6
        best_stars_df = pd.DataFrame()
        
        checks = 0
        for cand in candidates_list[:max_gs_checks]:
            c_ra, c_dec, c_pa, count = cand['ra'], cand['dec'], cand['pa'], cand['count']
            
            # Get guide star counts
            star_counts, stars_df = evaluate_guidestars_single(
                c_ra, c_dec, c_pa, df_gaia, obstime,
                min_mag=12.0, max_mag=21.5, minsep_arcsec=1.0
            )
            score, cams_ok = score_pointing(count, star_counts, min_stars_per_cam, min_cams_with_stars)
            
            if score > best_score:
                best_score = score
                best_pa = c_pa
                best_center_ra = c_ra
                best_center_dec = c_dec
                best_target_count = count
                best_cams_ok = cams_ok
                best_star_counts = star_counts
                best_stars_df = stars_df
                
            # If we fully satisfied the guide star constraint (penalty = 0), we can stop checking further down,
            # because candidates are sorted by target coverage, and any subsequent candidate will have <= target coverage
            if cams_ok >= min_cams_with_stars:
                print(f"  Coarse search short-circuited after {checks+1} checks.")
                break
                
            checks += 1
            
        print(f"  Coarse best: center=({best_center_ra:.4f}, {best_center_dec:.4f}), PA={best_pa:.1f}°, covered={best_target_count} targets")
        print(f"  Guide stars check: cams ok={best_cams_ok}/{min_cams_with_stars}, counts={best_star_counts}")
        
        # Fine search / Local refinement around the best coarse candidate
        if best_target_count > 0:
            fine_range = coarse_grid_step * 1.2
            fine_grid_ra = np.arange(best_center_ra - fine_range, best_center_ra + fine_range, fine_grid_step)
            fine_grid_dec = np.arange(best_center_dec - fine_range, best_center_dec + fine_range, fine_grid_step)
            
            fg_ra_mesh, fg_dec_mesh = np.meshgrid(fine_grid_ra, fine_grid_dec)
            fine_cand_ra = fg_ra_mesh.flatten()
            fine_cand_dec = fg_dec_mesh.flatten()
            
            # Refine PA locally
            fine_pas = np.arange(best_pa - pa_step, best_pa + pa_step + 0.1, 0.5)
            
            fine_candidates = []
            for f_pa in fine_pas:
                f_pa_norm = f_pa % 60.0
                counts = evaluate_candidates_chunk(
                    active_ra, active_dec, fine_cand_ra, fine_cand_dec, f_pa_norm
                )
                for idx, count in enumerate(counts):
                    if count > 0:
                        fine_candidates.append({
                            'ra': fine_cand_ra[idx],
                            'dec': fine_cand_dec[idx],
                            'pa': f_pa_norm,
                            'count': count
                        })
                        
            # Sort fine candidates by target count descending
            fine_candidates = sorted(fine_candidates, key=lambda x: x['count'], reverse=True)
            
            # Evaluate guide stars for the top fine candidates
            fine_checks = 0
            for cand in fine_candidates[:max_gs_checks]:
                f_ra, f_dec, f_pa, count = cand['ra'], cand['dec'], cand['pa'], cand['count']
                
                # We can skip evaluating if the count is smaller than a valid solution we already have
                # and our best score so far is already a fully valid solution (no penalty)
                if best_cams_ok >= min_cams_with_stars and count < best_target_count:
                    break
                    
                star_counts, stars_df = evaluate_guidestars_single(
                    f_ra, f_dec, f_pa, df_gaia, obstime,
                    min_mag=12.0, max_mag=19.0, minsep_arcsec=1.0
                )
                score, cams_ok = score_pointing(count, star_counts, min_stars_per_cam, min_cams_with_stars)
                
                if score > best_score:
                    print(f"  Local refinement improved score: target count {best_target_count} -> {count}, cams ok {best_cams_ok} -> {cams_ok}")
                    best_score = score
                    best_pa = f_pa
                    best_center_ra = f_ra
                    best_center_dec = f_dec
                    best_target_count = count
                    best_cams_ok = cams_ok
                    best_star_counts = star_counts
                    best_stars_df = stars_df
                    
                if cams_ok >= min_cams_with_stars and count >= best_target_count:
                    print(f"  Fine search short-circuited after {fine_checks+1} checks.")
                    break
                    
                fine_checks += 1
                
            # Re-evaluate candidate target coverage mask for chosen best parameters
            covered_counts = evaluate_candidates_chunk(
                active_ra, active_dec, np.array([best_center_ra]), np.array([best_center_dec]), best_pa
            )
            
            # Find which active targets are covered
            cos_dec_tel = np.cos(np.radians(best_center_dec))
            dx = (active_ra - best_center_ra) * cos_dec_tel
            dy = active_dec - best_center_dec
            
            rad = np.radians(best_pa)
            cos_pa = np.cos(rad)
            sin_pa = np.sin(rad)
            dx_rot = dx * cos_pa + dy * sin_pa
            dy_rot = -dx * sin_pa + dy * cos_pa
            
            cond1 = np.abs(dy_rot) <= R_HEX * np.sqrt(3)/2.0
            cond2 = np.abs(dy_rot) + np.sqrt(3) * np.abs(dx_rot) <= np.sqrt(3) * R_HEX
            best_mask = cond1 & cond2
            
        # Map active mask back to the full target list
        full_mask = np.zeros(N, dtype=bool)
        uncovered_indices = np.where(~covered)[0]
        full_mask[uncovered_indices[best_mask]] = True
        
        # Mark targets as covered
        covered[full_mask] = True
        
        pointings.append({
            'ppc_code': f"OPT_FOV_{fov_idx+1}",
            'ppc_ra': best_center_ra,
            'ppc_dec': best_center_dec,
            'ppc_pa': best_pa,
            'covered_count': best_target_count,
            'covered_target_ids': list(ids[full_mask]),
            'star_counts': best_star_counts,
            'stars_df': best_stars_df
        })
        
        print(f"  Final selected: center=({best_center_ra:.4f}, {best_center_dec:.4f}), PA={best_pa:.2f}°, newly covered={best_target_count} targets")
        print(f"  Guide stars in final pointing: cams ok={best_cams_ok}/{min_cams_with_stars}, counts={best_star_counts}")
        
    return pointings, covered

def plot_optimized_fovs(pointings, covered, df_filtered, max_priority, plot_path, obstime):
    print("Generating plot...")
    plt.style.use('default')
    fig, ax = plt.subplots(figsize=(10, 10), dpi=150)
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')
    
    # Plot uncovered targets
    uncovered_targets = df_filtered[~covered]
    ax.scatter(uncovered_targets['ra'], uncovered_targets['dec'], 
               color='#94a3b8', alpha=0.3, s=15, label='Uncovered Targets', marker='.')
               
    # Plot covered targets
    for p_idx, p in enumerate(pointings):
        covered_ids = p['covered_target_ids']
        df_cov = df_filtered[df_filtered['obj_id'].isin(covered_ids)]
        ax.scatter(df_cov['ra'], df_cov['dec'], 
                   alpha=0.8, s=25, label=f"Covered by FoV #{p_idx+1} ({len(df_cov)})")
                   
    # Plot Hexagons, Guide Camera footprints, and Guide Stars
    agcoord = guidecam_geometry()
    
    for idx, p in enumerate(pointings):
        ra_c, dec_c, pa = p['ppc_ra'], p['ppc_dec'], p['ppc_pa']
        
        # 1. Hexagon footprint
        v_ra, v_dec = get_vertices(ra_c, dec_c, R_HEX, pa)
        poly = patches.Polygon(
            np.column_stack((v_ra, v_dec)),
            closed=True, facecolor='#38bdf8', alpha=0.08, zorder=2
        )
        ax.add_patch(poly)
        ax.plot(v_ra, v_dec, color='#38bdf8', linewidth=1.5, linestyle='-', zorder=2)
        
        # 2. Guide camera footprints
        for cam_idx in range(agcoord.shape[0]):
            cam_pfi = agcoord[cam_idx].T # Shape (2, 4)
            tmp_sky = ctrans(
                xyin=cam_pfi,
                mode="pfi_sky",
                pa=pa,
                cent=np.array([ra_c, dec_c]).reshape((2, 1)),
                time=obstime,
                epoch=2016.0
            )
            v_cam_ra, v_cam_dec = tmp_sky[0, :], tmp_sky[1, :]
            
            cam_poly = patches.Polygon(
                np.column_stack((v_cam_ra, v_cam_dec)),
                closed=True, facecolor='#22c55e', alpha=0.15, zorder=3
            )
            ax.add_patch(cam_poly)
            ax.plot(v_cam_ra, v_cam_dec, color='#22c55e', linewidth=1.0, linestyle='--', zorder=3)
            
            # Label camera id
            label_ra, label_dec = np.mean(v_cam_ra), np.mean(v_cam_dec)
            dra = label_ra - ra_c
            ddec = label_dec - dec_c
            dist = np.hypot(dra, ddec)
            if dist > 0:
                label_ra += (dra / dist) * 0.06
                label_dec += (ddec / dist) * 0.06
            ax.text(
                label_ra, label_dec, f"AG{cam_idx}",
                color='#22c55e', fontsize=6, fontweight='bold',
                ha='center', va='center', zorder=4
            )
            
        # 3. Plot selected guide stars
        stars_df = p['stars_df']
        if len(stars_df) > 0:
            ax.scatter(stars_df['ra'], stars_df['dec'], color='#eab308', s=45, marker='*', edgecolors='black', linewidths=0.5, zorder=5, label='Selected Guide Stars' if idx==0 else "")
            
        # Center marker and text
        ax.scatter(ra_c, dec_c, color='#38bdf8', s=40, marker='o', zorder=4)
        ax.text(
            ra_c, dec_c + 0.03, f"FoV #{idx+1}\nPA={pa:.1f}°\nGS counts={p['star_counts']}",
            color='#0f172a', fontsize=8, fontweight='bold',
            ha='center', va='bottom', zorder=5
        )
        
    ax.set_xlabel('RA (deg)', fontsize=12)
    ax.set_ylabel('Dec (deg)', fontsize=12)

    # Calculate coverage statistics per priority level (priority <= max_priority)
    stats_str = ""
    priorities = sorted(df_filtered['priority'].unique())
    for prio in priorities:
        sub = df_filtered[df_filtered['priority'] == prio]
        sub_covered = covered[df_filtered['priority'] == prio]
        cov_count = np.sum(sub_covered)
        tot_count = len(sub)
        pct = (100.0 * cov_count / tot_count) if tot_count > 0 else 0.0
        stats_str += f"P{prio}: {cov_count}/{tot_count} ({pct:.1f}%)   "

    title_text = f'Optimized PFS FoV with Guide Star Constraints\n({len(pointings)} Fields, Priority <= {max_priority})\n{stats_str.strip()}'
    ax.set_title(title_text, fontsize=12, fontweight='bold', pad=15)
                 
    ax.grid(True, color='#e2e8f0', linestyle='--', alpha=0.5)
    ax.set_aspect('equal')
    ax.invert_xaxis()
    
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(plot_path, bbox_inches='tight', dpi=150, facecolor=fig.get_facecolor(), edgecolor='none')
    plt.close()
    print(f"Saved plot to {plot_path}")

def main():
    parser = argparse.ArgumentParser(description="Optimize PFS hexagon FoV pointings to maximize target coverage and satisfy guide star constraints")
    parser.add_argument("--input", default="cosmos/targets_all_20260514.csv", help="Input targets CSV file")
    parser.add_argument("--gaia-catalog", default="cosmos/gaia.ecsv", help="Path to local Gaia ECSV catalog for guide stars")
    parser.add_argument("--obstime", default="2026-05-09T06:00:00Z", help="Observing time in UTC")
    parser.add_argument("--min-stars-per-cam", type=int, default=2, help="Minimum guide stars per camera")
    parser.add_argument("--min-cams-with-stars", type=int, default=6, help="Minimum cameras satisfying the requirement")
    parser.add_argument("--max-priority", type=int, default=2, help="Filter targets with priority <= max_priority")
    parser.add_argument("--num-fovs", type=int, default=1, help="Number of FoVs to place")
    parser.add_argument("--pa-step", type=float, default=5.0, help="PA search step size in degrees")
    parser.add_argument("--max-gs-checks", type=int, default=500, help="Max candidates checked for guide stars")
    parser.add_argument("--output", default="optimized_pointings_with_gs.ecsv", help="Output ECSV file path")
    parser.add_argument("--plot", default="optimized_coverage_with_gs.png", help="Output PNG plot path")
    args = parser.parse_args()
    
    # 1. Read input CSV
    print(f"Reading target file: {args.input}...")
    if not os.path.exists(args.input):
        print(f"Error: target file not found at {args.input}")
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
        
    # 3. Read Gaia Catalog
    print(f"Reading Gaia catalog from: {args.gaia_catalog}...")
    if not os.path.exists(args.gaia_catalog):
        print(f"Error: Gaia catalog not found at {args.gaia_catalog}")
        return
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
    
    # 4. Run Optimization
    pointings, covered = optimize_fovs_with_guidestars(
        df_filtered, df_gaia, args.obstime, num_fovs=args.num_fovs,
        min_stars_per_cam=args.min_stars_per_cam, min_cams_with_stars=args.min_cams_with_stars,
        pa_step=args.pa_step, max_gs_checks=args.max_gs_checks
    )
    
    # 5. Save results to ECSV
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
    
    # 6. Generate Plot
    plot_optimized_fovs(pointings, covered, df_filtered, args.max_priority, args.plot, args.obstime)

if __name__ == "__main__":
    main()
