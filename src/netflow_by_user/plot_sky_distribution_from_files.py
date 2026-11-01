#!/usr/bin/env python3
"""
Plot the sky distribution of assigned targets, calibration stars, sky positions,
and unassigned high-priority targets from saved files, including the Rubin/LSST FoV.
"""
import os
import sys
import glob
import argparse
from collections import defaultdict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from astropy.table import Table

# Regular hexagon parameters for PFS
R_HEX = 0.69  # Outer radius of PFS FoV in degrees

def get_vertices(ra, dec, r, pa_deg):
    """Get the 6 vertices of a rotated PFS hexagon centered at (ra, dec)."""
    cos_dec = np.cos(np.radians(dec))
    angles = np.radians([60 * i for i in range(7)])
    v_dx = r * np.cos(angles)
    v_dy = r * np.sin(angles)
    
    rad = np.radians(pa_deg)
    cos_pa = np.cos(rad)
    sin_pa = np.sin(rad)
    
    v_dx_rot = v_dx * cos_pa - v_dy * sin_pa
    v_dy_rot = v_dx * sin_pa + v_dy * cos_pa
    
    v_ra = ra + v_dx_rot / cos_dec
    v_dec = dec + v_dy_rot
    return v_ra, v_dec

def get_circular_fov_boundary(center_ra, center_dec, radius_deg=1.75, n_points=360):
    """
    Generate (ra, dec) coordinates along a circular FoV of radius radius_deg on the celestial sphere.
    Default radius is 1.75 degrees (3.5 degree diameter for Rubin/LSST).
    """
    theta = np.linspace(0, 2 * np.pi, n_points)
    r_rad = np.radians(radius_deg)
    dec0_rad = np.radians(center_dec)
    ra0_rad = np.radians(center_ra)
    
    # Spherical trigonometry to compute circle coordinates
    sin_dec = np.sin(dec0_rad) * np.cos(r_rad) + np.cos(dec0_rad) * np.sin(r_rad) * np.cos(theta)
    dec_rad = np.arcsin(np.clip(sin_dec, -1.0, 1.0))
    
    cos_dec = np.cos(dec_rad)
    sin_dra = (np.sin(r_rad) * np.sin(theta)) / np.where(cos_dec == 0, 1e-10, cos_dec)
    cos_dra = (np.cos(r_rad) - np.sin(dec0_rad) * sin_dec) / (np.cos(dec0_rad) * np.where(cos_dec == 0, 1e-10, cos_dec))
    dra_rad = np.arctan2(sin_dra, cos_dra)
    
    ra_rad = ra0_rad + dra_rad
    return np.degrees(ra_rad), np.degrees(dec_rad)

def compute_pointing_mask(ra, dec, center_ra, center_dec, pa_deg, r=R_HEX):
    """Compute boolean mask indicating whether (ra, dec) points fall inside the rotated PFS hexagon."""
    cos_dec = np.cos(np.radians(center_dec))
    dx = (ra - center_ra) * cos_dec
    dy = dec - center_dec
    rad = np.radians(pa_deg)
    cos_pa = np.cos(rad)
    sin_pa = np.sin(rad)
    dx_rot = dx * cos_pa + dy * sin_pa
    dy_rot = -dx * sin_pa + dy * cos_pa
    cond1 = np.abs(dy_rot) <= r * np.sqrt(3) / 2.0
    cond2 = np.abs(dy_rot) + np.sqrt(3) * np.abs(dx_rot) <= np.sqrt(3) * r
    return cond1 & cond2

def sort_key(c):
    if c.startswith('sci_P'):
        return (0, int(c.split('_P')[1]))
    elif c == 'cal':
        return (1, 0)
    else:
        return (2, 0)

def plot_sky_distribution_from_files(
    targets_dir="targets",
    science_targets_file=None,
    pointing_file=None,
    output_plot="sky_distribution_with_rubin.png",
    show_rubin=True,
    rubin_center_ra=None,
    rubin_center_dec=None,
    rubin_radius=1.75,
    show_pfs_hex=True,
    title=None
):
    print("\n================ Loading Assigned Targets ================")
    plot_data = defaultdict(lambda: {'ra': [], 'dec': []})
    assigned_obj_ids = set()

    # If pointing_file is provided, determine target ppc_codes to avoid reading leftover files
    target_ppc_codes = None
    if pointing_file and os.path.exists(pointing_file):
        try:
            t_ptg = Table.read(pointing_file, format="ascii.ecsv").to_pandas()
            for col in ['ppc_code', 'pointing_id', 'code']:
                if col in t_ptg.columns:
                    target_ppc_codes = set(t_ptg[col].astype(str).tolist())
                    print(f"Filtering targets by {len(target_ppc_codes)} pointings specified in {pointing_file}.")
                    break
        except Exception as e:
            print(f"Warning reading pointing file {pointing_file}: {e}")

    # 1. Load assigned science targets
    sci_dir = os.path.join(targets_dir, "science")
    sci_files = sorted(glob.glob(os.path.join(sci_dir, "*.ecsv")))
    # Filter out .orig files
    sci_files = [f for f in sci_files if not f.endswith('.orig')]
    if target_ppc_codes:
        sci_files = [f for f in sci_files if os.path.splitext(os.path.basename(f))[0] in target_ppc_codes]
    print(f"Found {len(sci_files)} science target files in {sci_dir}")

    per_pointing_assigned_p0 = defaultdict(int)
    per_pointing_assigned_total = defaultdict(int)

    for fpath in sci_files:
        try:
            p_code = os.path.splitext(os.path.basename(fpath))[0]
            t = Table.read(fpath, format="ascii.ecsv")
            df = t.to_pandas()
            for _, row in df.iterrows():
                obj_id = str(row.get('obj_id', ''))
                assigned_obj_ids.add(obj_id)
                prio = int(row.get('priority', 9))
                cls = f"sci_P{prio}"
                plot_data[cls]['ra'].append(float(row['ra']))
                plot_data[cls]['dec'].append(float(row['dec']))
                
                per_pointing_assigned_total[p_code] += 1
                if prio == 0:
                    per_pointing_assigned_p0[p_code] += 1
        except Exception as e:
            print(f"  Warning reading {fpath}: {e}")

    # 2. Load assigned fluxstd targets
    cal_dir = os.path.join(targets_dir, "fluxstd")
    cal_files = sorted(glob.glob(os.path.join(cal_dir, "*.ecsv")))
    if target_ppc_codes:
        cal_files = [f for f in cal_files if os.path.splitext(os.path.basename(f))[0] in target_ppc_codes]
    print(f"Found {len(cal_files)} fluxstd files in {cal_dir}")
    for fpath in cal_files:
        try:
            t = Table.read(fpath, format="ascii.ecsv")
            df = t.to_pandas()
            for _, row in df.iterrows():
                plot_data['cal']['ra'].append(float(row['ra']))
                plot_data['cal']['dec'].append(float(row['dec']))
        except Exception as e:
            print(f"  Warning reading {fpath}: {e}")

    # 3. Load assigned sky targets
    sky_dir = os.path.join(targets_dir, "sky")
    sky_files = sorted(glob.glob(os.path.join(sky_dir, "*.ecsv")))
    if target_ppc_codes:
        sky_files = [f for f in sky_files if os.path.splitext(os.path.basename(f))[0] in target_ppc_codes]
    print(f"Found {len(sky_files)} sky files in {sky_dir}")
    for fpath in sky_files:
        try:
            t = Table.read(fpath, format="ascii.ecsv")
            df = t.to_pandas()
            for _, row in df.iterrows():
                plot_data['sky']['ra'].append(float(row['ra']))
                plot_data['sky']['dec'].append(float(row['dec']))
        except Exception as e:
            print(f"  Warning reading {fpath}: {e}")

    # 4. Find unassigned high priority (P0, P1, P2) science targets and compute P0 statistics
    unassigned_high_prio = {'ra': [], 'dec': []}
    df_p0 = None
    if science_targets_file and os.path.exists(science_targets_file):
        print(f"Loading input science targets from {science_targets_file}...")
        df_all_sci = pd.read_csv(science_targets_file)
        
        # Determine obj_id / id column
        id_col = None
        for c in ['obj_id', 'ob_code', 'target_id', 'ID', 'id']:
            if c in df_all_sci.columns:
                id_col = c
                break
                
        prio_col = None
        for c in ['priority', 'prio', 'Priority']:
            if c in df_all_sci.columns:
                prio_col = c
                break
                
        ra_col = 'ra' if 'ra' in df_all_sci.columns else 'RA'
        dec_col = 'dec' if 'dec' in df_all_sci.columns else 'DEC'
        
        if id_col and prio_col:
            df_high_prio = df_all_sci[df_all_sci[prio_col].isin([0, 1, 2])]
            for _, row in df_high_prio.iterrows():
                t_id = str(row[id_col])
                if t_id not in assigned_obj_ids:
                    unassigned_high_prio['ra'].append(float(row[ra_col]))
                    unassigned_high_prio['dec'].append(float(row[dec_col]))
            print(f"Found {len(unassigned_high_prio['ra'])} unassigned P0, P1, P2 targets (out of {len(df_high_prio)}).")
            
            df_p0 = df_all_sci[df_all_sci[prio_col] == 0].copy()
        else:
            print("  Warning: Required columns for ID or Priority not found in science targets catalog.")
    else:
        print("Input science targets catalog not specified or found; skipping unassigned target overlay.")

    # Sort classes by priority
    classes = sorted(list(plot_data.keys()), key=sort_key)
    print("\n--- Summary of Assigned Objects ---")
    total_assigned = 0
    for cls in classes:
        n = len(plot_data[cls]['ra'])
        total_assigned += n
        print(f"  {cls:<10}: {n:5d}")
    print(f"  Total Assigned: {total_assigned}")
    print("-----------------------------------")

    # Determine Rubin / LSST Center coordinates if not specified
    if show_rubin:
        if rubin_center_ra is None or rubin_center_dec is None:
            # Check if pointing file exists and inspect centers
            if pointing_file and os.path.exists(pointing_file):
                t_ptg = Table.read(pointing_file, format="ascii.ecsv").to_pandas()
                ra_c = float(t_ptg['ppc_ra'].mean())
                dec_c = float(t_ptg['ppc_dec'].mean())
            elif plot_data:
                # Use mean of assigned science targets
                all_ra = []
                all_dec = []
                for cls in classes:
                    all_ra.extend(plot_data[cls]['ra'])
                    all_dec.extend(plot_data[cls]['dec'])
                ra_c = float(np.mean(all_ra))
                dec_c = float(np.mean(all_dec))
            else:
                ra_c, dec_c = 35.70, -4.80  # Default XMM-LSS
                
            # If coordinates are close to XMM-LSS (RA ~ 35.7, Dec ~ -4.8)
            if abs(ra_c - 35.7) < 5.0 and abs(dec_c - (-4.8)) < 5.0:
                rubin_center_ra = 35.70 if rubin_center_ra is None else rubin_center_ra
                rubin_center_dec = -4.80 if rubin_center_dec is None else rubin_center_dec
            elif abs(ra_c - 150.1) < 5.0 and abs(dec_c - 2.2) < 5.0:
                # COSMOS DDF
                rubin_center_ra = 150.10 if rubin_center_ra is None else rubin_center_ra
                rubin_center_dec = 2.18 if rubin_center_dec is None else rubin_center_dec
            else:
                rubin_center_ra = ra_c if rubin_center_ra is None else rubin_center_ra
                rubin_center_dec = dec_c if rubin_center_dec is None else rubin_center_dec

    # --- Compute per-field Priority 0 (P0) Statistics ---
    per_fov_p0_stats = {}
    if df_p0 is not None and len(df_p0) > 0:
        p0_ra = df_p0[ra_col].values
        p0_dec = df_p0[dec_col].values
        p0_ids = set(df_p0[id_col].astype(str).values)
        total_p0_catalog = len(df_p0)
        
        pfs_covered_p0_mask = np.zeros(total_p0_catalog, dtype=bool)
        
        if pointing_file and os.path.exists(pointing_file):
            t_ptg_df = Table.read(pointing_file, format="ascii.ecsv").to_pandas()
            for _, r_ptg in t_ptg_df.iterrows():
                p_code = str(r_ptg.get('ppc_code', ''))
                ra_ptg = float(r_ptg['ppc_ra'])
                dec_ptg = float(r_ptg['ppc_dec'])
                pa_ptg = float(r_ptg.get('ppc_pa', 0.0))
                
                mask_in_fov = compute_pointing_mask(p0_ra, p0_dec, ra_ptg, dec_ptg, pa_ptg)
                pfs_covered_p0_mask |= mask_in_fov
                
                n_avail = int(mask_in_fov.sum())
                n_assigned = per_pointing_assigned_p0.get(p_code, 0)
                per_fov_p0_stats[p_code] = {
                    'available': n_avail,
                    'assigned': n_assigned,
                    'pct': (100.0 * n_assigned / n_avail) if n_avail > 0 else 0.0
                }
                
        # Total uniquely assigned P0
        total_unique_p0_assigned = len(assigned_obj_ids.intersection(p0_ids))
        total_pfs_covered_p0 = int(pfs_covered_p0_mask.sum())
        
        print("\n================ Priority 0 (P0) Target Statistics ================")
        print(f"Catalog Total P0 Targets: {total_p0_catalog}")
        print(f"PFS FoV Covered Unique P0: {total_pfs_covered_p0} ({100.0 * total_pfs_covered_p0 / total_p0_catalog:.1f}% of catalog)")
        print(f"PFS Total Unique P0 Assigned: {total_unique_p0_assigned} ({100.0 * total_unique_p0_assigned / total_p0_catalog:.1f}% of catalog, {100.0 * total_unique_p0_assigned / max(1, total_pfs_covered_p0):.1f}% of covered)")
        print("-------------------------------------------------------------------")
        print("PFS Pointing FoVs:")
        for p_code, s in per_fov_p0_stats.items():
            print(f"  {p_code:<22}: Available in FoV = {s['available']:3d} | Assigned = {s['assigned']:3d} ({s['pct']:5.1f}%)")
            
        if show_rubin and rubin_center_ra is not None:
            # Distance on sphere to Rubin center
            cos_rubin_dec = np.cos(np.radians(rubin_center_dec))
            dist_to_rubin = np.hypot((p0_ra - rubin_center_ra) * cos_rubin_dec, p0_dec - rubin_center_dec)
            in_rubin_mask = (dist_to_rubin <= rubin_radius)
            n_rubin_p0 = int(in_rubin_mask.sum())
            
            p0_ids_in_rubin = set(df_p0[id_col][in_rubin_mask].astype(str).values)
            n_rubin_p0_assigned = len(assigned_obj_ids.intersection(p0_ids_in_rubin))
            print("-------------------------------------------------------------------")
            print(f"Rubin/LSST Field of View (Center: RA={rubin_center_ra:.2f}°, Dec={rubin_center_dec:.2f}°, Radius={rubin_radius:.2f}°):")
            print(f"  Available P0 in Rubin FoV: {n_rubin_p0} ({100.0 * n_rubin_p0 / total_p0_catalog:.1f}% of catalog)")
            print(f"  Assigned P0 in Rubin FoV : {n_rubin_p0_assigned} ({100.0 * n_rubin_p0_assigned / max(1, n_rubin_p0):.1f}% of Rubin available)")
        print("===================================================================\n")

    # Color map for priorities (matching netflow_plot.py)
    colors = {
        'sci_P0': 'red',
        'sci_P1': 'orange',
        'sci_P2': 'gold',
        'sci_P3': 'yellowgreen',
        'sci_P4': 'green',
        'sci_P6': 'cyan',
        'sci_P7': 'dodgerblue',
        'sci_P9': 'navy',
        'cal': 'magenta',
        'sky': 'gray'
    }

    plt.style.use('default')
    fig, ax = plt.subplots(figsize=(12, 10), dpi=200)
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')

    # 1. Plot Rubin / LSST Field of View
    if show_rubin:
        r_ra, r_dec = get_circular_fov_boundary(rubin_center_ra, rubin_center_dec, radius_deg=rubin_radius)
        rubin_poly = patches.Polygon(
            np.column_stack((r_ra, r_dec)),
            closed=True,
            facecolor='#f43f5e',
            alpha=0.06,
            zorder=1
        )
        ax.add_patch(rubin_poly)
        ax.plot(
            r_ra, r_dec,
            color='#e11d48',
            linewidth=2.2,
            linestyle='--',
            label=f"Rubin/LSST FoV (r={rubin_radius}°)",
            zorder=4
        )
        ax.scatter(rubin_center_ra, rubin_center_dec, color='#e11d48', marker='+', s=60, linewidth=1.5, zorder=5)

    # 2. Plot PFS Hexagonal FoVs if pointing file is available
    if show_pfs_hex and pointing_file and os.path.exists(pointing_file):
        try:
            t_ptg = Table.read(pointing_file, format="ascii.ecsv").to_pandas()
            for idx, row in t_ptg.iterrows():
                p_code = str(row.get('ppc_code', f"FOV_{idx+1}"))
                ra_c = float(row['ppc_ra'])
                dec_c = float(row['ppc_dec'])
                pa_c = float(row.get('ppc_pa', 0.0))
                
                v_ra, v_dec = get_vertices(ra_c, dec_c, R_HEX, pa_c)
                hex_poly = patches.Polygon(
                    np.column_stack((v_ra, v_dec)),
                    closed=True,
                    facecolor='#0284c7',
                    alpha=0.04,
                    zorder=2
                )
                ax.add_patch(hex_poly)
                lbl = "PFS Hexagonal FoV (7 Fields)" if idx == 0 else None
                ax.plot(v_ra, v_dec, color='#0284c7', linewidth=1.2, linestyle='-', label=lbl, zorder=3)
                if p_code in per_fov_p0_stats:
                    s_p0 = per_fov_p0_stats[p_code]
                    label_text = f"{p_code}\nP0: {s_p0['assigned']}/{s_p0['available']}"
                else:
                    label_text = f"{p_code}"
                ax.text(
                    ra_c, dec_c, label_text,
                    color='#0369a1', fontsize=7, fontweight='bold',
                    ha='center', va='center', alpha=0.9, zorder=6,
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.6, edgecolor='none')
                )
        except Exception as e:
            print(f"  Warning plotting PFS hexagons: {e}")

    # 3. Plot Unassigned High Priority targets (background)
    if unassigned_high_prio['ra']:
        ax.scatter(
            unassigned_high_prio['ra'], unassigned_high_prio['dec'],
            label=f"Unassigned P0, P1, P2 ({len(unassigned_high_prio['ra'])})",
            color='black', marker='x', alpha=0.5, s=20, zorder=5
        )

    # 4. Plot lower priorities first (sky, cal, P9) so high priorities are on top
    for cls in reversed(classes):
        color = colors.get(cls, 'black')
        alpha = 0.25 if cls in ['sky', 'sci_P9'] else 0.85
        
        if cls == 'sky':
            marker_size = 5
        elif cls == 'cal':
            marker_size = 10
        elif cls.startswith('sci_P'):
            pri = int(cls.split('_P')[1])
            marker_size = max(6, 50 - pri * 5)
        else:
            marker_size = 15
        
        ax.scatter(
            plot_data[cls]['ra'], plot_data[cls]['dec'],
            label=f"Assigned {cls} ({len(plot_data[cls]['ra'])})",
            color=color, alpha=alpha, s=marker_size, edgecolors='none', zorder=7
        )

    ax.set_xlabel('RA (deg)', fontsize=12)
    ax.set_ylabel('Dec (deg)', fontsize=12)
    
    if title is None:
        title = "Sky Distribution of Targeted Objects & Rubin/LSST Field of View"
    ax.set_title(title, fontsize=13, fontweight='bold', pad=12)

    # Legend outside plot
    ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), title="Legend", frameon=True, fontsize=9)
    ax.grid(True, linestyle='--', alpha=0.5)

    # Invert x-axis to represent RA properly (east to the left)
    ax.invert_xaxis()
    ax.set_aspect('equal')

    plt.tight_layout()
    plt.savefig(output_plot, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"\nSuccessfully generated plot: {output_plot}")

def main():
    # 0. Pre-parse config if provided
    conf_parser = argparse.ArgumentParser(add_help=False)
    conf_parser.add_argument("--config", "-c", help="Path to YAML configuration file (e.g. netflow_config_2026Nov.yaml)")
    args_conf, remaining_argv = conf_parser.parse_known_args()

    defaults = {
        "targets_dir": "targets",
        "science_targets": None,
        "pointing_file": None,
        "output": "sky_distribution_with_rubin.png",
        "show_rubin": True,
        "rubin_ra": None,
        "rubin_dec": None,
        "rubin_radius": 1.75,
        "show_pfs_hex": True,
        "title": None,
    }

    if args_conf.config and os.path.exists(args_conf.config):
        import yaml
        print(f"Loading configurations from {args_conf.config}...")
        with open(args_conf.config, "r") as f:
            cfg = yaml.safe_load(f)
        if cfg:
            inputs = cfg.get("inputs", {})
            cat_dir = inputs.get("catalog_dir", ".")
            sci_tgt = inputs.get("science_targets")
            if sci_tgt:
                defaults["science_targets"] = os.path.join(cat_dir, sci_tgt) if not os.path.isabs(sci_tgt) else sci_tgt
                
            ptg = inputs.get("pointing_file")
            if ptg:
                defaults["pointing_file"] = ptg
                
            outputs = cfg.get("outputs", {})
            targets_dir = outputs.get("targets_dir")
            if targets_dir:
                defaults["targets_dir"] = targets_dir
                
            plot_file = outputs.get("plot_file")
            if plot_file:
                # Add suffix or use plot_file name
                base, ext = os.path.splitext(plot_file)
                defaults["output"] = f"{base}_with_rubin{ext}"

    parser = argparse.ArgumentParser(
        description="Plot sky distribution of assigned targets and unassigned high priority targets with Rubin/LSST FoV."
    )
    parser.add_argument("--config", "-c", help="Path to YAML configuration file")
    parser.add_argument("--targets-dir", default=defaults["targets_dir"], help="Path to directory containing partitioned targets (default: targets)")
    parser.add_argument("--science-targets", default=defaults["science_targets"], help="Path to input full science targets CSV (for identifying unassigned P0..P2)")
    parser.add_argument("--pointing-file", "-p", default=defaults["pointing_file"], help="Path to PFS pointing list ECSV (to draw PFS hexagonal FoVs)")
    parser.add_argument("--output", "-o", default=defaults["output"], help="Output PNG path (default: sky_distribution_with_rubin.png)")
    parser.add_argument("--rubin-ra", type=float, default=defaults["rubin_ra"], help="Center RA of Rubin/LSST FoV in deg (default: auto-detected or 35.70)")
    parser.add_argument("--rubin-dec", type=float, default=defaults["rubin_dec"], help="Center Dec of Rubin/LSST FoV in deg (default: auto-detected or -4.80)")
    parser.add_argument("--rubin-radius", type=float, default=defaults["rubin_radius"], help="Radius of Rubin/LSST FoV in deg (default: 1.75)")
    parser.add_argument("--no-rubin", action="store_true", help="Do not draw Rubin/LSST FoV circle")
    parser.add_argument("--no-pfs-hex", action="store_true", help="Do not draw PFS hexagonal FoVs")
    parser.add_argument("--title", default=defaults["title"], help="Custom plot title")

    args = parser.parse_args(remaining_argv)

    plot_sky_distribution_from_files(
        targets_dir=args.targets_dir,
        science_targets_file=args.science_targets,
        pointing_file=args.pointing_file,
        output_plot=args.output,
        show_rubin=(not args.no_rubin),
        rubin_center_ra=args.rubin_ra,
        rubin_center_dec=args.rubin_dec,
        rubin_radius=args.rubin_radius,
        show_pfs_hex=(not args.no_pfs_hex),
        title=args.title
    )

if __name__ == "__main__":
    main()
