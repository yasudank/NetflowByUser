#!/usr/bin/env python3
"""
Plot the overlap between science target positions and PFS hexagonal FoV pointings.
"""
import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from astropy.table import Table

# Regular hexagon parameters
R_HEX = 0.69  # Outer radius of PFS FoV in degrees

def get_vertices(ra, dec, r, pa_deg):
    """Get the 6 vertices of a rotated hexagon centered at (ra, dec) with radius r and rotation pa_deg."""
    cos_dec = np.cos(np.radians(dec))
    
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

def compute_pointing_mask(ra, dec, center_ra, center_dec, pa_deg, r=R_HEX):
    """Compute boolean mask indicating whether (ra, dec) points fall inside the rotated hexagon."""
    cos_dec = np.cos(np.radians(center_dec))
    dx = (ra - center_ra) * cos_dec
    dy = dec - center_dec
    
    rad = np.radians(pa_deg)
    cos_pa = np.cos(rad)
    sin_pa = np.sin(rad)
    
    # Rotate clockwise to align with hexagon axes
    dx_rot = dx * cos_pa + dy * sin_pa
    dy_rot = -dx * sin_pa + dy * cos_pa
    
    cond1 = np.abs(dy_rot) <= r * np.sqrt(3) / 2.0
    cond2 = np.abs(dy_rot) + np.sqrt(3) * np.abs(dx_rot) <= np.sqrt(3) * r
    return cond1 & cond2

def load_pointings(pointing_path):
    """Load pointings from ECSV or CSV."""
    if pointing_path.endswith('.ecsv'):
        t = Table.read(pointing_path, format='ascii.ecsv')
        df = t.to_pandas()
    else:
        df = pd.read_csv(pointing_path)
        
    col_map = {}
    for c in df.columns:
        cl = c.lower()
        if cl in ['ppc_ra', 'ra']: col_map['ra'] = c
        elif cl in ['ppc_dec', 'dec']: col_map['dec'] = c
        elif cl in ['ppc_pa', 'pa']: col_map['pa'] = c
        elif cl in ['ppc_code', 'pointing_id', 'code', 'id']: col_map['code'] = c
        
    if 'ra' not in col_map or 'dec' not in col_map:
        raise ValueError(f"Could not find RA/Dec columns in pointing file: {pointing_path}")
        
    pointings = []
    for i, row in df.iterrows():
        code = str(row[col_map['code']]) if 'code' in col_map else f"FoV_{i+1}"
        ra = float(row[col_map['ra']])
        dec = float(row[col_map['dec']])
        pa = float(row[col_map['pa']]) if 'pa' in col_map else 0.0
        pointings.append({
            'code': code,
            'ra': ra,
            'dec': dec,
            'pa': pa,
        })
    return pointings

def load_targets(target_path, max_priority=None):
    """Load targets from CSV or ECSV and optionally filter by max_priority."""
    if target_path.endswith('.ecsv'):
        t = Table.read(target_path, format='ascii.ecsv')
        df = t.to_pandas()
    else:
        df = pd.read_csv(target_path)
        
    col_map = {}
    for c in df.columns:
        cl = c.lower().replace('.', '').replace(' ', '_')
        if cl in ['ra']: col_map['ra'] = c
        elif cl in ['dec']: col_map['dec'] = c
        elif cl in ['priority', 'prio']: col_map['priority'] = c
        elif cl in ['obj_id', 'targetid', 'id']: col_map['obj_id'] = c
        
    if 'ra' not in col_map or 'dec' not in col_map:
        raise ValueError(f"Could not find RA/Dec columns in target file: {target_path}")
        
    res_df = pd.DataFrame()
    res_df['ra'] = df[col_map['ra']].astype(float)
    res_df['dec'] = df[col_map['dec']].astype(float)
    
    if 'priority' in col_map:
        res_df['priority'] = df[col_map['priority']].astype(int)
    else:
        res_df['priority'] = 0
        
    if 'obj_id' in col_map:
        res_df['obj_id'] = df[col_map['obj_id']]
    else:
        res_df['obj_id'] = np.arange(len(df))
        
    if max_priority is not None:
        res_df = res_df[res_df['priority'] <= max_priority].copy().reset_index(drop=True)
        
    return res_df

def plot_overlap(pointings, df_targets, output_plot, color_by='pointing', title=None):
    """Plot targets and pointing hexagons."""
    n_targets = len(df_targets)
    ra_arr = df_targets['ra'].values
    dec_arr = df_targets['dec'].values
    
    # Calculate coverage masks for each pointing
    coverage_masks = []
    for p in pointings:
        mask = compute_pointing_mask(ra_arr, dec_arr, p['ra'], p['dec'], p['pa'], R_HEX)
        coverage_masks.append(mask)
        p['covered_count'] = int(mask.sum())
        
    coverage_matrix = np.array(coverage_masks) # Shape: (n_pointings, n_targets)
    total_covered_mask = coverage_matrix.any(axis=0) if len(pointings) > 0 else np.zeros(n_targets, dtype=bool)
    n_covered_total = int(total_covered_mask.sum())
    overlap_counts = coverage_matrix.sum(axis=0) if len(pointings) > 0 else np.zeros(n_targets, dtype=int)
    
    print("\n================ Coverage Statistics ================")
    print(f"Total Targets Evaluated: {n_targets}")
    print(f"Total Unique Targets Covered: {n_covered_total} ({100.0 * n_covered_total / max(1, n_targets):.1f}%)")
    print(f"Uncovered Targets: {n_targets - n_covered_total} ({100.0 * (n_targets - n_covered_total) / max(1, n_targets):.1f}%)")
    print("-----------------------------------------------------")
    for i, p in enumerate(pointings):
        print(f"  [{i+1}] {p['code']:<20}: RA={p['ra']:9.4f}, Dec={p['dec']:9.4f}, PA={p['pa']:5.1f}° | Covered: {p['covered_count']:5d}")
    print("-----------------------------------------------------")
    for k in range(len(pointings) + 1):
        c_k = int((overlap_counts == k).sum())
        if c_k > 0:
            print(f"  Targets covered by exactly {k} FoV(s): {c_k:6d} ({100.0 * c_k / max(1, n_targets):.1f}%)")
    print("=====================================================\n")
    
    plt.style.use('default')
    fig, ax = plt.subplots(figsize=(11, 10), dpi=150)
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')
    
    # Color palettes
    fov_palette = [
        '#2563eb', '#dc2626', '#16a34a', '#d97706', '#9333ea',
        '#0891b2', '#ea580c', '#4f46e5', '#059669', '#be123c'
    ]
    prio_palette = {
        0: '#dc2626', 1: '#ea580c', 2: '#eab308', 3: '#84cc16',
        4: '#10b981', 5: '#06b6d4', 6: '#3b82f6', 7: '#6366f1',
        8: '#a855f7', 9: '#64748b'
    }
    
    # 1. Uncovered targets (background)
    uncovered = df_targets[~total_covered_mask]
    if len(uncovered) > 0:
        ax.scatter(uncovered['ra'], uncovered['dec'],
                   color='#cbd5e1', alpha=0.35, s=12, label=f"Uncovered ({len(uncovered)})", marker='.', zorder=1)
                   
    # 2. Covered targets depending on color_by mode
    if color_by == 'pointing':
        # Mark each pointing's covered targets
        for idx, p in enumerate(pointings):
            c_color = fov_palette[idx % len(fov_palette)]
            mask = coverage_masks[idx]
            sub_df = df_targets[mask]
            ax.scatter(sub_df['ra'], sub_df['dec'],
                       color=c_color, alpha=0.6, s=18,
                       label=f"{p['code']} ({len(sub_df)})", zorder=3)
    elif color_by == 'priority':
        # Group by priority
        covered_df = df_targets[total_covered_mask]
        for prio, group in covered_df.groupby('priority'):
            c_color = prio_palette.get(prio, '#334155')
            ax.scatter(group['ra'], group['dec'],
                       color=c_color, alpha=0.7, s=20,
                       label=f"Priority {prio} ({len(group)})", zorder=3)
    elif color_by == 'overlap':
        # Group by number of covering FoVs
        for k in range(1, len(pointings) + 1):
            k_mask = (overlap_counts == k)
            if not k_mask.any(): continue
            sub_df = df_targets[k_mask]
            color = fov_palette[(k - 1) % len(fov_palette)]
            ax.scatter(sub_df['ra'], sub_df['dec'],
                       color=color, alpha=0.7, s=20,
                       label=f"Covered by {k} FoV(s) ({len(sub_df)})", zorder=3)
                       
    # 3. Hexagon Footprints
    for idx, p in enumerate(pointings):
        ra_c, dec_c, pa = p['ra'], p['dec'], p['pa']
        v_ra, v_dec = get_vertices(ra_c, dec_c, R_HEX, pa)
        line_color = fov_palette[idx % len(fov_palette)] if color_by == 'pointing' else '#0284c7'
        
        # Hexagon fill
        poly = patches.Polygon(
            np.column_stack((v_ra, v_dec)),
            closed=True,
            facecolor=line_color,
            alpha=0.08,
            zorder=2
        )
        ax.add_patch(poly)
        
        # Hexagon outline
        ax.plot(v_ra, v_dec, color=line_color, linewidth=2.0, linestyle='-', zorder=4)
        
        # Center marker & text
        ax.scatter(ra_c, dec_c, color=line_color, s=35, marker='o', zorder=5)
        ax.text(
            ra_c, dec_c + 0.035, f"{p['code']}\nPA={pa:.1f}°",
            color='#0f172a', fontsize=9, fontweight='bold',
            ha='center', va='bottom', zorder=6,
            bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.75, edgecolor='none')
        )
        
    ax.set_xlabel('RA (deg)', fontsize=12)
    ax.set_ylabel('Dec (deg)', fontsize=12)
    
    if title is None:
        title = f"PFS Pointings & Target Overlap ({len(pointings)} FoVs, {n_covered_total}/{n_targets} Covered: {100.0*n_covered_total/max(1,n_targets):.1f}%)"
    ax.set_title(title, fontsize=13, fontweight='bold', pad=12)
    
    ax.grid(True, color='#e2e8f0', linestyle='--', alpha=0.6)
    ax.set_aspect('equal')
    ax.invert_xaxis()
    
    plt.legend(loc='upper right', bbox_to_anchor=(1.28, 1.0), borderaxespad=0.0, frameon=True, fontsize=9)
    plt.tight_layout()
    plt.savefig(output_plot, bbox_inches='tight', dpi=150, facecolor=fig.get_facecolor(), edgecolor='none')
    plt.close()
    print(f"Saved overlap plot to {output_plot}")

def main():
    # 0. Pre-parse config if provided
    conf_parser = argparse.ArgumentParser(add_help=False)
    conf_parser.add_argument("--config", "-c", help="Path to YAML configuration file (e.g. netflow_pipeline_config.yaml)")
    args_conf, remaining_argv = conf_parser.parse_known_args()
    
    defaults = {
        "pointings": None,
        "targets": None,
        "output": "pointing_target_overlap.png",
        "max_priority": None,
        "color_by": "pointing",
        "title": None,
    }
    
    if args_conf.config and os.path.exists(args_conf.config):
        import yaml
        print(f"Loading defaults from config: {args_conf.config}")
        with open(args_conf.config, "r") as f:
            cfg = yaml.safe_load(f)
        if cfg:
            inputs = cfg.get("inputs", {})
            cat_dir = inputs.get("catalog_dir", ".")
            pointing_file = inputs.get("pointing_file")
            if pointing_file:
                defaults["pointings"] = pointing_file
                
            sci_file = inputs.get("science_targets")
            if sci_file:
                defaults["targets"] = os.path.join(cat_dir, sci_file) if not os.path.isabs(sci_file) else sci_file
                
            netflow_cfg = cfg.get("netflow", {})
            if "max_priority" in netflow_cfg:
                defaults["max_priority"] = netflow_cfg["max_priority"]
                
            outputs = cfg.get("outputs", {})
            if "fov_plot_file" in outputs:
                defaults["output"] = outputs["fov_plot_file"]

    parser = argparse.ArgumentParser(
        description="Plot overlap between science targets and PFS hexagonal pointing fields (FoVs)."
    )
    parser.add_argument("--config", "-c", help="Path to YAML configuration file")
    parser.add_argument("--pointings", "-p", required=(defaults["pointings"] is None), default=defaults["pointings"],
                        help="Path to pointings ECSV or CSV file (must contain RA, Dec, and optional PA)")
    parser.add_argument("--targets", "-t", required=(defaults["targets"] is None), default=defaults["targets"],
                        help="Path to science targets CSV or ECSV file (must contain RA, Dec, and optional priority)")
    parser.add_argument("--output", "-o", default=defaults["output"], help="Path for output PNG plot (default: pointing_target_overlap.png)")
    parser.add_argument("--max-priority", type=int, default=defaults["max_priority"],
                        help="Filter targets with priority <= max_priority")
    parser.add_argument("--color-by", choices=['pointing', 'priority', 'overlap'], default=defaults["color_by"],
                        help="Target coloring scheme: 'pointing' (by FoV), 'priority' (by priority level), or 'overlap' (by coverage multiplicity)")
    parser.add_argument("--title", default=defaults["title"], help="Custom title for the plot")
    
    args = parser.parse_args(remaining_argv)
    
    print(f"Reading pointings from {args.pointings}...")
    pointings = load_pointings(args.pointings)
    print(f"Loaded {len(pointings)} pointings.")
    
    print(f"Reading targets from {args.targets}...")
    df_targets = load_targets(args.targets, max_priority=args.max_priority)
    print(f"Loaded {len(df_targets)} targets" + (f" with priority <= {args.max_priority}" if args.max_priority is not None else "") + ".")
    
    plot_overlap(pointings, df_targets, args.output, color_by=args.color_by, title=args.title)

if __name__ == "__main__":
    main()

