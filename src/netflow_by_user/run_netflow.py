#!/usr/bin/env python3
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import argparse
import yaml
import numpy as np
import warnings
from erfa import ErfaWarning

# Ignore ERFA warnings about distance overridden
warnings.filterwarnings("ignore", category=ErfaWarning, message=".*distance overridden.*")

def add_dummy_targets_for_unassigned_near_bright_stars(bench, telescopes, pipe_config, exposures_data):
    import pandas as pd
    import numpy as np
    from astropy.table import Table
    from pfs.utils.coordinates.CoordTransp import CoordinateTransform as ctrans
    from pfs.utils.fiberids import FiberIds

    # Find the local Gaia catalog path
    gaia_path = pipe_config["inputs"].get("gaia_catalog", "cosmos/gaia.ecsv")
    if not os.path.exists(gaia_path):
        print(f"Warning: Gaia catalog not found at {gaia_path}. Skipping dummy target insertion.")
        return

    print("Loading Gaia catalog for dummy target verification...")
    t_gaia = Table.read(gaia_path, format="ascii.ecsv")
    df_gaia = t_gaia.to_pandas()
    df_gaia["magnitude"] = df_gaia["phot_g_mean_mag"]

    # Bright star thresholds
    bright_mag_limit = pipe_config.get("netflow", {}).get("bright_star_mag_limit", 12.0)
    bright_radius_arcmin = pipe_config.get("netflow", {}).get("bright_star_radius_arcmin", 1.5)
    radius_deg = bright_radius_arcmin / 60.0

    df_bright = df_gaia[df_gaia["magnitude"] <= bright_mag_limit]
    if len(df_bright) == 0:
        print("No bright stars found in catalog. Skipping dummy target insertion.")
        return

    targets_dir = pipe_config["outputs"]["targets_dir"]
    pointing_file = pipe_config["inputs"]["pointing_file"]
    if pointing_file is None:
        pointing_file = "optimized_pointings.ecsv"
    ppcList = Table.read(pointing_file, format="ascii.ecsv")
    ppc_codes = ppcList['ppc_code'].tolist()

    # Load fiber mapping
    import pfs.utils
    pfs_utils_dir = os.path.dirname(pfs.utils.__file__)
    path = os.path.join(pfs_utils_dir, "data", "fiberids")
    fibId = FiberIds(path)

    # Detailed report data list
    report_lines = []
    report_lines.append("Dummy Target Distance Improvement Report\n")
    report_lines.append("="*80 + "\n")

    for ivis, tel in enumerate(telescopes):
        ppc_code = ppc_codes[ivis] if ivis < len(ppc_codes) else f"EXP_{ivis+1}"
        sci_path = os.path.join(targets_dir, "science", f"{ppc_code}.ecsv")
        if not os.path.exists(sci_path):
            print(f"Warning: Science target file {sci_path} not found. Skipping...")
            continue

        # Find assigned cobras in this exposure
        assigned_cobra_ids = set()
        for item in exposures_data:
            if item['ppc_code'] == ppc_code:
                assigned_cobra_ids.add(item['cobraId'] - 1)

        # Identify healthy unassigned cobras
        unassigned_healthy_cobras = []
        for cidx in range(bench.cobras.nCobras):
            if bench.cobras.isGood[cidx] and cidx not in assigned_cobra_ids:
                unassigned_healthy_cobras.append(cidx)

        if not unassigned_healthy_cobras:
            print(f"No unassigned healthy cobras for {ppc_code}.")
            continue

        # Load and project science targets in this exposure for plotting
        try:
            t_sci_full = Table.read(sci_path, format="ascii.ecsv")
            if len(t_sci_full) > 0:
                sci_sky = np.array([t_sci_full['ra'].tolist(), t_sci_full['dec'].tolist()])
                sci_pfi_coords = ctrans(
                    xyin=sci_sky,
                    mode="sky_pfi",
                    pa=tel._posang,
                    cent=np.array([tel._ra, tel._dec]).reshape((2, 1)),
                    pm=np.zeros_like(sci_sky),
                    par=np.zeros(len(t_sci_full)),
                    time=tel._time,
                )
                sci_pfi = sci_pfi_coords[0, :] + 1j * sci_pfi_coords[1, :]
            else:
                sci_pfi = np.array([], dtype=complex)
        except Exception as e:
            print(f"Warning: Could not project science targets for plotting: {e}")
            sci_pfi = np.array([], dtype=complex)

        # Coordinates of unassigned cobras
        unassigned_centers = bench.cobras.centers[unassigned_healthy_cobras]
        unassigned_pfi = np.array([unassigned_centers.real, unassigned_centers.imag])

        # Convert to sky coordinates
        unassigned_sky = ctrans(
            xyin=unassigned_pfi,
            mode="pfi_sky",
            pa=tel._posang,
            cent=np.array([tel._ra, tel._dec]).reshape((2, 1)),
            time=tel._time,
            epoch=2016.0
        )
        unassigned_ra = unassigned_sky[0, :]
        unassigned_dec = unassigned_sky[1, :]

        # Find unassigned cobras near bright stars
        cobras_to_dummy = []
        cobras_to_dummy_star_pos = []
        cobras_to_dummy_star_sky = []
        cobras_to_dummy_star_mag = []

        cos_dec_tel = np.cos(np.radians(tel._dec))
        dist_to_center = np.hypot((df_bright['ra'] - tel._ra) * cos_dec_tel, df_bright['dec'] - tel._dec)
        df_bright_near = df_bright[dist_to_center < 0.8]

        if len(df_bright_near) > 0:
            bright_ra = df_bright_near['ra'].values
            bright_dec = df_bright_near['dec'].values
            bright_sky_coords = np.array([bright_ra, bright_dec])
            bright_pfi_coords = ctrans(
                xyin=bright_sky_coords,
                mode="sky_pfi",
                pa=tel._posang,
                cent=np.array([tel._ra, tel._dec]).reshape((2, 1)),
                pm=np.zeros_like(bright_sky_coords),
                par=np.zeros(len(bright_ra)),
                time=tel._time,
            )
            bright_pfi = bright_pfi_coords[0, :] + 1j * bright_pfi_coords[1, :]

            cos_dec_un = np.cos(np.radians(unassigned_dec))[:, np.newaxis]
            d_ra = (unassigned_ra[:, np.newaxis] - bright_ra[np.newaxis, :]) * cos_dec_un
            d_dec = unassigned_dec[:, np.newaxis] - bright_dec[np.newaxis, :]
            dist_deg = np.hypot(d_ra, d_dec)

            for i, cidx in enumerate(unassigned_healthy_cobras):
                near_star_indices = np.where(dist_deg[i] < radius_deg)[0]
                if len(near_star_indices) > 0:
                    closest_idx = near_star_indices[np.argmin(dist_deg[i][near_star_indices])]
                    cobras_to_dummy.append(cidx)
                    cobras_to_dummy_star_pos.append(bright_pfi[closest_idx])
                    cobras_to_dummy_star_sky.append((bright_ra[closest_idx], bright_dec[closest_idx]))
                    cobras_to_dummy_star_mag.append(df_bright_near['magnitude'].values[closest_idx])

        if not cobras_to_dummy:
            print(f"No unassigned healthy fibers near bright stars for {ppc_code}.")
            continue

        print(f"Found {len(cobras_to_dummy)} unassigned fibers near bright stars for {ppc_code}. Generating dummy targets...")
        report_lines.append(f"\nExposure/Pointing: {ppc_code}\n")
        report_lines.append("-" * 40 + "\n")

        # Initialize current fiber positions for all cobras
        fiber_pos = np.zeros(bench.cobras.nCobras, dtype=complex)
        for item in exposures_data:
            if item['ppc_code'] == ppc_code:
                fiber_pos[item['cobraId'] - 1] = item['pfi_X'] + 1j * item['pfi_Y']
        for cidx in range(bench.cobras.nCobras):
            if cidx not in assigned_cobra_ids:
                fiber_pos[cidx] = bench.cobras.centers[cidx]

        collision_distance = pipe_config.get("netflow", {}).get("collision_distance", 2.0)
        safe_collision_distance = collision_distance + 0.1

        new_dummy_rows = []

        for cidx, P_star, (star_ra, star_dec), star_mag in zip(
            cobras_to_dummy, cobras_to_dummy_star_pos, cobras_to_dummy_star_sky, cobras_to_dummy_star_mag
        ):
            C = bench.cobras.centers[cidx]
            r_max = bench.cobras.rMax[cidx]
            r_min = bench.cobras.rMin[cidx]
            neighbors = bench.getCobraNeighbors(cidx)

            best_P = None
            max_dist_to_star = -1.0

            # Grid search - restrict to r_max - 0.05 to avoid floating point patrol region warnings in validation
            radii = np.linspace(r_max - 0.05, r_min + 0.05, 10)
            v_star = P_star - C
            angle_opposite = np.angle(-v_star)
            angles = angle_opposite + np.linspace(-np.pi, np.pi, 72, endpoint=False)

            for r_test in radii:
                for angle in angles:
                    P_test = C + r_test * np.exp(1j * angle)
                    has_collision = False
                    for n_cidx in neighbors:
                        if np.abs(P_test - fiber_pos[n_cidx]) < safe_collision_distance:
                            has_collision = True
                            break
                    if not has_collision:
                        d_star = np.abs(P_test - P_star)
                        if d_star > max_dist_to_star:
                            max_dist_to_star = d_star
                            best_P = P_test

            if best_P is None:
                # Retry with exact collision distance
                for r_test in radii:
                    for angle in angles:
                        P_test = C + r_test * np.exp(1j * angle)
                        has_collision = False
                        for n_cidx in neighbors:
                            if np.abs(P_test - fiber_pos[n_cidx]) < collision_distance:
                                has_collision = True
                                break
                        if not has_collision:
                            d_star = np.abs(P_test - P_star)
                            if d_star > max_dist_to_star:
                                max_dist_to_star = d_star
                                best_P = P_test

            if best_P is None:
                print(f"  Warning: could not find collision-free position for cobra {cidx+1}. Fallback to cobra center.")
                best_P = C
                max_dist_to_star = np.abs(C - P_star)

            # Update fiber_pos
            fiber_pos[cidx] = best_P

            # Convert best_P to sky coordinates
            dummy_sky = ctrans(
                xyin=np.array([[best_P.real], [best_P.imag]]),
                mode="pfi_sky",
                pa=tel._posang,
                cent=np.array([tel._ra, tel._dec]).reshape((2, 1)),
                time=tel._time,
                epoch=2016.0
            )
            ra_dummy = dummy_sky[0, 0]
            dec_dummy = dummy_sky[1, 0]

            # Calculate separations
            init_dist_pfi = np.abs(C - P_star)
            final_dist_pfi = np.abs(best_P - P_star)

            # Sky separation (arcsec)
            init_dist_sky = np.hypot(
                (unassigned_ra[unassigned_healthy_cobras.index(cidx)] - star_ra) * np.cos(np.deg2rad(star_dec)),
                unassigned_dec[unassigned_healthy_cobras.index(cidx)] - star_dec
            ) * 3600.0
            final_dist_sky = np.hypot(
                (ra_dummy - star_ra) * np.cos(np.deg2rad(star_dec)),
                dec_dummy - star_dec
            ) * 3600.0

            # Find fiber ID (pass as list to avoid 0D array indexing error)
            try:
                fib = fibId.cobraIdToFiberId([cidx + 1])[0]
            except Exception:
                fib = "N/A"

            # Log details
            detail_log = (f"Cobra {cidx+1} (Fiber {fib}):\n"
                          f"  Bright Star: RA={star_ra:.5f}, Dec={star_dec:.5f}, Mag={star_mag:.2f}\n"
                          f"  Cobra Center (Initial): PFI ({C.real:.2f}, {C.imag:.2f})\n"
                          f"  Dummy Target (Final):   PFI ({best_P.real:.2f}, {best_P.imag:.2f})\n"
                          f"  Distance (PFI): {init_dist_pfi:.3f} mm -> {final_dist_pfi:.3f} mm (diff: {final_dist_pfi - init_dist_pfi:+.3f} mm)\n"
                          f"  Distance (Sky): {init_dist_sky:.1f} arcsec -> {final_dist_sky:.1f} arcsec (diff: {final_dist_sky - init_dist_sky:+.1f} arcsec)\n")
            print(detail_log)
            report_lines.append(detail_log + "\n")

            # Generate target local plot showing fibers, science targets, and gaia stars
            try:
                import matplotlib
                matplotlib.use('Agg')
                import matplotlib.pyplot as plt

                fig, ax = plt.subplots(figsize=(7, 7))

                # 1. Plot Cobra center & patrol boundaries
                ax.plot(C.real, C.imag, 'b+', markersize=10, markeredgewidth=2, label=f'Cobra Center ({C.real:.2f}, {C.imag:.2f})')
                circle_outer = plt.Circle((C.real, C.imag), r_max, color='blue', fill=False, linestyle='--', alpha=0.5, label='Patrol Boundaries')
                circle_inner = plt.Circle((C.real, C.imag), r_min, color='blue', fill=False, linestyle='--', alpha=0.3)
                ax.add_patch(circle_outer)
                ax.add_patch(circle_inner)

                # 2. Plot neighbors and their collision boundary
                for idx_n, n_cidx in enumerate(neighbors):
                    C_n = bench.cobras.centers[n_cidx]
                    f_pos = fiber_pos[n_cidx]
                    lbl_n = 'Neighbor Centers' if idx_n == 0 else ""
                    lbl_f = 'Neighbor Fibers' if idx_n == 0 else ""
                    ax.plot(C_n.real, C_n.imag, 'x', color='gray', markersize=6, label=lbl_n)
                    ax.plot(f_pos.real, f_pos.imag, 'o', color='orange', markersize=5, label=lbl_f)
                    circle_coll = plt.Circle((f_pos.real, f_pos.imag), collision_distance, color='red', fill=False, alpha=0.15)
                    ax.add_patch(circle_coll)

                # 3. Plot nearby science targets
                if len(sci_pfi) > 0:
                    dist_sci = np.abs(sci_pfi - C)
                    near_sci_idx = np.where((dist_sci < 10.0) & (t_sci_full['cobraId'] != cidx + 1))[0]
                    first_sci = True
                    for s_idx in near_sci_idx:
                        pt = sci_pfi[s_idx]
                        lbl = 'Nearby Science Targets' if first_sci else ""
                        ax.plot(pt.real, pt.imag, 'o', color='purple', markersize=5, label=lbl)
                        first_sci = False

                # 4. Plot other Gaia stars near C
                if len(bright_pfi) > 0:
                    dist_gaia = np.abs(bright_pfi - C)
                    near_gaia_idx = np.where((dist_gaia < 10.0) & (bright_pfi != P_star))[0]
                    first_gaia = True
                    for g_idx in near_gaia_idx:
                        pt = bright_pfi[g_idx]
                        lbl = 'Other Gaia Stars' if first_gaia else ""
                        ax.plot(pt.real, pt.imag, '*', color='yellow', markersize=8, markeredgecolor='black', label=lbl)
                        first_gaia = False

                # 5. Plot the target bright Gaia star to avoid
                ax.plot(P_star.real, P_star.imag, '*', color='gold', markersize=14, markeredgecolor='black', markeredgewidth=1.5,
                        label=f'Target Bright Star (G={star_mag:.1f})')

                # 6. Plot initial fiber position (Cobra Center since it was unassigned)
                ax.plot(C.real, C.imag, 'bo', markersize=6, label='Initial Fiber Position')
                # 7. Plot final dummy target position
                ax.plot(best_P.real, best_P.imag, 'go', markersize=8, label='Final Position (Dummy Target)')

                # Arrow indicating movement
                ax.annotate("", xy=(best_P.real, best_P.imag), xytext=(C.real, C.imag),
                            arrowprops=dict(arrowstyle="->", color="green", lw=1.5, ls="--"))

                # Formatting
                ax.set_xlim(C.real - 10.0, C.real + 10.0)
                ax.set_ylim(C.imag - 10.0, C.imag + 10.0)
                ax.set_aspect('equal')
                ax.grid(True, linestyle=':', alpha=0.5)
                ax.set_xlabel('PFI X (mm)')
                ax.set_ylabel('PFI Y (mm)')
                ax.set_title(f'PPC: {ppc_code} | Cobra {cidx+1} (Fiber {fib})\nBright Star Avoidance for Unassigned Fiber')
                ax.legend(loc='upper right', fontsize='x-small', framealpha=0.9)

                # Save plot to targets_dir/science/dummy_plots_{ppc_code}
                plot_dir = os.path.join(targets_dir, "science", f"dummy_plots_{ppc_code}")
                os.makedirs(plot_dir, exist_ok=True)
                plot_path = os.path.join(plot_dir, f"cobra_{cidx+1}.png")
                plt.savefig(plot_path, dpi=150, bbox_inches='tight')
                plt.close()
                print(f"  Generated local plot: {plot_path}")
            except Exception as e:
                print(f"  Warning: Failed to generate plot for Cobra {cidx+1}: {e}")

            new_dummy_rows.append({
                'ob_code': f"dummy_{cidx+1}_{ppc_code}",
                'obj_id': 999000000 + cidx,
                'ra': ra_dummy,
                'dec': dec_dummy,
                'exptime': 3600.0,
                'priority': 4,
                'resolution': 'L',
                'r2_hsc': np.nan,
                'i2_hsc': np.nan,
                'z_hsc': np.nan,
                'reference_arm': "",
                'g_hsc': np.nan,
                'cobraId': cidx + 1,
                'pfi_X': best_P.real,
                'pfi_Y': best_P.imag
            })

        if new_dummy_rows:
            t_sci = Table.read(sci_path, format="ascii.ecsv")
            orig_path = sci_path + ".orig"
            if not os.path.exists(orig_path):
                t_sci.write(orig_path, format="ascii.ecsv", overwrite=True)
                print(f"  Saved original science target file to {orig_path}")

            for row in new_dummy_rows:
                row_dict = {}
                for col in t_sci.colnames:
                    row_dict[col] = row.get(col, None)
                t_sci.add_row(row_dict)
            t_sci.write(sci_path, format="ascii.ecsv", overwrite=True)
            print(f"  Saved updated science target file with dummy targets to {sci_path}")

    # Write report file
    report_path = os.path.join(targets_dir, "science", "dummy_target_improvements.txt")
    with open(report_path, "w") as f_rep:
        f_rep.writelines(report_lines)
    print(f"Written detailed improvements report to {report_path}")

def main():
    # 1. Parse configuration file path
    parser = argparse.ArgumentParser(description="Run netflow fiber assignment.")
    parser.add_argument("-c", "--config", default="netflow_pipeline_config.yaml", help="Path to pipeline configuration YAML file.")
    parser.add_argument("--config-yaml", default="config.yaml", help="Path to config.yaml to override Gurobi and PFS parameters.")
    args, unknown = parser.parse_known_args()

    if not os.path.exists(args.config):
        print(f"Error: Configuration file {args.config} not found.")
        sys.exit(1)

    print(f"Loading pipeline configuration from {args.config}...")
    with open(args.config, "r") as f:
        pipe_config = yaml.safe_load(f)

    def merge_config_yaml(pipe_config, config_yaml_path):
        """Merges Gurobi and PFS parameters from config.yaml into the pipeline config."""
        if os.path.exists(config_yaml_path):
            print(f"Merging parameters from {config_yaml_path}...")
            try:
                with open(config_yaml_path, 'r') as f:
                    config_yaml = yaml.safe_load(f)
                if config_yaml:
                    # Merge gurobi parameters
                    if "gurobi" in config_yaml and "param" in config_yaml["gurobi"]:
                        if "gurobi" not in pipe_config:
                            pipe_config["gurobi"] = {}
                        pipe_config["gurobi"].update(config_yaml["gurobi"]["param"])
                        print("  Merged Gurobi parameters.")
                    
                    # Merge pfs parameters
                    if "pfs" in config_yaml:
                        if "pfs" not in pipe_config:
                            pipe_config["pfs"] = {}
                        pipe_config["pfs"].update(config_yaml["pfs"])
                        print("  Merged PFS parameters.")
            except Exception as e:
                print(f"Warning: Failed to merge config.yaml: {e}")
        else:
            print(f"Warning: config.yaml not found at {config_yaml_path}. Skipping merge.")
        return pipe_config

    # Merge config.yaml if present
    pipe_config = merge_config_yaml(pipe_config, args.config_yaml)

    # Set random seed
    np.random.seed(pipe_config["netflow"]["random_seed"])

    # Helper imports from refactored modules
    import netflow_io
    import netflow_instrument
    import netflow_solver
    import netflow_plot

    # 2. Load input targets
    catalog_path = pipe_config["inputs"]["catalog_dir"]
    def resolve_path(base_dir, path):
        if os.path.isabs(path):
            return path
        return os.path.join(base_dir, path)

    fscience_targets = resolve_path(catalog_path, pipe_config["inputs"]["science_targets"])
    fcal_stars = resolve_path(catalog_path, pipe_config["inputs"]["fluxstd_targets"])
    fsky_pos = resolve_path(catalog_path, pipe_config["inputs"]["sky_targets"])

    tgt = netflow_io.load_all_targets(
        fscience_targets, 
        fcal_stars, 
        fsky_pos, 
        fluxstd_mag_min=pipe_config.get("netflow", {}).get("fluxstd", {}).get("mag_min", 17.0),
        fluxstd_mag_max=pipe_config.get("netflow", {}).get("fluxstd", {}).get("mag_max", 19.0)
    )

    # 3. Setup Instrument (Bench) and Telescopes (Pointing Center)
    black_dot_radius_margin = pipe_config.get("pfs", {}).get("black_dot_radius_margin", 1.65)
    bench = netflow_instrument.getBench(black_dot_radius_margin)

    nvisit = pipe_config["netflow"]["nvisit"]
    posang = pipe_config["netflow"].get("posang", 0.0)
    otime = pipe_config["obstime"]
    pointing_file = pipe_config["inputs"]["pointing_file"]

    if pointing_file is None:
        print("\npointing_file is null. Running automatic FoV optimization with guide star constraints...")
        num_fovs = pipe_config["netflow"].get("num_fields", 1)
        max_priority = pipe_config["netflow"].get("max_priority", 2)
        min_stars_per_cam = pipe_config["netflow"].get("min_stars_per_cam", 2)
        min_cams_with_stars = pipe_config["netflow"].get("min_cams_with_stars", 6)
        
        # 視野最適化関数の読み込み
        from optimize_hex_fov_with_guidestars import optimize_fovs_with_guidestars, plot_optimized_fovs
        import pandas as pd
        from astropy.table import Table
        
        print(f"Reading targets from {fscience_targets} for optimization...")
        df_tgt = pd.read_csv(fscience_targets)
        
        # 列名のクリーニング (optimize_hex_fov.py より移植)
        col_mapping = {}
        for col in ['ra', 'dec', 'priority', 'obj_id']:
            found = False
            for c in df_tgt.columns:
                if c.lower() == col or c.lower().replace('.', '').replace(' ', '_') == col:
                    col_mapping[col] = c
                    found = True
                    break
            if not found:
                if col == 'obj_id' and 'ob_code' in df_tgt.columns:
                    col_mapping['obj_id'] = 'ob_code'
                elif col == 'obj_id' and 'ID' in df_tgt.columns:
                    col_mapping['obj_id'] = 'ID'
                else:
                    raise KeyError(f"Could not find equivalent column for '{col}' in input CSV.")
                    
        df_clean = pd.DataFrame({
            'ra': df_tgt[col_mapping['ra']],
            'dec': df_tgt[col_mapping['dec']],
            'priority': df_tgt[col_mapping['priority']],
            'obj_id': df_tgt[col_mapping['obj_id']]
        })
        
        print(f"Filtering targets with priority <= {max_priority}...")
        df_filtered = df_clean[df_clean['priority'] <= max_priority].copy()
        print(f"Filtered target count: {len(df_filtered)}")
        
        if len(df_filtered) == 0:
            raise ValueError(f"No targets found with priority <= {max_priority} for optimization.")
            
        fgaia_catalog = pipe_config["inputs"].get("gaia_catalog", "cosmos/gaia.ecsv")
        print(f"Reading Gaia catalog from: {fgaia_catalog} for optimization...")
        t_gaia = Table.read(fgaia_catalog, format="ascii.ecsv")
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
        
        gs_config = pipe_config["netflow"].get("guidestars", {})
        gs_mag_min = gs_config.get("mag_min", 17.0)
        gs_mag_max = gs_config.get("mag_max", 21.5)
        print(f"Optimizing for {num_fovs} fields with guide star magnitude range: {gs_mag_min} - {gs_mag_max} mag...")
        bright_mag_limit = pipe_config.get("netflow", {}).get("bright_star_mag_limit", 12.0)
        bright_radius_arcmin = pipe_config.get("netflow", {}).get("bright_star_radius_arcmin", 1.5)
        
        pointings, covered = optimize_fovs_with_guidestars(
            df_filtered, df_gaia, otime, num_fovs=num_fovs,
            min_stars_per_cam=min_stars_per_cam, min_cams_with_stars=min_cams_with_stars,
            pa_step=5.0, min_mag=gs_mag_min, max_mag=gs_mag_max,
            bench=bench,
            bright_star_mag_limit=bright_mag_limit,
            bright_star_radius_arcmin=bright_radius_arcmin
        )
        
        # 結果を ECSV に保存する
        opt_table = Table(
            names=('ppc_code', 'ppc_ra', 'ppc_dec', 'ppc_pa', 'covered_count'),
            dtype=('S30', 'f8', 'f8', 'f8', 'i4')
        )
        for p in pointings:
            opt_table.add_row((p['ppc_code'], p['ppc_ra'], p['ppc_dec'], p['ppc_pa'], p['covered_count']))
            
        pointing_file = "optimized_pointings.ecsv"
        opt_table.write(pointing_file, format="ascii.ecsv", overwrite=True)
        print(f"FoV optimization completed. Pointings saved to {pointing_file}")
        
        # Generate FoV optimization plot
        fov_plot_file = pipe_config["outputs"].get("fov_plot_file", "fov_coverage.png")
        plot_optimized_fovs(pointings, covered, df_filtered, max_priority, fov_plot_file, otime)

    telescopes = netflow_instrument.getPointingCenter(pointing_file, nvisit, posang, otime)

    # 4. Solve the fiber assignment netflow problem
    res, tpos, exposures_data, all_classes, stats_per_exp = netflow_solver.solve_assignment(
        bench, tgt, telescopes, pipe_config
    )

    # 5. Output text summary / stats to standard output
    print("\n--- Observation Stats per Exposure ---")
    def sort_cls(c):
        if c.startswith('sci_P'):
            return (0, int(c.split('_P')[1]))
        elif c == 'cal': return (1, 0)
        elif c == 'sky': return (2, 0)
        return (3, c)

    sorted_classes = sorted(list(all_classes), key=sort_cls)
    header = ["Exp", "Total"] + sorted_classes
    row_format = "{:>5} | {:>7} | " + " | ".join(["{:>8}"] * len(sorted_classes))
    print(row_format.format(*header))
    print("-" * (18 + 11 * len(sorted_classes)))
    for i, total, tdict in stats_per_exp:
        row = [i, total] + [tdict[cls] for cls in sorted_classes]
        print(row_format.format(*row))
    print("-" * (18 + 11 * len(sorted_classes)))

    # 6. Match targets and save ECSV files under target directories
    netflow_io.save_targets_ecsv(
        exposures_data, fscience_targets, fcal_stars, fsky_pos, pipe_config["outputs"]["targets_dir"]
    )

    # 6b. Add dummy targets to unassigned fibers near bright stars
    add_dummy_targets_for_unassigned_near_bright_stars(
        bench, telescopes, pipe_config, exposures_data
    )

    # 7. Plot sky distribution
    netflow_plot.plot_sky_distribution(res, tgt, pipe_config["outputs"]["plot_file"])

    # 8. Generate pfsDesign FITS files and OPE templates
    from make_pfs_design import generate_pfs_designs
    generate_pfs_designs(
        pipeline_config_path=args.config,
        pipeline_config=pipe_config
    )
if __name__ == '__main__':
    main()
