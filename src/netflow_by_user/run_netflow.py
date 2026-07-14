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
