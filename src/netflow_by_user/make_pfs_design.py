#!/usr/bin/env python3
import os
import sys
import argparse
import numpy as np
import pandas as pd
from astropy.table import Table
from astropy.coordinates import Angle
from astropy import units as u
from astropy.time import Time
import yaml
from collections import defaultdict

# Astropy and PFS imports
import pfs.instdata
pfs_instdata_dir = os.path.dirname(pfs.instdata.__file__)
os.environ["PFS_INSTDATA_DIR"] = pfs_instdata_dir

import ets_fiber_assigner.netflow as nf
import pfs_design_tool.pointing_utils.designutils as designutils
import pfs_design_tool.pointing_utils.nfutils as nfutils
import pfs_design_tool.pointing_utils.dbutils as dbutils

def load_config_yaml(config_path="config.yaml"):
    if os.path.exists(config_path):
        with open(config_path, 'r') as f:
            return yaml.safe_load(f)
    return {}

def load_config_toml(config_path):
    if config_path and os.path.exists(config_path):
        try:
            import tomllib
            with open(config_path, "rb") as f:
                return tomllib.load(f)
        except Exception as e:
            print(f"Warning: Could not load TOML config {config_path}: {e}")
    return {}

def prepare_science_df(df, config_yaml):
    # Add expected columns if missing
    df = df.copy()
    
    # Proposal ID
    if "proposal_id" not in df.columns:
        df["proposal_id"] = config_yaml.get("proposal_id", "S26A-OT02")
        
    # Input Catalog ID (Always overwrite/use input_catalog.science.id for science targets)
    df["input_catalog_id"] = config_yaml.get("input_catalog", {}).get("science", {}).get("id", 25040)
        
    # Target Type ID (1 = Science)
    if "target_type_id" not in df.columns:
        df["target_type_id"] = 1
        
    # Epoch
    if "epoch" not in df.columns:
        df["epoch"] = "J2000.0"
        
    # pmra, pmdec, parallax
    for col, def_val in [("pmra", 0.0), ("pmdec", 0.0), ("parallax", 1e-7)]:
        if col not in df.columns:
            df[col] = def_val
            
    # Effective exptime
    if "effective_exptime" not in df.columns:
        df["effective_exptime"] = df["exptime"] if "exptime" in df.columns else 3600.0
        
    # qa_reference_arm
    if "qa_reference_arm" not in df.columns:
        df["qa_reference_arm"] = df["reference_arm"] if ("reference_arm" in df.columns and not df["reference_arm"].isna().all()) else "r"
        
    # Flux mappings
    flux_map = {
        "g": "g_hsc",
        "r": "r2_hsc",
        "i": "i2_hsc",
        "z": "z_hsc",
        "y": "y_hsc"
    }
    
    for band, orig_col in flux_map.items():
        # filter names
        filt_col = f"filter_{band}"
        if filt_col not in df.columns:
            if orig_col in df.columns and not df[orig_col].isna().all():
                df[filt_col] = orig_col
            else:
                df[filt_col] = "None"
                
        # total flux
        tot_col = f"total_flux_{band}"
        if tot_col not in df.columns:
            if orig_col in df.columns:
                df[tot_col] = df[orig_col].astype(float)
            else:
                df[tot_col] = np.nan
                
        # total flux error
        err_col = f"total_flux_error_{band}"
        if err_col not in df.columns:
            df[err_col] = np.nan
            
    # Cobra index (0-based)
    if "cidx" not in df.columns and "cobraId" in df.columns:
        df["cidx"] = df["cobraId"] - 1
        
    # netflow ID
    df["netflow_id"] = df["obj_id"].astype(str) + "_" + df["input_catalog_id"].astype(str)
    
    return df

def prepare_fluxstd_df(df, config_yaml):
    df = df.copy()
    
    if "input_catalog_id" not in df.columns:
        df["input_catalog_id"] = config_yaml.get("input_catalog", {}).get("fluxstd", {}).get("id", 3012)
        
    if "target_type_id" not in df.columns:
        df["target_type_id"] = 3
        
    for col, def_val in [("pmra", 0.0), ("pmdec", 0.0), ("parallax", 1e-7), ("prob_f_star", 0.0)]:
        if col not in df.columns:
            df[col] = def_val
            
    if "epoch" not in df.columns:
        df["epoch"] = "J2000.0"
        
    for band in ["g", "r", "i", "z", "y"]:
        filt_col = f"filter_{band}"
        if filt_col not in df.columns:
            df[filt_col] = f"{band}_ps1"
            
        flux_col = f"psf_flux_{band}"
        if flux_col not in df.columns:
            df[flux_col] = np.nan
            
        err_col = f"psf_flux_error_{band}"
        if err_col not in df.columns:
            df[err_col] = np.nan
            
    if "cidx" not in df.columns and "cobraId" in df.columns:
        df["cidx"] = df["cobraId"] - 1
        
    df["netflow_id"] = df["obj_id"].astype(str) + "_" + df["input_catalog_id"].astype(str)
    
    return df

def prepare_sky_df(df, config_yaml):
    df = df.copy()
    
    if "input_catalog_id" not in df.columns:
        df["input_catalog_id"] = config_yaml.get("input_catalog", {}).get("sky_ps1", {}).get("id", 1007)
        
    if "target_type_id" not in df.columns:
        df["target_type_id"] = 2
        
    if "cidx" not in df.columns and "cobraId" in df.columns:
        df["cidx"] = df["cobraId"] - 1
        
    df["netflow_id"] = df["obj_id"].astype(str) + "_" + df["input_catalog_id"].astype(str)
    
    return df

def generate_pfs_designs(
    pointing_file="optimized_pointings.ecsv",
    targets_dir="targets",
    outdir="design",
    config_toml_path=None,
    gaia_catalog="cosmos/gaia.ecsv",
    obstime="2026-05-09T06:00:00Z",
    config_yaml_path="config.yaml",
    ope_template=None,
    ope_outdir="ope",
    exptime_per_frame=900.0,
    n_frames=4,
    pipeline_config_path=None,
    pipeline_config=None,
):
    """Generate pfsDesign FITS files from fiber assignment results.
    
    Parameters
    ----------
    pointing_file : str
        Path to the pointing list ECSV file.
    targets_dir : str
        Directory where target ECSV files are saved (science/, fluxstd/, sky/).
    outdir : str
        Output directory for generated pfsDesign FITS files.
    config_toml_path : str or None
        Path to config.toml for GaiaDB configuration.
    gaia_catalog : str or None
        Path to local Gaia ECSV catalog for guide stars.
    obstime : str
        Observing time in UTC.
    config_yaml_path : str
        Path to config.yaml.
    ope_template : str or None
        Path to OPE file template (.ope). If None, OPE generation is skipped.
    ope_outdir : str
        Output directory for generated OPE files.
    exptime_per_frame : float
        Exposure time per sub-frame in seconds (default: 900.0).
    n_frames : int
        Number of sub-frames per pointing (default: 4).
    pipeline_config_path : str or None
        Path to pipeline configuration YAML file.
    """
    # Load config files
    if pipeline_config:
        p_cfg = pipeline_config
    elif pipeline_config_path and os.path.exists(pipeline_config_path):
        print(f"Loading pipeline configuration from {pipeline_config_path}...")
        with open(pipeline_config_path, "r") as f:
            p_cfg = yaml.safe_load(f)
            
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
            
        p_cfg = merge_config_yaml(p_cfg, config_yaml_path)
    else:
        p_cfg = None

    if p_cfg:
        pointing_file = p_cfg.get("inputs", {}).get("pointing_file", pointing_file)
        if pointing_file is None:
            pointing_file = "optimized_pointings.ecsv"
        targets_dir = p_cfg.get("outputs", {}).get("targets_dir", targets_dir)
        # Always output to "design" and "ope" directories as forced by the pfs_obsproc_planning library
        outdir = "design"
        gaia_catalog = p_cfg.get("inputs", {}).get("gaia_catalog", gaia_catalog)
        obstime = p_cfg.get("obstime", obstime)
        ope_template = p_cfg.get("pfs_design", {}).get("ope_template", ope_template)
        ope_outdir = "ope"
        exptime_per_frame = p_cfg.get("pfs_design", {}).get("exptime_per_frame", exptime_per_frame)
        n_frames = p_cfg.get("pfs_design", {}).get("n_frames", n_frames)
        config_toml_path = p_cfg.get("inputs", {}).get("gaiadb_config", config_toml_path)
        
        config_yaml = {
            "proposal_id": p_cfg.get("proposal_id", "S26A-OT02"),
            "input_catalog": {
                "science": {"id": p_cfg.get("pfs_design", {}).get("science_catalog_id", 25040)},
                "fluxstd": {"id": 3012},
                "sky_ps1": {"id": 1007}
            }
        }
    else:
        config_yaml = load_config_yaml(config_yaml_path)
        
    config_toml = load_config_toml(config_toml_path)
    
    # 1. Read pointing list
    print(f"\n=== Generating pfsDesign files ===")
    print(f"Reading pointing list from {pointing_file}...")
    if not os.path.exists(pointing_file):
        print(f"Error: Pointing file {pointing_file} does not exist.")
        return
    pointings = Table.read(pointing_file, format="ascii.ecsv")
    
    # Create output directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        print(f"Created output directory: {outdir}")

    # List to collect per-pointing info for OPE file generation
    ope_info_rows = []
    # List to collect summary info for validation
    summary_rows = []

    # Iterate through each pointing center
    for r_ptg in pointings:
        ppc_code = r_ptg["ppc_code"]
        ra_tel = r_ptg["ppc_ra"]
        dec_tel = r_ptg["ppc_dec"]
        pa_tel = r_ptg["ppc_pa"]
        
        print(f"\nProcessing Pointing: {ppc_code} (RA: {ra_tel:.5f}, Dec: {dec_tel:.5f}, PA: {pa_tel:.1f})")
        
        # Load targets ECSV
        sci_path = os.path.join(targets_dir, "science", f"{ppc_code}.ecsv")
        cal_path = os.path.join(targets_dir, "fluxstd", f"{ppc_code}.ecsv")
        sky_path = os.path.join(targets_dir, "sky", f"{ppc_code}.ecsv")
        
        if not (os.path.exists(sci_path) and os.path.exists(cal_path) and os.path.exists(sky_path)):
            print(f"  Warning: Target files for {ppc_code} not found. Skipping...")
            continue
            
        df_sci = Table.read(sci_path, format="ascii.ecsv").to_pandas()
        df_flux = Table.read(cal_path, format="ascii.ecsv").to_pandas()
        df_sky = Table.read(sky_path, format="ascii.ecsv").to_pandas()
        
        # Prepare dataframes
        df_sci_prep = prepare_science_df(df_sci, config_yaml)
        df_flux_prep = prepare_fluxstd_df(df_flux, config_yaml)
        df_sky_prep = prepare_sky_df(df_sky, config_yaml)
        
        # Concatenate targets & register objects
        def vis_generator(df):
            vis_ = {}
            for idx, r in df.iterrows():
                vis_[idx] = int(r["cidx"])
            return vis_
            
        # We need apply_nir_flag check
        apply_nir = True
        
        target1 = nfutils.register_objects(df_sci_prep, target_class="sci", apply_nir_flag=apply_nir)
        vis1 = vis_generator(df_sci_prep)
        
        target2 = nfutils.register_objects(df_flux_prep, target_class="cal")
        vis2 = vis_generator(df_flux_prep)
        
        target3 = nfutils.register_objects(df_sky_prep, target_class="sky")
        vis3 = vis_generator(df_sky_prep)
        
        targets = target1 + target2 + target3
        
        vis2_update = {k + len(vis1): v for k, v in vis2.items()}
        vis3_update = {k + len(vis1) + len(vis2): v for k, v in vis3.items()}
        
        vis = {}
        vis.update(vis1)
        vis.update(vis2_update)
        vis.update(vis3_update)
        
        # Target class mapping
        target_class_dict = {
            **{f"sci_P{i}": 1 for i in range(10001)},
            "sky": 2,
            "cal": 3,
        }
        
        # Telescopes
        tele = nf.Telescope(ra_tel, dec_tel, pa_tel, obstime)
        tpos = tele.get_fp_positions(targets)
        
        # Spectrograph arm
        resolution = df_sci['resolution'].values[0] if 'resolution' in df_sci.columns else 'L'
        arm_ = 'brn' if resolution == 'L' else 'bmn'
        
        # Generate pfsDesign
        print("  Generating pfsDesign object...")
        design = designutils.generate_pfs_design(
            df_sci_prep,
            df_flux_prep,
            df_sky_prep,
            vis,
            tpos,
            tele,
            targets,
            target_class_dict,
            bench=None,
            arms=arm_,
            design_name=ppc_code,
            obs_time=obstime
        )
        
        # Query guide stars (from local catalog or database)
        guidestars = None
        if gaia_catalog and os.path.exists(gaia_catalog):
            try:
                print(f"  Selecting guide stars from local catalog {gaia_catalog}...")
                # Load ECSV table and convert to pandas DataFrame
                t_gaia = Table.read(gaia_catalog, format="ascii.ecsv")
                df_gaia = t_gaia.to_pandas()
                
                # Add columns magnitude and color to prevent KeyError in ets_pointing
                df_gaia["magnitude"] = df_gaia["phot_g_mean_mag"]
                df_gaia["color"] = df_gaia["bp_rp"]
                
                # Pre-populate all columns expected by get_gs_flag to bypass ets_pointing bugs if missing
                if "catalog" not in df_gaia.columns:
                    df_gaia["catalog"] = "gaia_dr3"
                for col in ["pmra_error", "pmdec_error", "parallax_over_error", 
                            "astrometric_excess_noise", "astrometric_excess_noise_sig", 
                            "ruwe", "r_cmodel_mag", "r_cmodel_magerr", 
                            "r_extendedness_value", "phot_g_mean_flux_over_error", 
                            "in_galaxy_candidates"]:
                    if col not in df_gaia.columns:
                        df_gaia[col] = np.nan
                
                # Save to a temporary CSV file
                import tempfile
                with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as tmp_f:
                    df_gaia.to_csv(tmp_f.name, index=False)
                    tmp_csv_path = tmp_f.name
                
                try:
                    gs_mag_min = p_cfg.get("netflow", {}).get("guidestars", {}).get("mag_min", 17.0)
                    gs_mag_max = p_cfg.get("netflow", {}).get("guidestars", {}).get("mag_max", 21.5)
                    guidestars = designutils.generate_guidestars_from_csv(
                        ra_tel,
                        dec_tel,
                        pa_tel,
                        obstime,
                        conf=config_toml if config_toml else {"sfa": {}},
                        guidestar_mag_min=gs_mag_min,
                        guidestar_mag_max=gs_mag_max,
                        guidestar_neighbor_mag_min=21.0,
                        guidestar_minsep_deg=1.0 / 3600,
                        gs_csv=tmp_csv_path
                    )
                    design.guideStars = guidestars
                    print(f"  Added {len(guidestars.objId)} guide stars from local catalog.")
                finally:
                    if os.path.exists(tmp_csv_path):
                        os.remove(tmp_csv_path)
            except Exception as e:
                print(f"  Warning: Failed to retrieve guide stars from local catalog: {e}")
                
        elif config_toml and "gaiadb" in config_toml:
            try:
                print("  Querying guide stars from GaiaDB...")
                gs_mag_min = p_cfg.get("netflow", {}).get("guidestars", {}).get("mag_min", config_toml.get("sfa", {}).get("guidestar_mag_min", 17.0))
                gs_mag_max = p_cfg.get("netflow", {}).get("guidestars", {}).get("mag_max", config_toml.get("sfa", {}).get("guidestar_mag_max", 21.5))
                guidestars = designutils.generate_guidestars_from_gaiadb(
                    ra_tel,
                    dec_tel,
                    pa_tel,
                    obstime,
                    conf=config_toml,
                    guidestar_mag_min=gs_mag_min,
                    guidestar_mag_max=gs_mag_max,
                    guidestar_neighbor_mag_min=config_toml.get("sfa", {}).get("guidestar_neighbor_mag_min", 21.0),
                    guidestar_minsep_deg=config_toml.get("sfa", {}).get("guidestar_minsep_deg", 1.0 / 3600),
                )
                design.guideStars = guidestars
                print(f"  Added {len(guidestars.objId)} guide stars.")
            except Exception as e:
                print(f"  Warning: Failed to retrieve guide stars from GaiaDB: {e}")
        else:
            print("  Skipping guide star retrieval (Gaiadb configuration / local catalog not provided).")
            
        # Save design fits
        out_filename = design.filename
        design.write(dirName=outdir, fileName=out_filename)
        print(f"  Saved pfsDesign fits: {os.path.join(outdir, out_filename)}")
        print(f"  pfsDesignId = 0x{design.pfsDesignId:016x}")

        # Collect info for OPE file
        summary_rows.append({
            "design_filename": out_filename,
            "observation_time": Time(obstime).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "ppc_code": ppc_code,
            "ppc_ra": ra_tel,
            "ppc_dec": dec_tel,
            "ppc_pa": pa_tel
        })
        ra_str = Angle(ra_tel * u.deg).to_string(unit=u.hourangle, sep="", precision=3, pad=True)
        dec_str = Angle(dec_tel * u.deg).to_string(unit=u.deg, sep="", alwayssign=True, precision=2, pad=True)
        obsdate_str = Time(obstime).strftime("%Y-%m-%d")
        ope_info_rows.append([
            ppc_code,         # val[0]: ppc_code
            obsdate_str,      # val[1]: obsdate_in_hst (使用 UTC date as approximation)
            obstime,          # val[2]: obstime_in_utc
            design.pfsDesignId,  # val[3]: pfs_design_id
            ra_str,           # val[4]: ppc_ra_str (HHMMSS)
            dec_str,          # val[5]: ppc_dec_str (DDMMSS)
            pa_tel,           # val[6]: ppc_pa
            obstime,          # val[7]: obstime_in_hst
            exptime_per_frame,  # val[8]: single_exptime
            n_frames,         # val[9]: n_split_frame
        ])

    # --- Generate OPE file ---
    if ope_template and ope_info_rows:
        try:
            from pfs_obsproc_planning.utils.make_opefile import OpeFile
            print(f"\n=== Generating OPE file ===")

            # Minimal conf dict required by OpeFile
            conf = {
                "ope": {
                    "template": ope_template,
                    "outfilePath": ope_outdir,
                    "designPath": os.path.abspath(outdir),
                    "backup_dates": [],
                },
                "ssp": {"ssp": False},
            }

            os.makedirs(ope_outdir, exist_ok=True)
            ope = OpeFile(conf=conf, workDir=os.getcwd())
            ope.loadTemplate(filename=ope_template)

            # Use the date from obstime (UTC)
            obsdate_utc = Time(obstime).strftime("%Y-%m-%d")
            ope.update_obsdate(obsdate_utc, utc=False)
            ope.update_design(ope_info_rows)
            ope.write()
            print(f"  OPE file saved: {ope.outfile}")
        except ImportError:
            print("  Warning: pfs_obsproc_planning not installed. Skipping OPE file generation.")
        except Exception as e:
            print(f"  Warning: OPE file generation failed: {e}")
    elif not ope_template:
        print("\n  (OPE file generation skipped: no --ope-template specified)")

    # --- Run validation ---
    if summary_rows:
        try:
            summary_df = pd.DataFrame(summary_rows)
            summary_csv_path = "summary_reconfigure_ppp-ppp+qplan_output.csv"
            summary_df.to_csv(summary_csv_path, index=False)
            print(f"  Saved summary CSV for validation: {summary_csv_path}")

            # Create a mock ppp/obList.ecsv to avoid KeyError: 'ob_code' in validation.py
            ppp_dir = "ppp"
            os.makedirs(ppp_dir, exist_ok=True)
            oblist_path = os.path.join(ppp_dir, "obList.ecsv")
            with open(oblist_path, "w") as f:
                f.write("# %ECSV 1.0\n# ---\n# datatype:\n# - {name: ob_code, datatype: string}\n# - {name: qa_reference_arm, datatype: string}\n# schema: astropy-2.0\nob_code qa_reference_arm\n")

            from pfs_obsproc_planning.utils.validation import validation as run_validation_tool
            
            print("\n=== Running Validation ===")
            # Prepare packages dirs
            import pfs.utils
            pfs_utils_dir = os.path.dirname(pfs.utils.__file__)
            
            import ets_fiber_assigner
            ets_fiber_assigner_dir = os.path.dirname(ets_fiber_assigner.__file__)
            
            # Prepare mock conf
            conf_validation = {
                "packages": {
                    "pfs_instdata_dir": pfs_instdata_dir,
                    "pfs_utils_dir": pfs_utils_dir,
                },
                "sfa": {
                    "cobra_coach_dir": ets_fiber_assigner_dir,
                    "dot_margin": p_cfg.get("pfs", {}).get("black_dot_radius_margin", 1.65) if p_cfg else 1.65,
                    "fill_unassign_radius_check": 1.5 / 60.0,
                },
                "validation": {
                    "save_unassign_toobright": False,
                },
                "ppp": {
                    "mode": "queue",
                    "proposalIds": [config_yaml.get("proposal_id", "S26A-OT02")],
                }
            }
            
            # Load local Gaia catalog for mock DB queries during validation
            df_gaia_all = pd.DataFrame()
            if gaia_catalog and os.path.exists(gaia_catalog):
                try:
                    print(f"  Mocking Gaia DB queries using local catalog: {gaia_catalog}")
                    t_gaia = Table.read(gaia_catalog, format="ascii.ecsv")
                    df_gaia_all = t_gaia.to_pandas()
                except Exception as e:
                    print(f"  Warning: Failed to load gaia catalog for mock: {e}")

            class DummyGuideStars:
                def __init__(self):
                    self.agId = []
                    self.objId = []
                    self.ra = []
                    self.dec = []
                    self.magnitude = []
                    self.passband = []

            def mock_generate_targets_from_gaiadb(ra, dec, conf=None, fp_radius_degree=260.0 * 10.2 / 3600, fp_fudge_factor=1.5, search_radius=None, band_select="phot_g_mean_mag", mag_min=0.0, mag_max=99.0, good_astrometry=False, write_csv=False):
                if df_gaia_all.empty:
                    return pd.DataFrame()
                if search_radius is None:
                    search_radius = fp_radius_degree * fp_fudge_factor
                cos_dec = np.cos(np.deg2rad(dec))
                dist = np.hypot((df_gaia_all["ra"] - ra) * cos_dec, df_gaia_all["dec"] - dec)
                mask = (dist <= search_radius)
                mag_col = band_select if band_select in df_gaia_all.columns else ("phot_g_mean_mag" if "phot_g_mean_mag" in df_gaia_all.columns else "magnitude")
                if mag_col in df_gaia_all.columns:
                    mask &= (df_gaia_all[mag_col] >= mag_min) & (df_gaia_all[mag_col] <= mag_max)
                return df_gaia_all[mask].copy()

            def mock_generate_guidestars_from_gaiadb(ra, dec, pa, observation_time, telescope_elevation=None, conf=None, guidestar_mag_min=12.0, guidestar_mag_max=19.0, guidestar_neighbor_mag_min=21.0, guidestar_minsep_deg=1.0/3600, fp_radius_degree=260.0*10.2/3600, fp_fudge_factor=1.5, search_radius=None, gaiadb_input_catalog_id=4, guide_star_id_exclude=[], good_astrometry=False, gs_snr_thresh=5.0):
                if df_gaia_all.empty:
                    return DummyGuideStars()
                
                import tempfile
                df_gaia_tmp = df_gaia_all.copy()
                df_gaia_tmp["magnitude"] = df_gaia_tmp["phot_g_mean_mag"] if "phot_g_mean_mag" in df_gaia_tmp.columns else df_gaia_tmp["magnitude"]
                df_gaia_tmp["color"] = df_gaia_tmp["bp_rp"] if "bp_rp" in df_gaia_tmp.columns else df_gaia_tmp["color"]
                
                if "catalog" not in df_gaia_tmp.columns:
                    df_gaia_tmp["catalog"] = "gaia_dr3"
                for col in ["pmra_error", "pmdec_error", "parallax_over_error", 
                            "astrometric_excess_noise", "astrometric_excess_noise_sig", 
                            "ruwe", "r_cmodel_mag", "r_cmodel_magerr", 
                            "r_extendedness_value", "phot_g_mean_flux_over_error", 
                            "in_galaxy_candidates"]:
                    if col not in df_gaia_tmp.columns:
                        df_gaia_tmp[col] = np.nan
                        
                with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as tmp_f:
                    df_gaia_tmp.to_csv(tmp_f.name, index=False)
                    tmp_csv_path = tmp_f.name
                
                try:
                    res = designutils.generate_guidestars_from_csv(
                        ra, dec, pa, observation_time,
                        conf=conf if conf else {"sfa": {}},
                        guidestar_mag_min=guidestar_mag_min,
                        guidestar_mag_max=guidestar_mag_max,
                        guidestar_neighbor_mag_min=guidestar_neighbor_mag_min,
                        guidestar_minsep_deg=guidestar_minsep_deg,
                        gs_csv=tmp_csv_path
                    )
                    return res
                except Exception as e:
                    print(f"  Warning: Mock guidestars generation failed: {e}")
                    return DummyGuideStars()
                finally:
                    if os.path.exists(tmp_csv_path):
                        os.remove(tmp_csv_path)

            orig_generate_targets = dbutils.generate_targets_from_gaiadb
            orig_generate_guidestars = designutils.generate_guidestars_from_gaiadb

            dbutils.generate_targets_from_gaiadb = mock_generate_targets_from_gaiadb
            designutils.generate_guidestars_from_gaiadb = mock_generate_guidestars_from_gaiadb

            
            try:
                run_validation_tool(
                    parentPath=".",
                    figpath=outdir,
                    save=True,
                    show=False,
                    ssp=False,
                    conf=conf_validation
                )
                print(f"  Validation completed successfully. HTML report saved in {outdir}")
            finally:
                # Restore original database methods
                dbutils.generate_targets_from_gaiadb = orig_generate_targets
                designutils.generate_guidestars_from_gaiadb = orig_generate_guidestars
                
                # Clean up temporary validation files
                if os.path.exists(summary_csv_path):
                    os.remove(summary_csv_path)
                if os.path.exists(oblist_path):
                    os.remove(oblist_path)
                if os.path.exists(ppp_dir) and not os.listdir(ppp_dir):
                    os.rmdir(ppp_dir)
                
        except ImportError:
            print("\n  Warning: pfs_obsproc_planning not installed. Skipping validation.")
        except Exception as e:
            print(f"\n  Warning: Validation failed: {e}")

def main():
    parser = argparse.ArgumentParser(description="Create pfsDesign files from run_netflow.py fiber assignment results.")
    parser.add_argument("--pointing", default="optimized_pointings.ecsv", help="Pointing list ecsv file (default: optimized_pointings.ecsv)")
    parser.add_argument("--targets-dir", default="targets", help="Directory where target ECSV files are saved (default: targets)")
    parser.add_argument("--outdir", default="design", help="Output directory for generated pfsDesign fits files (default: design)")
    parser.add_argument("--config-toml", default=None, help="Path to config.toml for Gaiadb configuration")
    parser.add_argument("--gaia-catalog", default="cosmos/gaia.ecsv", help="Path to local Gaia ECSV catalog for guide stars (default: cosmos/gaia.ecsv)")
    parser.add_argument("--obstime", default="2026-05-09T06:00:00Z", help="Observing time in UTC (default: 2026-05-09T06:00:00Z)")
    parser.add_argument("--ope-template", default=None, help="Path to OPE file template. If not specified, OPE file generation is skipped.")
    parser.add_argument("--ope-outdir", default="ope", help="Output directory for generated OPE files (default: ope)")
    parser.add_argument("--exptime-per-frame", type=float, default=900.0, help="Exposure time per sub-frame in seconds (default: 900.0)")
    parser.add_argument("--n-frames", type=int, default=4, help="Number of sub-frames per pointing (default: 4)")
    parser.add_argument("--pipeline-config", default=None, help="Path to pipeline configuration YAML file. If specified, overrides arguments with values from the configuration file.")
    
    args = parser.parse_args()
    
    generate_pfs_designs(
        pointing_file=args.pointing,
        targets_dir=args.targets_dir,
        outdir=args.outdir,
        config_toml_path=args.config_toml,
        gaia_catalog=args.gaia_catalog,
        obstime=args.obstime,
        ope_template=args.ope_template,
        ope_outdir=args.ope_outdir,
        exptime_per_frame=args.exptime_per_frame,
        n_frames=args.n_frames,
        pipeline_config_path=args.pipeline_config,
        pipeline_config=None,
    )


if __name__ == "__main__":
    main()

