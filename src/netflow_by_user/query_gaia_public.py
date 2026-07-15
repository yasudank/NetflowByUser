#!/usr/bin/env python3
import os
import sys
import toml
import pandas as pd
import numpy as np
from astropy.table import Table
from tqdm import tqdm
from astroquery.gaia import Gaia

# Suppress astroquery / astropy warnings if needed
import warnings
warnings.filterwarnings("ignore")

def get_search_radius(
    fp_radius_degree = 260.0 * 10.2 / 3600,  # "Radius" of PFS FoV in degree
    fp_fudge_factor = 1.5  # fudge factor for search widths
):
    search_radius = fp_radius_degree * fp_fudge_factor
    print("search_radius is %f degree." % search_radius)
    return search_radius

def get_config(config_fn):
    with open(config_fn, "r") as f:
        config = toml.load(f)
    print(config)
    return config

def get_centerList(config):
    fn = os.path.join(config["input"]["dir"], config["input"]["fn_ppcList"]) 
    ppcList = Table.read(fn)
    print("There are %d pointings." % len(ppcList))
    print("ppcList read from %s" % fn)
    return ppcList

def get_guidestar_candidates(config, ra, dec, search_radius, mag_min=12.0, mag_max=21.5):
    """Query guide star candidates from public ESA Gaia DR3 Archive"""
    
    # Required columns in the output format
    COLUMNS = [
        "source_id", "ra", "dec", "parallax", "pmra", "pmdec", "ref_epoch", "phot_g_mean_mag", "bp_rp",
        "pmra_error", "pmdec_error", "parallax_error", "astrometric_excess_noise",
        "astrometric_excess_noise_sig", "ruwe", "phot_g_mean_flux_over_error"
    ]

    query_string = f"""
    SELECT source_id, ra, dec, parallax, pmra, pmdec, ref_epoch, phot_g_mean_mag, bp_rp,
           pmra_error, pmdec_error, parallax_error, astrometric_excess_noise,
           astrometric_excess_noise_sig, ruwe, phot_g_mean_flux_over_error
    FROM gaiadr3.gaia_source
    WHERE 1=CONTAINS(
      POINT('ICRS', ra, dec),
      CIRCLE('ICRS', {ra}, {dec}, {search_radius})
    )
    AND phot_g_mean_mag BETWEEN {mag_min} AND {mag_max}
    """
    
    try:
        job = Gaia.launch_job_async(query_string)
        r = job.get_results()
        df = r.to_pandas()
    except Exception as e:
        print(f"Error querying Gaia Archive: {e}")
        # Return an empty dataframe with correct columns
        df = pd.DataFrame(columns=COLUMNS)

    # Ensure all required columns are present and in the exact correct order
    df = df.reindex(columns=COLUMNS)
    
    # Cast source_id to nullable integer type Int64 so it writes cleanly to ECSV
    if "source_id" in df.columns:
        df["source_id"] = df["source_id"].astype("Int64")
        
    return df

def main():
    if len(sys.argv) > 1:
        config_file = sys.argv[1]
    else:
        config_file = 'config_targetdb_cosmos.toml'

    # set search radius in degree
    search_radius = get_search_radius()

    # read config from config_file
    config = get_config(config_file)

    # read center list
    ppcList = get_centerList(config)
    ppc_code_list, ra_list, dec_list = ppcList['ppc_code'], ppcList['ppc_ra'], ppcList['ppc_dec']

    gaia_dir = os.path.join(config['output']['dir'], "gaia")
    if not os.path.exists(gaia_dir):
        os.makedirs(gaia_dir, exist_ok=True)
    
    mag_min = config['gaiadb']['mag_min']
    mag_max = config['gaiadb']['mag_max']

    for ppc_code, ra, dec in tqdm(
        zip(ppc_code_list, ra_list, dec_list),
        total=len(ppc_code_list),
        desc="Querying Gaia pointings"
    ):
        df = get_guidestar_candidates(config, ra, dec, search_radius, mag_min, mag_max)
        outfn = os.path.join(gaia_dir, f"{ppc_code}.ecsv")

        table = Table.from_pandas(df)
        table.write(outfn, format="ascii.ecsv", overwrite=True)

        tqdm.write("%s: %d gaia stars selected." % (ppc_code, len(df)))
        tqdm.write('write to %s'%outfn)

if __name__ == "__main__":
    main()
