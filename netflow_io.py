import os
import numpy as np
import pandas as pd
from astropy.table import Table
import ets_fiber_assigner.netflow as nf

def readScientificFromCSV(file, prefix):
    """Read a set of scientific targets from a CSV file"""
    FEATURE_FLAG_IDX_N2 = 0  # we use bit #0 for this feature

    t = Table.read(file, format="ascii.csv")
    res = []
    for r in t:
        try:
            stg = int(r["stage"])
        except KeyError:
            stg = 0

        req_flags = 0  # require all features by default
        try:
            req_flags = (1 << FEATURE_FLAG_IDX_N2) if r["reference_arm"] != "n" else 0
        except KeyError:
            try:
                req_flags = (1 << FEATURE_FLAG_IDX_N2) if r["qa_reference_arm"] != "n" else 0
            except KeyError:
                req_flags = 0

        id_val = r["ob_code"] if "ob_code" in t.colnames else (r["obj_id"] if "obj_id" in t.colnames else r["ID"])
        ra_val = r["ra"] if "ra" in t.colnames else r["R.A."]
        dec_val = r["dec"] if "dec" in t.colnames else r["Dec."]
        exptime_val = r["exptime"] if "exptime" in t.colnames else r["Exposure Time"]
        priority_val = r["priority"] if "priority" in t.colnames else r["Priority"]

        st = nf.ScienceTarget(id_val, ra_val, dec_val,
                              exptime_val, priority_val, prefix, stage=stg,
                              req_flags=req_flags)

        for col in t.colnames:
            try:
                setattr(st, col, r[col])
            except AttributeError:
                setattr(st, "_" + col, r[col])
        res.append(st)
    return res


def readCalibrationFromECSV(file, targetclass):
    """Read a set of calibration targets from an ECSV file"""
    t = Table.read(file, format="ascii.ecsv")
    res = []
    for r in t:
        try:
            stg = int(r["stage"])
        except KeyError:
            stg = 0

        if "fluxstd_id" in t.colnames:
            id_val = '_fluxstd_' + str(r["fluxstd_id"])
        elif "sky_id" in t.colnames:
            id_val = '_sky_' + str(r["sky_id"])
        elif "obj_id" in t.colnames:
            id_val = r["obj_id"]
        elif "ID" in t.colnames:
            id_val = r["ID"]
        else:
            raise KeyError("No target ID column found in ECSV file.")

        ra_val = r["ra"] if "ra" in t.colnames else r["R.A."]
        dec_val = r["dec"] if "dec" in t.colnames else r["Dec."]

        pmra_val = r["pmra"] if "pmra" in t.colnames else 0.0
        pmdec_val = r["pmdec"] if "pmdec" in t.colnames else 0.0
        parallax_val = r["parallax"] if "parallax" in t.colnames else 0.0

        epoch_val = r["epoch"] if "epoch" in t.colnames else 2000.0
        if isinstance(epoch_val, str) and epoch_val.startswith("J"):
            try:
                epoch_val = float(epoch_val[1:])
            except ValueError:
                epoch_val = 2000.0

        penalty_val = r["penalty"] if "penalty" in t.colnames else 0.0

        st = nf.CalibTarget(id_val, ra_val, dec_val, targetclass,
                            penalty=penalty_val, pmra=pmra_val, pmdec=pmdec_val,
                            parallax=parallax_val, epoch=epoch_val, stage=stg)

        for col in t.colnames:
            val = r[col]
            if col == "epoch" and isinstance(val, str) and val.startswith("J"):
                try:
                    val = float(val[1:])
                except ValueError:
                    pass
            try:
                setattr(st, col, val)
            except AttributeError:
                setattr(st, "_" + col, val)
        res.append(st)
    return res


def _fluxstd_cuts(df_fluxstd, mag_min=17.0, mag_max=19.0, prob_threshold=0.5):
    mag_mask = (df_fluxstd["psf_mag_g"] > mag_min) & (df_fluxstd["psf_mag_g"] < mag_max)
    df_fluxstd = df_fluxstd.loc[mag_mask].copy()

    if df_fluxstd.empty:
        return df_fluxstd

    prob_mask = (df_fluxstd["prob_f_star"] > prob_threshold) | (df_fluxstd['is_fstar_gaia'] == True)
    df_selected = df_fluxstd.loc[prob_mask].copy()

    df_selected = df_selected.drop_duplicates(subset=['fluxstd_id'])
    return df_selected


def readFluxstdFromECSV(file, targetclass, mag_min=17.0, mag_max=19.0):
    t = Table.read(file, format="ascii.ecsv")
    df = t.to_pandas()
    df_selected = _fluxstd_cuts(df, mag_min=mag_min, mag_max=mag_max)

    res = []
    for idx, r in df_selected.iterrows():
        try:
            stg = int(r["stage"])
        except KeyError:
            stg = 0

        id_val = r["fluxstd_id"] if "fluxstd_id" in df_selected.columns else (r["obj_id"] if "obj_id" in df_selected.columns else r["ID"])
        ra_val = r["ra"] if "ra" in df_selected.columns else r["R.A."]
        dec_val = r["dec"] if "dec" in df_selected.columns else r["Dec."]

        pmra_val = r["pmra"] if "pmra" in df_selected.columns else 0.0
        pmdec_val = r["pmdec"] if "pmdec" in df_selected.columns else 0.0
        parallax_val = r["parallax"] if "parallax" in df_selected.columns else 0.0

        epoch_val = r["epoch"] if "epoch" in df_selected.columns else 2000.0
        if isinstance(epoch_val, str) and epoch_val.startswith("J"):
            try:
                epoch_val = float(epoch_val[1:])
            except ValueError:
                epoch_val = 2000.0

        penalty_val = r["penalty"] if "penalty" in df_selected.columns else 0.0

        st = nf.CalibTarget(id_val, ra_val, dec_val, targetclass,
                            penalty=penalty_val, pmra=pmra_val, pmdec=pmdec_val,
                            parallax=parallax_val, epoch=epoch_val, stage=stg)

        for col in df_selected.columns:
            val = r[col]
            if isinstance(val, float) and np.isnan(val):
                val = None
            if col == "epoch" and isinstance(val, str) and val.startswith("J"):
                try:
                    val = float(val[1:])
                except ValueError:
                    pass
            try:
                setattr(st, col, val)
            except AttributeError:
                setattr(st, "_" + col, val)
        res.append(st)
    return res


def load_all_targets(fscience_targets, fcal_stars, fsky_pos, fluxstd_mag_min=17.0, fluxstd_mag_max=19.0):
    """Load and merge science, fluxstd, and sky targets."""
    print("Reading science targets...")
    tgt = readScientificFromCSV(fscience_targets, "sci")
    print(f"Done. {len(tgt)} targets are read")

    print("Reading fluxstd targets...")
    tgt += readFluxstdFromECSV(fcal_stars, "cal", mag_min=fluxstd_mag_min, mag_max=fluxstd_mag_max)
    print(f"Done. {len(tgt)} targets are read")

    print("Reading sky targets...")
    tgt += readCalibrationFromECSV(fsky_pos, "sky")
    print(f"Done. {len(tgt)} targets are read")

    return tgt


def save_targets_ecsv(exposures_data, fscience_targets, fcal_stars, fsky_pos, targets_dir):
    """Match solved assignments back to raw target catalogs and save partitioned ECSV files."""
    if not exposures_data:
        return

    print("Matching targets and saving ECSV files...")
    df_out = pd.DataFrame(exposures_data)

    df_sci = pd.read_csv(fscience_targets)
    if 'ob_code' in df_sci.columns:
        df_sci['target_id'] = df_sci['ob_code'].astype(str)
    elif 'obj_id' in df_sci.columns:
        df_sci['target_id'] = df_sci['obj_id'].astype(str)
    else:
        df_sci['target_id'] = df_sci['ID'].astype(str)

    t_cal = Table.read(fcal_stars, format="ascii.ecsv").to_pandas()
    if 'fluxstd_id' in t_cal.columns:
        t_cal['target_id'] = t_cal['fluxstd_id'].astype(str)
    else:
        t_cal['target_id'] = t_cal['obj_id'].astype(str)

    t_sky = Table.read(fsky_pos, format="ascii.ecsv").to_pandas()
    if 'sky_id' in t_sky.columns:
        t_sky['target_id'] = '_sky_' + t_sky['sky_id'].astype(str)
    else:
        t_sky['target_id'] = t_sky['obj_id'].astype(str)

    # Create target directories
    for ttype in ['science', 'fluxstd', 'sky']:
        dirpath = os.path.join(targets_dir, ttype)
        if not os.path.exists(dirpath):
            os.makedirs(dirpath)

    for ppc_code, group in df_out.groupby('ppc_code'):
        group_merge = group.drop(columns=['ppc_code'])

        # Science
        sci_merge = pd.merge(df_sci, group_merge, on='target_id', how='inner')
        if not sci_merge.empty:
            sci_merge = sci_merge.drop(columns=['target_id'])
            cols = [c for c in sci_merge.columns if c not in ['cobraId', 'pfi_X', 'pfi_Y']] + ['cobraId', 'pfi_X', 'pfi_Y']
            sci_table = Table.from_pandas(sci_merge[cols])
            outpath = os.path.join(targets_dir, "science", f"{ppc_code}.ecsv")
            sci_table.write(outpath, format="ascii.ecsv", overwrite=True)
            print(f"Saved {outpath} ({len(sci_table)} rows)")

        # Fluxstd
        cal_merge = pd.merge(t_cal, group_merge, on='target_id', how='inner')
        if not cal_merge.empty:
            cal_merge['ob_code'] = cal_merge['target_id']
            cal_merge = cal_merge.drop(columns=['target_id'])
            cols = [c for c in cal_merge.columns if c not in ['cobraId', 'pfi_X', 'pfi_Y']] + ['cobraId', 'pfi_X', 'pfi_Y']
            cal_table = Table.from_pandas(cal_merge[cols])
            outpath = os.path.join(targets_dir, "fluxstd", f"{ppc_code}.ecsv")
            cal_table.write(outpath, format="ascii.ecsv", overwrite=True)
            print(f"Saved {outpath} ({len(cal_table)} rows)")

        # Sky
        sky_merge = pd.merge(t_sky, group_merge, on='target_id', how='inner')
        if not sky_merge.empty:
            sky_merge['ob_code'] = sky_merge['target_id']
            sky_merge = sky_merge.drop(columns=['target_id'])
            cols = [c for c in sky_merge.columns if c not in ['cobraId', 'pfi_X', 'pfi_Y']] + ['cobraId', 'pfi_X', 'pfi_Y']
            sky_table = Table.from_pandas(sky_merge[cols])
            outpath = os.path.join(targets_dir, "sky", f"{ppc_code}.ecsv")
            sky_table.write(outpath, format="ascii.ecsv", overwrite=True)
            print(f"Saved {outpath} ({len(sky_table)} rows)")
