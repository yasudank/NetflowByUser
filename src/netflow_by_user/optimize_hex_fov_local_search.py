#!/usr/bin/env python3
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import argparse
import numpy as np
import pandas as pd
from astropy.table import Table
from tqdm import tqdm

from optimize_hex_fov_with_guidestars import evaluate_guidestars_single, check_bright_stars_near_broken_fibers

def is_valid_pointing(cam_star_counts, min_stars_per_cam=2, min_cams_with_stars=6):
    """Check if the pointing satisfies the guide star constraints."""
    # Check for saturated cameras (-999)
    if any(count < 0 for count in cam_star_counts):
        return False
    # Check minimum number of cameras with enough stars
    cams_ok = sum(1 for count in cam_star_counts if count >= min_stars_per_cam)
    return cams_ok >= min_cams_with_stars

def local_search(ra_center, dec_center, pa_center, df_gaia, obstime,
                 min_stars_per_cam=2, min_cams_with_stars=6,
                 ra_range=0.1, dec_range=0.1, pa_range=5.0,
                 pos_step=0.02, pa_step=1.0,
                 avoid_gaps=False, all_pointings=None, current_idx=-1,
                 bench=None, bright_star_mag_limit=12.0, bright_star_radius_arcmin=1.5):
    
    cos_dec = np.cos(np.radians(dec_center))
    ra_step_adj = pos_step / cos_dec if cos_dec > 0.1 else pos_step
    ra_range_adj = ra_range / cos_dec if cos_dec > 0.1 else ra_range
    
    ra_offsets = np.arange(-ra_range_adj, ra_range_adj + 1e-5, ra_step_adj)
    dec_offsets = np.arange(-dec_range, dec_range + 1e-5, pos_step)
    pa_offsets = np.arange(-pa_range, pa_range + 1e-5, pa_step)
    
    neighbors = []
    if avoid_gaps and all_pointings is not None:
        for j, row in enumerate(all_pointings):
            if j == current_idx: continue
            n_ra = row["ppc_ra"]
            n_dec = row["ppc_dec"]
            dist = np.sqrt(((n_ra - ra_center) * cos_dec)**2 + (n_dec - dec_center)**2)
            if dist < 1.5:
                neighbors.append({'ra': n_ra, 'dec': n_dec, 'initial_dist': dist})
    
    # Create a list of all offsets and sort them by distance from the center
    candidates = []
    for d_ra in ra_offsets:
        for d_dec in dec_offsets:
            for d_pa in pa_offsets:
                dist_deg = np.sqrt((d_ra * cos_dec)**2 + d_dec**2)
                
                gap_penalty = 0.0
                if avoid_gaps and neighbors:
                    cand_ra = ra_center + d_ra
                    cand_dec = dec_center + d_dec
                    for n in neighbors:
                        cand_dist = np.sqrt(((n['ra'] - cand_ra) * cos_dec)**2 + (n['dec'] - cand_dec)**2)
                        # Penalty if distance to neighbor increases beyond its initial distance
                        gap_penalty += max(0, cand_dist - n['initial_dist'])
                        
                dist_metric = dist_deg + abs(d_pa) * 0.0001 + 10.0 * gap_penalty
                candidates.append((dist_metric, d_ra, d_dec, d_pa))
                
    candidates.sort(key=lambda x: x[0])
    
    for _, d_ra, d_dec, d_pa in candidates:
        cand_ra = ra_center + d_ra
        cand_dec = dec_center + d_dec
        cand_pa = pa_center + d_pa
        
        counts, _ = evaluate_guidestars_single(
            cand_ra, cand_dec, cand_pa, df_gaia, obstime,
            min_mag=12.0, max_mag=21.5, minsep_arcsec=1.0
        )
        
        if is_valid_pointing(counts, min_stars_per_cam, min_cams_with_stars):
            # Check for bright star near broken fibers
            if not check_bright_stars_near_broken_fibers(
                cand_ra, cand_dec, cand_pa, df_gaia, obstime, bench,
                radius_deg=bright_star_radius_arcmin/60.0,
                max_mag=bright_star_mag_limit
            ):
                return (cand_ra, cand_dec, cand_pa, counts)
            
    return None

def main():
    # 0. Pre-parse config option
    conf_parser = argparse.ArgumentParser(add_help=False)
    conf_parser.add_argument("--config", "-c", help="Path to YAML configuration file")
    args_conf, remaining_argv = conf_parser.parse_known_args()

    defaults = {
        "input": None,
        "output": None,
        "gaia": "cosmos/gaia.ecsv",
        "obstime": "2026-05-09T06:00:00Z",
        "min_stars": 2,
        "min_cams": 6,
        "search_radius": 0.1,
        "search_step": 0.02,
        "pa_radius": 5.0,
        "pa_step": 1.0,
        "avoid_gaps": False,
        "bright_star_mag_limit": 12.0,
        "bright_star_radius_arcmin": 1.5,
    }

    if args_conf.config and os.path.exists(args_conf.config):
        import yaml
        print(f"Loading configurations from YAML: {args_conf.config}")
        with open(args_conf.config, "r") as f:
            cfg = yaml.safe_load(f)
        if cfg:
            pointing_file = cfg.get("inputs", {}).get("pointing_file", None)
            if pointing_file:
                defaults["input"] = pointing_file
            
            gaia_cat = cfg.get("inputs", {}).get("gaia_catalog", "")
            if gaia_cat:
                defaults["gaia"] = gaia_cat
                
            obstime = cfg.get("obstime", "")
            if obstime:
                defaults["obstime"] = obstime
                
            netflow_cfg = cfg.get("netflow", {})
            if "min_stars_per_cam" in netflow_cfg:
                defaults["min_stars"] = netflow_cfg["min_stars_per_cam"]
            if "min_cams_with_stars" in netflow_cfg:
                defaults["min_cams"] = netflow_cfg["min_cams_with_stars"]
            if "bright_star_mag_limit" in netflow_cfg:
                defaults["bright_star_mag_limit"] = netflow_cfg["bright_star_mag_limit"]
            if "bright_star_radius_arcmin" in netflow_cfg:
                defaults["bright_star_radius_arcmin"] = netflow_cfg["bright_star_radius_arcmin"]

    parser = argparse.ArgumentParser(description="Local search to satisfy guide star constraints for a pointing list.")
    parser.add_argument("--config", "-c", help="Path to YAML configuration file")
    parser.add_argument("--input", "-i", required=(defaults["input"] is None), default=defaults["input"], help="Input ECSV file with ppc_ra, ppc_dec, ppc_pa")
    parser.add_argument("--output", "-o", required=(defaults["output"] is None), default=defaults["output"], help="Output ECSV file")
    parser.add_argument("--gaia", "-g", default=defaults["gaia"], help="Gaia catalog ECSV file")
    parser.add_argument("--obstime", "-t", default=defaults["obstime"], help="Observation time")
    parser.add_argument("--min_stars", type=int, default=defaults["min_stars"], help="Minimum guide stars per camera")
    parser.add_argument("--min_cams", type=int, default=defaults["min_cams"], help="Minimum guide cameras with stars")
    parser.add_argument("--bright-star-mag-limit", type=float, default=defaults["bright_star_mag_limit"], help="Bright star magnitude limit to avoid around broken fibers")
    parser.add_argument("--bright-star-radius-arcmin", type=float, default=defaults["bright_star_radius_arcmin"], help="Radius in arcminutes to avoid bright stars around broken fibers")
    
    # Search parameters
    parser.add_argument("--search_radius", type=float, default=defaults["search_radius"], help="Spatial search radius (deg)")
    parser.add_argument("--search_step", type=float, default=defaults["search_step"], help="Spatial search step (deg)")
    parser.add_argument("--pa_radius", type=float, default=defaults["pa_radius"], help="PA search radius (deg)")
    parser.add_argument("--pa_step", type=float, default=defaults["pa_step"], help="PA search step (deg)")
    parser.add_argument("--avoid-gaps", action="store_true", default=defaults["avoid_gaps"], help="Avoid creating gaps between adjacent pointings")
    parser.add_argument("--add-columns", action="store_true", help="Add guide star count columns (ag0..ag5, n_guidestars) to the output ECSV")
    parser.add_argument("--report-file", default=None, help="Save guide star report to a text file")
    
    args = parser.parse_args(remaining_argv)
    
    print(f"Reading pointings from {args.input}...")
    t_in = Table.read(args.input, format="ascii.ecsv")
    
    if "ppc_pa" not in t_in.colnames:
        t_in["ppc_pa"] = 0.0
        
    print(f"Reading Gaia catalog from {args.gaia}...")
    t_gaia = Table.read(args.gaia, format="ascii.ecsv")
    df_gaia = t_gaia.to_pandas()
    df_gaia["magnitude"] = df_gaia["phot_g_mean_mag"]
    df_gaia["color"] = df_gaia["bp_rp"]
    df_gaia = df_gaia.fillna({"parallax": 1.0e-07, "pmra": 0.0, "pmdec": 0.0})
    
    for col in ["pmra_error", "pmdec_error", "parallax_over_error"]:
        if col not in df_gaia.columns:
            df_gaia[col] = np.nan
    # Load instrument model (bench) to identify broken fibers
    try:
        import netflow_instrument
        bench = netflow_instrument.getBench()
        broken_centers = bench.cobras.centers[~bench.cobras.isGood]
        print(f"Local search: using bench with {len(broken_centers)} broken cobras.")
    except Exception as e:
        print(f"Warning: could not initialize bench dynamically in local search: {e}")
        bench = None
            
    success_count = 0
    fail_count = 0
    adjusted_count = 0
    report_rows = []
    
    print("Evaluating pointings...")
    for i, row in enumerate(tqdm(t_in)):
        ppc_code = str(row["ppc_code"]) if "ppc_code" in t_in.colnames else f"POINTING_{i+1}"
        ra = float(row["ppc_ra"])
        dec = float(row["ppc_dec"])
        pa = float(row["ppc_pa"])
        
        counts, _ = evaluate_guidestars_single(
            ra, dec, pa, df_gaia, args.obstime,
            min_mag=12.0, max_mag=21.5, minsep_arcsec=1.0
        )
        
        is_valid = is_valid_pointing(counts, args.min_stars, args.min_cams)
        if is_valid:
            # Check for bright star near broken fibers
            is_valid = not check_bright_stars_near_broken_fibers(
                ra, dec, pa, df_gaia, args.obstime, bench,
                radius_deg=args.bright_star_radius_arcmin/60.0,
                max_mag=args.bright_star_mag_limit
            )
            
        if is_valid:
            success_count += 1
            status = "Valid"
            final_counts = counts
            cams_ok = sum(1 for c in final_counts if c >= args.min_stars)
            total_stars = sum(c for c in final_counts if c > 0)
            tqdm.write(
                f"[{i+1}/{len(t_in)}] {ppc_code}: Valid (No adjustment needed) -> "
                f"RA={ra:.5f}, Dec={dec:.5f}, PA={pa:5.1f}° | "
                f"AG0..AG5: {final_counts} (Total: {total_stars}, Cams OK: {cams_ok}/{len(final_counts)})"
            )
            report_rows.append({
                "idx": i + 1,
                "ppc_code": ppc_code,
                "status": status,
                "ra": ra, "dec": dec, "pa": pa,
                "counts": final_counts,
                "total": total_stars,
                "cams_ok": cams_ok,
            })
        else:
            best_cand = local_search(
                ra, dec, pa, df_gaia, args.obstime,
                min_stars_per_cam=args.min_stars,
                min_cams_with_stars=args.min_cams,
                ra_range=args.search_radius,
                dec_range=args.search_radius,
                pa_range=args.pa_radius,
                pos_step=args.search_step,
                pa_step=args.pa_step,
                avoid_gaps=args.avoid_gaps,
                all_pointings=t_in,
                current_idx=i,
                bench=bench,
                bright_star_mag_limit=args.bright_star_mag_limit,
                bright_star_radius_arcmin=args.bright_star_radius_arcmin
            )
            
            if best_cand is not None:
                new_ra, new_dec, new_pa, final_counts = best_cand
                row["ppc_ra"] = new_ra
                row["ppc_dec"] = new_dec
                row["ppc_pa"] = new_pa
                success_count += 1
                adjusted_count += 1
                status = "Adjusted"
                cams_ok = sum(1 for c in final_counts if c >= args.min_stars)
                total_stars = sum(c for c in final_counts if c > 0)
                tqdm.write(
                    f"[{i+1}/{len(t_in)}] {ppc_code}: Adjusted -> "
                    f"RA={new_ra:.5f}, Dec={new_dec:.5f}, PA={new_pa:5.1f}° | "
                    f"AG0..AG5: {final_counts} (Total: {total_stars}, Cams OK: {cams_ok}/{len(final_counts)})"
                )
                report_rows.append({
                    "idx": i + 1,
                    "ppc_code": ppc_code,
                    "status": status,
                    "ra": new_ra, "dec": new_dec, "pa": new_pa,
                    "counts": final_counts,
                    "total": total_stars,
                    "cams_ok": cams_ok,
                })
            else:
                fail_count += 1
                status = "Failed"
                final_counts = counts
                cams_ok = sum(1 for c in final_counts if c >= args.min_stars)
                total_stars = sum(c for c in final_counts if c > 0)
                tqdm.write(
                    f"[{i+1}/{len(t_in)}] {ppc_code}: FAILED (No valid pointing found within search radius) | "
                    f"Initial AG0..AG5: {final_counts} (Cams OK: {cams_ok}/{len(final_counts)})"
                )
                report_rows.append({
                    "idx": i + 1,
                    "ppc_code": ppc_code,
                    "status": status,
                    "ra": ra, "dec": dec, "pa": pa,
                    "counts": final_counts,
                    "total": total_stars,
                    "cams_ok": cams_ok,
                })
                
    # Build guide star report table
    code_width = max(15, max((len(str(r["ppc_code"])) for r in report_rows), default=15))
    cams_ok_hdr = f"Cams OK(>={args.min_stars})"
    
    header = (
        f"{'#':<3}  {'Pointing Code':<{code_width}}  {'Status':<8}  "
        f"{'RA [deg]':>10}  {'Dec [deg]':>10}  {'PA [deg]':>8}  "
        f"{'AG0':>5}  {'AG1':>5}  {'AG2':>5}  {'AG3':>5}  {'AG4':>5}  {'AG5':>5}  "
        f"{'Total':>6}  {cams_ok_hdr:>14}"
    )
    separator = "-" * len(header)
    double_separator = "=" * len(header)
    
    table_lines = [
        "",
        double_separator,
        "Guide Stars Report per Pointing:",
        separator,
        header,
        separator,
    ]
    
    for r in report_rows:
        def fmt_cam(c):
            if c < 0:
                return "  SAT"
            return f"{c:5d}"
            
        c_strs = "  ".join(fmt_cam(c) for c in r["counts"])
        cams_ok_str = f"{r['cams_ok']}/{len(r['counts'])}"
        line = (
            f"{r['idx']:<3d}  {r['ppc_code']:<{code_width}}  {r['status']:<8}  "
            f"{r['ra']:10.5f}  {r['dec']:10.5f}  {r['pa']:8.2f}  "
            f"{c_strs}  {r['total']:6d}  {cams_ok_str:>14}"
        )
        table_lines.append(line)
        
    table_lines.append(double_separator)
    report_text = "\n".join(table_lines)
    print(report_text)
    
    print(f"\nSummary:")
    print(f"Total Pointings: {len(t_in)}")
    print(f"Initially Valid or Adjusted Successfully: {success_count} (Adjusted: {adjusted_count})")
    print(f"Failed to find valid pointing nearby: {fail_count}")
    
    if args.report_file:
        with open(args.report_file, "w") as f:
            f.write(report_text + f"\n\nTotal: {len(t_in)}, Success: {success_count}, Adjusted: {adjusted_count}, Failed: {fail_count}\n")
        print(f"Saved guide star report to {args.report_file}")
        
    if args.add_columns:
        for cam_idx in range(6):
            t_in[f"ag{cam_idx}"] = [r["counts"][cam_idx] for r in report_rows]
        t_in["n_guidestars"] = [r["total"] for r in report_rows]
        
    t_in.write(args.output, format="ascii.ecsv", overwrite=True)
    print(f"Saved optimized pointings to {args.output}")

if __name__ == "__main__":
    main()
