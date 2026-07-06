#!/usr/bin/env python
"""Check sky distribution across Laszlo 20-sectors using self-contained logic."""
import os, yaml
import numpy as np
import pandas as pd
from astropy.table import Table
import matplotlib
matplotlib.use('Agg')  # non-interactive backend

import netflow_io
import netflow_instrument

# Reproduce Laszlo sector logic from plot_pfsDesign.get_field_sector2
def get_sector_reference_xy():
    points = []
    points.append([[0.0, 0.0]])
    phi = np.radians(np.linspace(0, 360, 6, endpoint=False) + 15 + 90)
    points.append(100.0 * np.stack([np.cos(phi), np.sin(phi)], axis=-1))
    phi = np.radians(np.linspace(0, 360, 13, endpoint=False) - 20 + 90)
    points.append(175.0 * np.stack([np.cos(phi), np.sin(phi)], axis=-1))
    return np.concatenate(points, axis=0)

def get_field_sector2(pfi_x, pfi_y):
    xy = get_sector_reference_xy()
    uv = np.stack([pfi_x, pfi_y], axis=-1)
    diff = uv[None, :, :] - xy[:, None, :]
    d2 = np.sum(diff * diff, axis=2)
    tag = np.argmin(d2, axis=0)
    return tag

# Load config
config_file = "netflow_pipeline_config.yaml"
with open(config_file, "r") as f:
    pipe_config = yaml.safe_load(f)
if os.path.exists("config.yaml"):
    with open("config.yaml", "r") as f:
        config_yaml = yaml.safe_load(f)
        pipe_config["gurobi"] = config_yaml.get("gurobi", {}).get("param", pipe_config["gurobi"])
        pipe_config["pfs"] = config_yaml.get("pfs", {})

# Load bench
black_dot_radius_margin = pipe_config.get("pfs", {}).get("black_dot_radius_margin", 1.65)
bench = netflow_instrument.getBench(black_dot_radius_margin)

# Cobra sectors
cobra_sectors = get_field_sector2(bench.cobras.centers.real, bench.cobras.centers.imag)

print(f"\n=== Cobra sector distribution (20 Laszlo sectors) ===")
unique_sectors = np.arange(20)
for s in unique_sectors:
    cnt = np.sum(cobra_sectors == s)
    print(f"  Sector {s:2d}: {cnt:3d} cobras")

# Also compute what the solver's cobraLocationGroup looks like
ncobras = bench.cobras.nCobras
solverRegions = np.zeros(ncobras, dtype=np.int32)
solverRegions_ = np.array_split(solverRegions, 20)
for i in range(20):
    solverRegions_[i] += i
solverRegions = np.concatenate(solverRegions_)

print(f"\n=== Solver cobraLocationGroup (get_field_sector2) ===")
# Check: are they the same?
match = np.sum(cobra_sectors == solverRegions)
print(f"  Matches between Laszlo sectors and solver regions: {match}/{ncobras}")
if match != ncobras:
    print("  WARNING: Solver regions do NOT match Laszlo sectors!")
    # Show mismatches
    mismatches = np.where(cobra_sectors != solverRegions)[0]
    print(f"  First 10 mismatches (cobra_idx: laszlo_sector vs solver_region):")
    for m in mismatches[:10]:
        print(f"    cobra {m}: laszlo={cobra_sectors[m]} vs solver={solverRegions[m]}")

# Check the solver's actual current cobraRegions (get_field_sector2-based)
print(f"\n=== Solver cobraRegions (from get_field_sector2 on cobras) ===")
solver_sectors_geo = get_field_sector2(bench.cobras.centers.real, bench.cobras.centers.imag)
for s in unique_sectors:
    cnt = np.sum(solver_sectors_geo == s)
    print(f"  Sector {s:2d}: {cnt:3d} cobras")

# Now check saved sky assignments
pointing_file = "optimized_pointings.ecsv"
pointings = Table.read(pointing_file, format="ascii.ecsv")

for r_ptg in pointings:
    ppc_code = r_ptg["ppc_code"]
    sky_path = os.path.join("targets", "sky", f"{ppc_code}.ecsv")
    if not os.path.exists(sky_path):
        print(f"\n{ppc_code}: sky file not found")
        continue
    
    df_sky = Table.read(sky_path, format="ascii.ecsv").to_pandas()
    if "cobraId" not in df_sky.columns:
        print(f"\n{ppc_code}: no cobraId column")
        continue
    
    cobra_ids = df_sky["cobraId"].values - 1  # 0-based
    sky_sectors = cobra_sectors[cobra_ids]
    sky_counts = np.bincount(sky_sectors, minlength=20)
    
    arr = sky_counts[sky_counts > 0]
    
    print(f"\n=== {ppc_code}: {len(df_sky)} sky fibers ===")
    print(f"  Sectors with sky > 0: {len(arr)} / 20")
    if len(arr) >= 12:
        print(f"  sky_min={np.min(arr)}, sky_max={np.max(arr)}")
    else:
        print(f"  WARNING: fewer than 12 sectors have sky → sky_min forced to 0")
    
    for s in unique_sectors:
        cnt = sky_counts[s]
        marker = " <-- ZERO!" if cnt == 0 else (" <-- LOW" if cnt < 12 else "")
        print(f"  Sector {s:2d}: {cnt:3d} sky{marker}")
    
    # Also check fluxstd
    cal_path = os.path.join("targets", "fluxstd", f"{ppc_code}.ecsv")
    if os.path.exists(cal_path):
        df_cal = Table.read(cal_path, format="ascii.ecsv").to_pandas()
        if "cobraId" in df_cal.columns:
            cal_ids = df_cal["cobraId"].values - 1
            cal_sectors = cobra_sectors[cal_ids]
            cal_counts = np.bincount(cal_sectors, minlength=20)
            cal_arr = cal_counts[cal_counts > 0]
            print(f"  --- fluxstd: {len(df_cal)} total ---")
            print(f"  Sectors with std > 0: {len(cal_arr)} / 20")
            if len(cal_arr) >= 12:
                print(f"  std_min={np.min(cal_arr)}, std_max={np.max(cal_arr)}")
            else:
                print(f"  WARNING: fewer than 12 sectors have std → std_min forced to 0")
