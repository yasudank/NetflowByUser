import numpy as np
import pprint
from collections import defaultdict
from astropy.table import Table
import ets_fiber_assigner.netflow as nf

def compute_costs(priority_list):
    # P0=3^10, P1=3^9, ..., P9=3
    return [3 ** (10 - p) for p in priority_list]

def cobraMoveCost(dist):
    return 0.*dist

def solve_assignment(bench, tgt, telescopes, pipe_config):
    """Solve the fiber assignment netflow problem.
    
    Returns
    -------
    res : list of dict
        Assigned target-to-cobra mappings per exposure.
    tpos : list of ndarray
        Focal plane positions of all targets per exposure.
    exposures_data : list of dict
        Extracted assignment data for saving.
    all_classes : set
        Set of targeted object classes.
    stats_per_exp : list
        List of statistics (index, count, class dict) per exposure.
    """
    print("Getting focal plane positions")
    tpos = []
    for tele in telescopes:
        fp_pos = tele.get_fp_positions(tgt)
        # Calculate angular distance to the telescope pointing center
        cos_dec = np.cos(np.radians(tele._dec))
        tgt_ra = np.array([t.ra for t in tgt])
        tgt_dec = np.array([t.dec for t in tgt])
        
        dx = (tgt_ra - tele._ra) * cos_dec
        dy = tgt_dec - tele._dec
        dist = np.sqrt(dx**2 + dy**2)
        
        # Filter out targets that are too far (dist > 0.8 degrees) by setting to a dummy value
        fp_pos[dist > 0.8] = 9999.0 + 9999.0j
        tpos.append(fp_pos)
    print("Done")

    # Create the dictionary containing the costs and constraints for all classes of targets
    classdict = {}
    priority_list = np.unique([t._pri for t in tgt if t.targetclass.startswith('sci')])
    priority_list.sort()
    print(f"Target priorities (in total {len(priority_list)}): {priority_list}.")

    costs = compute_costs(priority_list)

    for x, cost in zip(priority_list, costs):
        classdict[f"sci_P{x}"] = {
            "nonObservationCost": cost,
            "partialObservationCost": cost,
            "calib": False,
        }

    classdict["sky"] = {
        "numRequired": pipe_config["netflow"]["sky"]["num_required"],
        "nonObservationCost": pipe_config["netflow"]["sky"]["non_observation_cost"],
        "calib": True
    }
    classdict["cal"] = {
        "numRequired": pipe_config["netflow"]["fluxstd"]["num_required"],
        "nonObservationCost": pipe_config["netflow"]["fluxstd"]["non_observation_cost"],
        "calib": True
    }

    pprint.pprint(classdict)

    # Optional: slightly increase the cost for later observations, to observe as early as possible
    vis_cost = [i*0. for i in range(len(telescopes))]

    # Duration of one observation in seconds
    t_obs = pipe_config["netflow"]["t_obs"]

    # 露出時間を退避して t_obs に上書き
    orig_obs_times = {}
    for t in tgt:
        if hasattr(t, "targetclass") and t.targetclass.startswith("sci"):
            orig_obs_times[t.ID] = getattr(t, "_obs_time", None)
            t._obs_time = float(t_obs)

    try:
        # Gurobi configuration from pipeline config
        gurobi_config = pipe_config['gurobi']
        gurobiOptions = dict(seed=gurobi_config['seed'], 
                             presolve=gurobi_config['presolve'], 
                             method=gurobi_config['method'], 
                             degenmoves=gurobi_config['degenmoves'],
                             heuristics=gurobi_config['heuristics'], 
                             mipfocus=gurobi_config['mipfocus'], 
                             mipgap=gurobi_config['mipgap'],
                             PreSOS2Encoding=gurobi_config['PreSOS2Encoding'],
                             PreSOS1Encoding=gurobi_config['PreSOS1Encoding'],
                             threads=gurobi_config['threads'])

        alreadyObserved={}

        # Retrieve PFS software parameters from config
        pfs_config = pipe_config.get("pfs", {})
        
        # brokenCobrasMargin
        broken_cobras_margin = pfs_config.get("brokenCobrasMargin")
        if broken_cobras_margin is None:
            broken_cobras_margin = pipe_config["netflow"].get("broken_cobras_margin", 1.0)
            
        # fiducialsAvoidDistance
        fiducials_avoid_distance = pfs_config.get("fiducialsAvoidDistance", 0.0)
        
        # dot_penalty -> blackDotPenalty
        black_dot_penalty = pfs_config.get("dot_penalty")

        # cobraSafetyMargin
        cobra_safety_margin = pfs_config.get("cobraSafetyMargin")
        if cobra_safety_margin is None:
            cobra_safety_margin = pfs_config.get("cobra_safety_margin", 0.0)

        # numReservedFibers
        num_reserved_fibers = pfs_config.get("numReservedFibers")
        if num_reserved_fibers is None:
            num_reserved_fibers = pfs_config.get("num_reserved_fibers", 0)

        # fiberNonAllocationCost
        fiber_non_allocation_cost = pfs_config.get("fiberNonAllocationCost")
        if fiber_non_allocation_cost is None:
            fiber_non_allocation_cost = pfs_config.get("fiber_non_allocation_cost", 0.0)

        # SFA-like per-location constraints for sky fibers to avoid validation warnings
        # The validation tool uses 20 specific geometric sectors (Laszlo regions).
        # We must group the cobras exactly according to these sectors.
        from pfs_obsproc_planning.utils import plot_pfsDesign as pldes
        import pandas as pd
        
        df_fib_mock = pd.DataFrame({
            'pfi_x': bench.cobras.centers.real,
            'pfi_y': bench.cobras.centers.imag
        })
        cobraRegions = pldes.get_field_sector2(df_fib_mock)

        # Compute observation strategy
        penalty = pipe_config["netflow"].get("locationGroupPenalty", 1e11)
        prob = nf.buildProblem(bench, tgt, tpos, classdict, t_obs,
                               vis_cost, cobraMoveCost=cobraMoveCost,
                               collision_distance=pipe_config["netflow"]["collision_distance"],
                               elbow_collisions=pipe_config["netflow"]["elbow_collisions"],
                               gurobi=True, gurobiOptions=gurobiOptions,
                               alreadyObserved=alreadyObserved,
                               brokenCobrasMargin=broken_cobras_margin,
                               fiducialsAvoidDistance=fiducials_avoid_distance,
                               blackDotPenalty=black_dot_penalty,
                               cobraSafetyMargin=cobra_safety_margin,
                               numReservedFibers=num_reserved_fibers,
                               fiberNonAllocationCost=fiber_non_allocation_cost,
                               cobraLocationGroup=cobraRegions,
                               minSkyTargetsPerLocation=12,
                               locationGroupPenalty=1e11)

        print("solving the problem")
        prob.solve()

        # Extract solution
        res = [{} for _ in range(len(telescopes))]
        for k1, v1 in prob._vardict.items():
            if k1.startswith("Tv_Cv_"):
                visited = prob.value(v1) > 0
                if visited:
                    _, _, tidx, cidx, ivis = k1.split("_")
                    res[int(ivis)][int(tidx)] = int(cidx)

        # Extract and format exposures data
        all_classes = set()
        stats_per_exp = []
        exposures_data = []

        pointing_file = pipe_config["inputs"]["pointing_file"]
        if pointing_file is None:
            pointing_file = "optimized_pointings.ecsv"
        ppcList = Table.read(pointing_file, format="ascii.ecsv")
        ppc_codes = ppcList['ppc_code'].tolist()

        with open(pipe_config["outputs"]["text_output"], "w") as f:
            for i, (vis, tp, tel) in enumerate(zip(res, tpos, telescopes)):
                tdict = defaultdict(int)
                ppc_code = ppc_codes[i] if i < len(ppc_codes) else f"EXP_{i+1}"
                
                f.write("# Exposure {}: duration {}s, RA: {}, Dec: {}, PA: {}\n".
                        format(i+1, t_obs, tel._ra, tel._dec, tel._posang))
                f.write("# Target    Fiber          X          Y         RA        DEC\n")
                for tidx, cidx in vis.items():
                    cls = tgt[tidx].targetclass
                    tdict[cls] += 1
                    all_classes.add(cls)
                    f.write("{:} {:6d} {:10.5f} {:10.5f} {:10.5f} {:10.5f}\n"
                            .format(tgt[tidx].ID, cidx+1, tp[tidx].real, tp[tidx].imag,
                                    tgt[tidx].ra, tgt[tidx].dec))
                    
                    exposures_data.append({
                        'ppc_code': ppc_code,
                        'target_id': str(tgt[tidx].ID),
                        'cobraId': cidx + 1,
                        'pfi_X': tp[tidx].real,
                        'pfi_Y': tp[tidx].imag
                    })
                stats_per_exp.append((i, len(vis), tdict))

        return res, tpos, exposures_data, all_classes, stats_per_exp

    finally:
        # 露出時間を元の値に復帰
        for t in tgt:
            if hasattr(t, "targetclass") and t.targetclass.startswith("sci"):
                if t.ID in orig_obs_times and orig_obs_times[t.ID] is not None:
                    t._obs_time = orig_obs_times[t.ID]

