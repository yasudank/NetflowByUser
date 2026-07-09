import os
import numpy as np
from astropy.table import Table
import pfs.instdata
import ets_fiber_assigner
import ets_fiber_assigner.netflow as nf

def getBench(black_dot_radius_margin=1.65):
    from ics.cobraOps.Bench import Bench
    from ics.cobraCharmer.cobraCoach.cobraCoach import CobraCoach
    pfs_instdata_dir = os.path.dirname(pfs.instdata.__file__)
    os.environ["PFS_INSTDATA_DIR"] = pfs_instdata_dir
    ets_fiber_assigner_dir = os.path.dirname(ets_fiber_assigner.__file__)
    cobraCoach = CobraCoach(
        loadModel=True, trajectoryMode=True, rootDir=ets_fiber_assigner_dir)
    if black_dot_radius_margin is None:
        black_dot_radius_margin = 1.0
    bench = Bench(cobraCoach, blackDotsMargin=black_dot_radius_margin)
    print("Number of cobras:", bench.cobras.nCobras)
    return bench

def getPointingCenter(file, nvisit, posang=0.0, otime="2020-01-01"):
    ppcList = Table.read(file, format="ascii.ecsv")
    raTel = ppcList['ppc_ra']
    decTel = ppcList['ppc_dec']
    if 'ppc_pa' in ppcList.columns:
        paTel = ppcList['ppc_pa']
    else:
        paTel = np.zeros(len(ppcList)) + posang

    telescopes = []
    for _ in range(nvisit):
        for _ra, _dec, _pa in zip(raTel, decTel, paTel):
            telescopes.append(nf.Telescope(_ra, _dec, _pa, otime))

    return telescopes
