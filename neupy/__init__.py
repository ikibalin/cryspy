from .f_common.cl_fitable import Fitable

from .scripts.rhochi.cl_rhochi import RhoChi, rhochi_refinement, rhochi_read_file

from .f_crystal import (
    SpaceGroup,
    Cell,
    AtomType,
    AtomSite,
    Crystal,
    calc_mRmCmRT,
    AtomSiteMagnetism,
    AtomSiteMagnetismAniso,
    AtomSiteAniso
    )

from .f_experiment import (
    BeamPolarization,
    OrientMatrix,
    Extinction,
    DiffrnRefln, 
    Diffrn,
    PdBackground,
    PdExclude2Theta,
    PdInstrReflexAsymmetry,
    PdInstrResolution,
    PdMeas,
    PdPeak,
    PdPhase,
    PdProc,
    Pd,
    Pd2dBackground,
    Pd2dExclude2Theta,
    Pd2dInstrReflexAsymmetry,
    Pd2dInstrResolution,
    Pd2dMeas,
    Pd2dPhase,
    Pd2dPeak,
    Pd2dProc,
    Pd2d,
    Pd2dt
    )


def crystal_from_file(f_name):
    rhochi = RhoChi()
    with open(f_name, 'r') as fid:
        string = fid.read()
    rhochi.from_cif(string)
    if len(rhochi.crystals) != 0:
        crystal = rhochi.crystals[0]
    else:
        crystal = Crystal()
        print("Information about crystal is not found in the file: '{:}'".format(f_name))
    return crystal
