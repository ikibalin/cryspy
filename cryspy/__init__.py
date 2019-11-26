"""
CrysPy
======

To run rhochi refinement from a command-line:
>>> python -m cryspy.run_rhochi

To create a template in the folder run from a command-line:
>>> python -m cryspy.rhochi_diffrn # template for single-crystal diffraction experiment
>>> python -m cryspy.rhochi_pd # template for 1d powder diffraction experiment
>>> python -m cryspy.rhochi_pd2d # template for 2d powder diffraction experiment
>>> python -m cryspy.rhochi_pd2dt # template for 2d powder diffraction with one-axial texture model experiment

Objects to describe crystal structure
---------------------------------------
- Crystal
    the main object to describe crystal structure
    - SpaceGroup
    - Cell
    - AtomType
    - AtomSite
    - AtomSiteMagnetism
    - AtomSiteMagnetismAniso
    - AtomSiteAniso

Objects to describe experiments
---------------------------------------
- BeamPolarization
    common object for all polaraized neutron diffraction experiments
- Diffrn
    the main object to describe single crystal diffraction experiment
    - OrientMatrix 
    - Extinction 
    - DiffrnRefln 
- Pd
    the main object to describe 1d polaraized neutrhon powder diffraction experiment
    - PdMeas
        measurements
    - PdBackground
    - PdExclude2Theta
    - PdInstrReflexAsymmetry
    - PdInstrResolution
    - PdPhase
    - PdPeak
    - PdProc
- Pd2d
    the main object to describe 2d polaraized neutrhon powder diffraction experiment
    - Pd2dMeas
        measurements
    - Pd2dBackground
    - Pd2dExclude2Theta
    - Pd2dInstrReflexAsymmetry
    - Pd2dInstrResolution
    - Pd2dPhase
    - Pd2dPeak
    - Pd2dProc
- Pd2dt
    the main object to describe 2d polaraized neutrhon powder diffraction experiment in case of one-axial texture
    
Use help function to get information for each of them:
>>> import cryspy
>>> help(cryspy.Cell)
    class Cell(builtins.object)
     |  Cell(a=1.0, b=1.0, c=1.0, alpha=90.0, beta=90.0, gamma=90.0, bravais_lattice='triclinic')
     |
     |  Data items in the Cell class record details about the
     |  crystallographic cell parameters and their measurement.
     ...
     
Specify the input parameters of the object to get detailed description
>>> help(cryspy.Cell.bravais_lattice)
"""
name = "cryspy"

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
    AtomSiteAniso,
    AtomSiteMoment
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
    Pd2dExclude,
    Pd2dInstrReflexAsymmetry,
    Pd2dInstrResolution,
    Pd2dMeas,
    Pd2dPhase,
    Pd2dPeak,
    Pd2dProc,
    Pd2d,
    Pd2dt
    )

#from .scripts.rhochi.rhochi_viewer import main as rhochi_gui

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

