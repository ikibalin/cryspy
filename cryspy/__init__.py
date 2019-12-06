"""
CrysPy Library
=================

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
from .cl_crystal import (
    Crystal, 
    Cell, 
    SpaceGroup, 
    Fitable, 
    ItemConstr, 
    LoopConstr, 
    DataConstr)

from .symcif.cl_space_group_wyckoff import SpaceGroupWyckoff, SpaceGroupWyckoffL
from cryspy.symcif.cl_space_group_symop import SpaceGroupSymop, SpaceGroupSymopL
    
from .corecif.cl_atom_site import AtomSite, AtomSiteL