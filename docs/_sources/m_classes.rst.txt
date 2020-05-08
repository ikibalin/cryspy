Built-in classes
===================

CrysPy provides several built-in objects describing crystal structure or 
polarized neutron diffraction experiments.

Objects to describe the crystal structure:

   * Cell
   * SpaceGroup
   * AtomSite
   * AtomSiteAniso
   * Crystal

Objects to describe the magnetic crystal structure:

   * AtomSiteMoment
   * AtomSiteSusceptibility
   * AtomTypeScat
   * AtomSiteScat


Objects describing incident beam of polarized neutrons

   * DiffrnRadiation
   * Setup

Objects describing single-crystal diffraction experiment: 

   * DiffrnRefln
   * DiffrnOrientMatrix
   * Extinction
   * Diffrn

Objects describing powder diffraction experiment: 

   * Exclude
   * Range
   * Texture
   * Pd
   * Pd2d
   * PdInstrResolution
   * PdInstrReflexAsymmetry
   * PdBackground
   * PdMeas
   * Pd2dMeas
   * Pd2dBackground
   * Pd2dInstrReflexAsymmetry
   * Pd2dInstrResolution
   * Pd2dPeak

Objects used for data refinement: 
   
   * RhoChi

For further information see:

.. toctree::
   :maxdepth: 1

   cryspy.corecif
   cryspy.symcif
   cryspy.magneticcif
   cryspy.cif_like
   cryspy.pd1dcif_like
   cryspy.pd2dcif_like
   cryspy.scripts
