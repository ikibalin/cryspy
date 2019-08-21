# Work with sublibrary experiment of library NeuPy

from neupy import (ResolutionPowder1D,
                    AsymmetryPowder1D,
                    FactorLorentzPowder1D,
                    BackgroundPowder1D,
                    BeamPolarization,
                    SetupPowder1D,
                    ExperimentPowder1D)



resolution_powder_1d = ResolutionPowder1D(u=16, v=-5, w=0.5)
asymmetry_powder_1d = AsymmetryPowder1D()
factor_lorentz = FactorLorentzPowder1D()
background_powder_1d = BackgroundPowder1D()
beam_polarization = BeamPolarization(p_u=0.95, flipper_efficiency=1.0)

setup_powder_1d = SetupPowder1D(label="powder_1d",wave_length=1.4, zero_shift=0.1, 
    resolution=resolution_powder_1d,factor_lorentz=factor_lorentz, asymmetry=asymmetry_powder_1d,
    beam_polarization=beam_polarization, background=background_powder_1d)


experiment_powder_1d = ExperimentPowder1D(name="powder_1d", setup=setup_powder_1d)




