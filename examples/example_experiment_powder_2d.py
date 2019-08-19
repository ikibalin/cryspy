# Work with sublibrary experiment of library NeuPy

from libNeuPy.f_experiment import ResolutionPowder2D

resolution_powder_2d = ResolutionPowder2D(u=16, v=-5, w=0.5)


from libNeuPy.f_experiment import AsymmetryPowder2D

asymmetry_powder_2d = AsymmetryPowder2D()

from libNeuPy.f_experiment import FactorLorentzPowder2D

factor_lorentz = FactorLorentzPowder2D()


from libNeuPy.f_experiment import BackgroundPowder2D

background_powder_2d = BackgroundPowder2D()


from libNeuPy.f_experiment import BeamPolarization

beam_polarization = BeamPolarization(p_u=0.95, flipper_efficiency=1.0)


from libNeuPy.f_experiment import SetupPowder2D

setup_powder_2d = SetupPowder2D(label="powder_2d",wave_length=1.4, zero_shift=0.1, 
    resolution=resolution_powder_2d,factor_lorentz=factor_lorentz, asymmetry=asymmetry_powder_2d,
    beam_polarization=beam_polarization, background=background_powder_2d)


from libNeuPy.f_experiment import ExperimentPowder2D

experiment_powder_2d = ExperimentPowder2D(name="powder_2d", setup=setup_powder_2d)




