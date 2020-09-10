# lambda/2 contamination

For experiments with single crystal the correction of lambda/2 contamination is introduced by coeff. in equation for flip ratio:

$$
FlipRatio = (Int^{UP}_{hkl} + coeff * Int^{UP}_{2h2k2l}) / (Int^{DOWN}_{hkl} + coeff * Int^{DOWN}_{2h2k2l})
$$

In the calculations of  the extinction correction y for  2h2k2l reflections the half of the wavelength is used.