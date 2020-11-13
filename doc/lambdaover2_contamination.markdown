# lambda/2 contamination

For experiments with single crystal the correction of lambda/2 contamination is
introduced by a nonzero coefficient in front of unpolarized integrated
intensity with doubled Miller indeces:

$$
Int_{2h2k2l} = \frac{1}{2} \cdot ( Int^{UP}_{2h2k2l} + Int^{DOWN}_{2h2k2l} )
$$

$$
FlipRatio = (Int^{UP}_{hkl} + coeff * Int_{2h2k2l}) /
            (Int^{DOWN}_{hkl} + coeff * Int_{2h2k2l})
$$

In the calculations of  the extinction correction y for  2h2k2l reflections the
half of the wavelength is used.
