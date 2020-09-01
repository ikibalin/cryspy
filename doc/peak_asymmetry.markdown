# Peak asymmetry

The default asymmetry correction is a multiplier term to the peak shape. The asymmetry correction
adopts the form:

$$
A_s(z) = 1+ \fract{P_1 F_a(z)+P_2 F_b(z)}{tanh \theta_h} + \fract{P_3 F_a(z)+P_4 F_b(z)}{tanh 2\theta_h}
$$

where

$$
z = \fract{2\theta_i - 2\theta_h-S_{shf}}{FMHM}
$$

$S_{shf}$ includes the zero-point and other shifting terms,

$$
F_a(z) = 2 z \exp\left( -z2 \right)
$$

$$
F_b(z) = 2 (2 z^2-3) F_a(z)
$$

The asymmetry correction has four independent parameters $P_1$, $P_2$, $P_3$, $P_4$.