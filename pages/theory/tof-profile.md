---
layout: page
permalink: /theory/tof-profile/
---
<html>
<!-- Mathjax Support -->
  <head>
    <meta charset="utf-8">
    <title>{{ page.title }}</title>
<script type="text/javascript" async
  src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-MML-AM_CHTML">
</script>    
  </head>
</html>

[Back to expressions][expressions]

# The Gauss peak shape function

$$
\text{profile}_G = \frac{\alpha \beta}{\alpha + \beta} \left[ \exp(u) \cdot \text{erfc}(y) + \exp(v) \cdot \text{erfc}(z) \right]
$$

$$ u = \frac{\alpha}{2} \left( \alpha \sigma^2 + 2 \Delta \right) $$

$$ v = \frac{\beta}{2} \left( \beta \sigma^2 - 2 \Delta \right) $$

$$ y = \frac{\alpha \sigma^2 + \Delta}{\sqrt{2} \sigma} $$

$$ z = \frac{\beta \sigma^2 + \Delta}{\sqrt{2} \sigma} $$

$$ \Delta = t - t_{hkl} $$

Peak-shape TOF dependence is:

$$ \alpha = \alpha_0 + \alpha_1 / d $$

$$ \beta = \beta_0 + \beta_1 / d^4 $$

$$ \sigma = \sigma_0 + \sigma_1 \cdot d $$

The above expressions are taken from R.B. Von Dreele, J.D. Jorgensen and C.G. Windsor, "Rietveld refinement with spallation Neutron Powder Diffraction Data", J.Appl.Cryst.(1982). 15, 581-589.

# The pseudo-Voigt peak shape function

$$ \text{profile}_{pV} = (1 - \eta) \cdot \text{profile}_G + \eta \cdot \text{profile}_L $$

The calculation of the Gauss part is given above. The  Lorentz part of the pseudo-Voigt peak shape function is given as:

$$ \text{profile}_L = \frac{1}{2} \frac{\alpha \beta}{\alpha + \beta} \cdot (\Omega_a + \Omega_b) $$

The expressions for $$\alpha$$ and $$\beta$$ parameters are given above. The other parameters are defined as:

$$ \Omega_a = -2 \pi \cdot \text{Im}(fz_1) $$

$$ \Omega_b = -2 \pi \cdot \text{Im}(fz_2) $$

in which

$$ fz_1 = \text{exp1}(z_1) $$

$$ fz_2 = \text{exp1}(z_2) $$

where $$z_1, z_2$$:

$$ z_1 = \alpha \cdot \Delta + 0.5j \cdot \alpha \cdot \gamma $$

$$ z_2 = -\beta \cdot \Delta + 0.5j  \cdot \beta \cdot \gamma $$

The parameter $$\Delta$$ is defined above. The $$\gamma$$ value is calculated as

$$ \gamma = (\gamma_2+\text{size}_L) \cdot d^2 + (\gamma_1+\text{strain}_L) \cdot d + \gamma_0 $$


The parameter $$\sigma$$ which is defined above in the case of the Gauss shape profile functions is replaced by next one:

$$ \sigma^2 = (\sigma_2^2+\text{size}_G) \cdot d^4 + (\sigma_1^2+\text{strain}_G) \cdot d^2 + \sigma_0^2 $$


The parameters $$\alpha_0, \alpha_1, \beta_0, \beta_1, \sigma_0, \sigma_1, \sigma_2, \gamma_0, \gamma_1, \gamma_2$$ are defined
in [tof profile][tof-profile]. The parameter **peak_shape** defines the shape of the peak: "Gauss" or "pseudo-Voigt".

[expressions]: /cryspy/theory
[tof-profile]: /cryspy/rcif-format/tof-profile
