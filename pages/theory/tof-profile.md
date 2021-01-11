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

Rietveld refinement with spallation Neutron Powder Diffraction Data
(according J.Appl.Cryst.(1982). 15, 581-589)

# Peak shape function

For Gauss profile the peak shape function is calculated as:

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


For pseudo-Voigt profile the Lorentz part of the peak shape function is calculated as:

$$ \sigma^2 = (\sigma_2^2+\text{size}_G) \cdot d^4 + (\sigma_1^2+\text{strain}_G) \cdot d^2 + \sigma_0^2 $$

$$ \gamma = (\gamma_2+\text{size}_L) \cdot d^2 + (\gamma_1+\text{strain}_L) \cdot d + \gamma_0 $$

$$ z_1 = \alpha \cdot \Delta + 0.5j \cdot \alpha \cdot \gamma $$

$$ z_2 = -\beta \cdot \Delta + 0.5j  \cdot \beta \cdot \gamma $$

$$ fz_1 = \text{exp1}(z_1) $$

$$ fz_2 = \text{exp1}(z_2) $$

$$ \Omega_a = -2 \pi \cdot \text{Im}(fz_1) $$

$$ \Omega_b = -2 \pi \cdot \text{Im}(fz_2) $$

$$ \text{profile}_L = \frac{1}{2} \frac{\alpha \beta}{\alpha + \beta} \cdot (\Omega_a + \Omega_b) $$

$$ \text{profile}_{pV} = (1 - \eta) \cdot \text{profile}_G + \eta \cdot \text{profile}_L $$



The parameters $$\alpha_0, \alpha_1, \beta_0, \beta_1, \sigma_0, \sigma_1, \gamma_0, \gamma_1$$ are defined
in [tof profile][tof-profile]. The parameter **peak_shape** defines the shape of the peak: "Gauss" or "pseudo-Voigt".

[expressions]: /cryspy/theory
[tof-profile]: /cryspy/rcif-format/tof-profile
