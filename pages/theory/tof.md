---
layout: page
permalink: /theory/tof/
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

# Intensity function

$$
I_c(t) = k I_s(t) y_B(t) + y_b(t)
$$

$$I_s$$ represent empirically the effective spectrum by exponential series:

$$
I_s(t) = A_0 + A_1 \exp(-A_2 t) +
 A_3 \exp(-A_4 t^2) + A_5 \exp(-A_6 t^3) + A_7 \exp(-A_8 t^4)
$$

or Maxwellian spectrum:

$$
I_s(t) = A_0 + \frac{A_1}{t^5} \exp(-\frac{A_2}{t^2}) +
 A_3 \exp(-A_4 t^2) + A_5 \exp(-A_6 t^3) + A_7 \exp(-A_8 t^4)
$$

The parameters $$A_{1, .., 8}$$ are defined in [TOF Intensity Incident][tof-intensity-incident]. 
The parameter **spectrum** defines the spectrum type: "Empirical-Exponents" or "Maxwell".


# TOF and d spacing

TOF of a Bragg peak on a given instrument is proportional to the d spacing and depends from the neutron spectra. 
For thermal neutrons:
    
$$ t = \text{zero} + \text{dtt1} \cdot d + \text{dtt2} \cdot d^2. $$

For epithermal neutrons:

$$ t_e = \text{zero} + \text{dtt1} \cdot d $$    

$$ t_t = \text{zerot} + \text{dtt1t} \cdot d - \text{dtt2t} / d $$    
        
$$ n_{\text{cross}} = 0.5 \cdot \text{erfc}(\text{width} \cdot (x_{\text{cross}} - 1/d)) $$    
        
$$ t = n_{\text{cross}} \cdot t_e + (1-n_{\text{cross}}) \cdot t_t $$    

The parameters $$\text{zero}, \text{dtt1}, \text{dtt2}, \text{zerot}, \text{dtt1t}, \text{dtt2t}, \text{width}, x_{\text{cross}}, n_{\text{cross}}$$ are defined
in [tof parameters][tof-parameters]. The parameter **neutrons** defines the neutron spectra: "thermal" or "epithermal".

# Extinction

During the analysis it was observed that the reflections with the largest d spacings displayed considerable extinction. To account for this effect a simple correction factor was applied to the calculated Bragg intensities:

$$
I_i^{'} = I_i \cdot (1-E \cdot I_i m_i d_{i}^{4})
$$

The extinction parameter $$E$$ is defined in [tof parameters][tof-parameters] as extinction.

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


# Background $$y_b(t)$$

Background is described by the series of coefficients in front of cosines functions:

$$ y_b = \sum_{i=1}^{18} \text{coeff}_{i} \cdot \cos \left(\pi \cdot i \cdot t / t_{\max} \right)  $$

The parameters $$\text{coeff}_{1, .., 18}, t_{\max}$$ are defined
in [tof background][tof-background]. 


[expressions]: /cryspy/theory

[tof-parameters]: /cryspy/rcif-format/tof-parameters
[tof-intensity-incident]: /cryspy/rcif-format/tof-intensity-incident
[tof-profile]: /cryspy/rcif-format/tof-profile
[tof-background]: /cryspy/rcif-format/tof-background