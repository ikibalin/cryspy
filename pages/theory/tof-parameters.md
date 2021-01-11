---
layout: page
permalink: /theory/tof-parameters/
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


[expressions]: /cryspy/theory

[tof-parameters]: /cryspy/rcif-format/tof-parameters
