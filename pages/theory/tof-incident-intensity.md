---
layout: page
permalink: /theory/tof-incident-intensity/
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



[expressions]: /cryspy/theory
[tof-intensity-incident]: /cryspy/rcif-format/tof-intensity-incident
