---
layout: page
permalink: /theory/tof-background/
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

# Background $$y_b(t)$$

Background is described by the series of coefficients in front of cosines functions:

$$ y_b = \sum_{i=1}^{18} \text{coeff}_{i} \cdot \cos \left(\pi \cdot i \cdot t / t_{\max} \right)  $$

The parameters $$\text{coeff}_{1, .., 18}, t_{\max}$$ are defined
in [tof background][tof-background]. 


[expressions]: /cryspy/theory
[tof-background]: /cryspy/rcif-format/tof-background
