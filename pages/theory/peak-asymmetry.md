---
layout: page
permalink: /theory/peak-asymmetry/
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

[Back to expressions][notes]

# Peak asymmetry

The default asymmetry correction is a multiplier term to the peak shape. The asymmetry correction adopts the form:

$$ A_s(z) = 1+ \frac{P_1 F_a(z)+P_2 F_b(z)}{\tanh \theta_h} + \frac{P_3 F_a(z)+P_4 F_b(z)}{\tanh 2\theta_h} $$

where

$$ z = \frac{2\theta_i - 2\theta_h-S_{shf}}{FMHM} $$

$$S_{shf}$$ includes the zero-point and other shifting terms,

$$ F_a(z) = 2 z \cdot \exp\left( -z^2 \right) $$

$$ F_b(z) = 2 (2 z^2-3) F_a(z) $$

The asymmetry correction has four independent parameters $$P_1$$, $$P_2$$, $$P_3$$, $$P_4$$.


[notes]: /cryspy/theory