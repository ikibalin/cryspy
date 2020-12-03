---
layout: page
permalink: /theory/preferred-orientation/
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

[Back to notes][notes]

# Preferred orientation

Modified March's function:

$$ P_h = G_2 + \left(1-G_2\right) \cdot \left(\left(G_1 \cos \alpha_h\right)^2 + \frac{\sin^2 \alpha_h} {G_1}\right)^{-3/2} $$

Where $$G_1$$ and $$G_2$$ are refinable parameters and is the acute angle between the scattering vector and the normal to the crystallites (platy habit).

Description of texture in cif file:

    _texture_g_1 0.1239
    _texture_g_2 0.94211
    _texture_h_ax -0.66119
    _texture_k_ax -0.0541
    _texture_l_ax 3.0613
    
[notes]: /cryspy/theory