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

[Back to expressions][notes]

# Preferred orientation

Modified March's function:

$$ P_h = G_2 + \left(1-G_2\right) \cdot \left(\left(G_1 \cos \alpha_h\right)^2 + \frac{\sin^2 \Delta\alpha} {G_1}\right)^{-3/2} $$

Where $$G_1$$ and $$G_2$$ are refinable parameters, $$\Delta\alpha$$ is the acute angle between the scattering vector and the normal to the crystallites (platy habit).

In 2D diffraction pattern the angle $$\Delta\alpha$$ is calculated as:

$$ \Delta\alpha = \alpha_{CSC} - \alpha_{det.},$$

where

$$ cos \alpha_{det.} = cos \theta_{det.} \cdot sin \phi_{det.} $$

$$ cos \alpha_{CSC} = \fract{\left( \vec{q}_{hkl} \cdot \vec{q}_{ax.} \right)}{q_{hkl} \cdot q_{ax.}} $$

in which $$2 \theta_{det.}, \phi_{det.}$$ are the angles describing the scattered beam, $$q_{ax.}$$ and $$q_{hkl}$$ is calculated as:

$$q_{ax.} = B \cdot \vec{hkl}_{ax}$$

$$q_{hkl} = B \cdot \vec{hkl}$$

$$\vec{hkl}_{ax}$$ describes the texture axis whic is along the vertical direction. $$\vec{hkl}$$ is the vector composed from the Miller indeces.

The parameters  $$G_1, G_1, h_{\text{ax.}}, k_{\text{ax.}}, l_{\text{ax.}}$$ are given in [Texture][texture] item.

![alpha](/fig_alpha.png){:class="img-responsive"}


[notes]: /cryspy/theory
[texture]: /cryspy/rcif-format/texture
