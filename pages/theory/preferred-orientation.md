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

$$ P_h = G_2 + \left(1-G_2\right) \cdot \left(\left(G_1 \cos \Delta\alpha \right)^2 + \frac{\sin^2 \Delta\alpha} {G_1}\right)^{-3/2} $$

Where $$G_1$$ and $$G_2$$ are refinable parameters, $$\Delta\alpha$$ is the acute angle between the scattering vector and the normal to the crystallites (platy habit).

![alpha](/fig_alpha.png){:class="img-responsive"}

In 2D diffraction pattern the angle $$\Delta\alpha$$ is calculated as:

$$ \Delta\alpha = \alpha_{\text{CSC}} - \alpha_{\text{det.}},$$

where

$$ \cos \alpha_{\text{det.}} = \cos \theta_{\text{det.}} \cdot \sin \phi_{\text{det.}},$$

$$ \cos \alpha_{\text{CSC}} = \frac{\left( \vec{q}_{hkl} \cdot \vec{q}_{\text{ax.}} \right)}{q_{hkl} \cdot q_{\text{ax.}}},$$

in which $$2 \theta_{\text{det.}}, \phi_{\text{det.}}$$ are the angles describing the scattered beam, $$q_{\text{ax.}}$$ and $$q_{hkl}$$ is calculated as:

$$q_{\text{ax.}} = B \cdot \vec{hkl}_{\text{ax.}},$$

$$q_{hkl} = B \cdot \vec{hkl}.$$

$$\vec{hkl}_{\text{ax.}}$$ describes the texture axis which is along the vertical direction. $$\vec{hkl}$$ is the vector composed from the Miller indeces. The $$B$$ matrix is estimated from the unit cell parameters.

The parameters  $$G_1, G_1, h_{\text{ax.}}, k_{\text{ax.}}, l_{\text{ax.}}$$ are given in [Texture][texture] item.


[notes]: /cryspy/theory
[texture]: /cryspy/rcif-format/texture
