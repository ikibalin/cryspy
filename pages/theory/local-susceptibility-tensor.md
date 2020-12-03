---
layout: page
permalink: /theory/local-susceptibility-tensor/
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

# Local susceptibility tensor

The local susceptibility tensor is defined in unit less reciprocal unit cell (a\*/\|a\*\|, b\*/\|b\*\|, c\*/\|c\*\|). In fact, it is almost the same system as for atomic vibration parameters U_ij, where the reciprocal unit cell is used (a\*, b\*, c\*).

# Representation of magnetization ellipsoid

Equations transformation from the local susceptibility tensor to the the thermal ellipsoid. The components of susceptibility tensor in Cartesian coordinate system have to be finded. After all negative eigenvalues of $$\chi_{ij}$$ are replaced by positive ones. In Cartesian coordinate system it can be done through the eigendecomposition:

$$ \chi_{ij}^{orto} = Q D Q^{-1} $$

where D is diagonal matrix whose diagonal elements are eigenvalues, Q is composed by [eigenvectors][wiki-decomposition]. Therefore

$$ U_{ij}^{orto}= Q \vert D \vert Q^{-1} $$

The logic of this mathematical trick is to keep the directions (given by $$Q$$) and absolute values along main axes of the susceptibility ellipsoid.

The components of the susceptibility tensor can be negatives. But, in principle, the positive eigenvalues of the susceptibility tensor should be expected.

The principle axis susceptibilities are the eigenvalues of the susceptibility tensor written in Cartesian coordinate system. To estimate it from the local susceptibility tensor defined in reciprocal unit cell the following equation should be used:

$$ \chi^{orto} = nM \cdot \chi_{loc} \cdot nM^T $$


[notes]: /cryspy/theory
[wiki-decomposition]: https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix