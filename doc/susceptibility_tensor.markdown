
# Local suSceptibility tensor

## Local coordinate system

The local susceptibility tensor is defined in unit less reciprocal unit cell (a*/|a*|, b*/|b*|, c*/|c*|). In fact, it is almost the same system as for atomic vibration parameters U_ij, where the reciprocal unit cell is used (a*, b*, c*).

## Representation of magnetization ellipsoid

Equations transformation from the local susceptibility tensor to the the thermal ellipsoid. The components of susceptibility tensor in Cartesian coordinate system have to be finded. After all negative eigenvalues of chi_ij are replaced by positive ones. In Cartesian coordinate system it can be done through the eigendecomposition:

$$ 
chi_ij(orto) = Q  D  Q^{-1} 
$$

where D is diagonal matrix whose diagonal elements are eigenvalues,  Q is composed by eigenvectors (https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix). Therefore
 
$$
U_ij(orto)= Q  |D|  Q^{-1}
$$

The logic of this mathematical trick is to keep the directions (given by Q) and absolute values along main axes of the susceptibility ellipsoid. 

The components of the susceptibility tensor can be negatives. But, in principle, the positive eigenvalues of the susceptibility tensor  should be expected.

The principle axis susceptibilities are the eigenvalues of the susceptibility tensor written in Cartesian coordinate system. To estimate it from the local susceptibility tensor defined in reciprocal unit cell the following equation should be used:

$$
\chi_{orto} = nM \cdot \chi_{loc} \cdot nM^T
$$
