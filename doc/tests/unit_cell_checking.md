# Checking unit cell:

$$
a = 3.4 \AA, b = 7.5 \AA, c = 12 \AA \\
\alpha = 87^\circ, \beta = 102^\circ, \gamma = 124^\circ
$$



### Results

Reciprocal unit cell parameters
$$
a^{*} = 0.36332942\AA^{-1}, b^{*} = 0.16133114\AA^{-1}, c^{*} = 0.08546101 \AA^{-1}\\
\alpha^{*} = 85.47855015^\circ, \beta^{*} = 77.53862078^\circ, \gamma^{*} = 55.85213672
^\circ
$$
Volume of unit cells:
$$
V = 247.36961584064142 \AA^3\\
V^* = 0.00404 \AA^{-3}
$$
$$
\phi^2 = 0.653506416775028
$$



Metric tensors:
$$
\mathbf{G_\text{direct}} = 
\left[\begin{matrix}
 11.56000 &-14.25942 & -8.48280\\
-14.25942 & 56.25000 &  4.71024\\
 -8.48280 &  4.71024 &144.00000
\end{matrix}\right]
$$
$$
\mathbf{G_\text{direct norm.}} = 
\left[\begin{matrix}
 1. & -0.5591929 & -0.20791169\\
-0.5591929 & 1. &   0.05233596\\
 -0.20791169 &   0.05233596 & 1.
\end{matrix}\right]
$$


$$
\mathbf{G_\text{recip.}} = 
\left[\begin{matrix}
 0.13201 & 0.03290 & 0.00670\\
 0.03290 & 0.02603 & 0.00109\\
 0.00670 & 0.00109 & 0.00730
\end{matrix}\right]
$$

$$
\mathbf{G_\text{recip.norm.}} = 
\left[\begin{matrix}
 1. & 0.56133054 & 0.21578149\\
 0.56133054 & 1. & 0.07883231\\
 0.21578149 & 0.07883231 & 1.
\end{matrix}\right]
$$



Matrices:
$$
\mathbf{B} = 
\left[\begin{matrix}
  0.36332942 & 0.0905601 &  0.0184409\\
 0  & 0.13351631 & -0.00436731\\
 0  & 0  & 0.08333333\\
\end{matrix}\right] \\

\mathbf{B^{-1}} = 
\left[\begin{matrix}
 2.75232 &-1.86682 &-0.70690\\
 0.00000 & 7.48972 & 0.39252\\
 0.00000 & 0.00000 &12.00000\\
\end{matrix}\right]
$$

$$
\mathbf{B_\text{norm}} = 
\left[\begin{matrix}
 1.00000  &0.56133 & 0.21578\\
 0.00000  &0.82759 &-0.05110\\
 0.00000  &0.00000 & 0.97510
\end{matrix}\right] \\

\mathbf{B_\text{norm}^{-1}} = 
\left[\begin{matrix}
 1.00000 &-0.67827& -0.25684\\
 0.00000 & 1.20833&  0.06333\\
 0.00000 & 0.00000&  1.02553
\end{matrix}\right]
$$

$$
\mathbf{M} = 
\left[\begin{matrix}
 2.75232&  0.00000&  0.00000\\
-1.86682&  7.48972&  0.00000\\
-0.70690&  0.39252& 12.00000
\end{matrix}\right] \\

\mathbf{M^{-1}} = 
\left[\begin{matrix}
 0.36333 & 0.00000&  0.00000\\
 0.09056 & 0.13352&  0.00000\\
 0.01844 &-0.00437&  0.08333\\
\end{matrix}\right]
$$

$$
\mathbf{M_\text{norm}} = 
\left[\begin{matrix}
 0.80951 & 0.00000 & 0.00000\\
-0.54906 & 0.99863 & 0.00000\\
-0.20791 & 0.05234 & 1.00000
\end{matrix}\right] \\

\mathbf{M_\text{norm}^{-1}} = 
\left[\begin{matrix}
 1.23532&  0.00000&  0.00000\\
 0.67920&  1.00137&  0.00000\\
 0.22129& -0.05241&  1.00000\\
\end{matrix}\right]
$$

$$
\mathbf{M^{-1}} \cdot \mathbf{B_\text{norm.}} =
\begin{bmatrix}
0.36332942 & 0.203947897 & 0.0783997613\\
0.0905600957&  0.16133114 & 0.0127181067\\
0.0184409045 &  0.00673708925 &  0.08546101
\end{bmatrix}
$$

$$
\mathbf{B_\text{norm.}^{-1}} \cdot \mathbf{M} =
\begin{bmatrix}
4.20008805 & -5.1808664 & -3.08204968 \\
-2.30048833 & 9.07487662 & 0.75990775 \\
-0.72494843 & 0.40254155 & 12.30638596\\
\end{bmatrix}
$$



# Results for transformations of quadratic form

Quadratic form in different coordinate systems:
$$
\mathbf{Q_\text{recip.norm}} =
\begin{bmatrix}
3.1 & 0.4 & -0.7 \\
0.4 & 2.9 & 5.2 \\
-0.7 & 5.2 & 5.7
\end{bmatrix}
$$

$$
\mathbf{Q_\text{XYZ}} =
\begin{bmatrix}
   3.10000&   -1.61930656&   -1.48873842\\
  -1.61930656&    5.00464377&   7.55124895\\
  -1.48873842&   7.55124895&    7.24204832
\end{bmatrix}
$$

$$
\mathbf{Q_\text{recip.}} =
\begin{bmatrix}
 0.40922562 &  0.02344654 &  -0.02173535 \\
 0.02344654 &   0.07548044 &   0.07169512 \\
 -0.02173535 &  0.07169512 &  0.04163043
\end{bmatrix}
$$

$$
\mathbf{Q_\text{direct}} =
\begin{bmatrix}
  86.90688& -152.48631& -279.76439\\
-152.48631&  326.25521&  712.79278\\
-279.76439&  712.79278& 1042.85496
\end{bmatrix}
$$

$$
\mathbf{Q_\text{direct norm.}} =
\begin{bmatrix}
 7.51789858 & -5.9798562 & -6.8569725\\
-5.9798562&  5.80008894 & 7.91992089\\
-6.8569725 &  7.91992089 &  7.24205
\end{bmatrix}
$$


For reflections:
$$
\text{index_hkl} = \begin{bmatrix}
0 & 0 & 2 \\
0 & 2 & 0 \\
2 & 0 & 0 \\
0 & 0 & 0 \\
1 & 2 & 3 \\
\end{bmatrix}
$$

$$
\text{q_ccs} = \begin{bmatrix}
0.03688181 & -0.00873463 & 0.16666667 \\
0.18112019 & 0.26703263 & 0. \\
0.72665883 & 0. & 0. \\
0. & 0. & 0. \\
0.59977232 & 0.25393068 & 0.25 \\
\end{bmatrix}
$$


$$
\text{eq_ccs} = \begin{bmatrix}
0.21578149 & -0.05110301 & 0.9751035 \\
0.56133054 & 0.8275917 & 0. \\
1. & 0. & 0. \\
0. & 0. & 0. \\
0.85971072 & 0.363983 & 0.35834878 \\
\end{bmatrix}
$$

$$
d^{-1} = \begin{bmatrix}
0.17092203 \\
0.32266228 \\
0.72665883 \\ 
0. \\
0.68715675
\end{bmatrix}
$$

$$
\frac{\sin \theta}{\lambda} = \begin{bmatrix}
0.08546101 \\ 
0.16133114 \\ 
0.36332942 \\ 
0. \\ 
0.34357837
\end{bmatrix}
$$

