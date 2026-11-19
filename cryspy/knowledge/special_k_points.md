---
title: Special points of the Brillouin zone (propagation vectors)
keywords: propagation vector, k vector, special point, high-symmetry point, Brillouin zone, zone boundary, commensurate, TRIM, Lifshitz condition, star of k
---
# Special points of the Brillouin zone

Magnetic propagation vectors are most often special points of the Brillouin
zone of the paramagnetic group: commensurate vectors whose components are
0, 1/2, 1/3 or 1/4 in the conventional reciprocal basis. A vector that is
measured close to such a point (within the accuracy of the indices) is
usually that point: CrysPy takes a component within 0.001 of a fraction
with a denominator up to 12 as that fraction, and names the point when the
lattice is in the table below.

- k ≡ −k (2k is a reciprocal lattice vector: components 0 or 1/2 in a
  primitive lattice) are the time-reversal invariant points: the moments
  are real (collinear or coplanar with fixed phases), the magnetic cell is
  at most doubled along each axis, and every representation of the group of
  k gives a magnetic space group.
- Points on the zone boundary with k not ≡ −k (K, H, P, W) and lines
  (Λ, Δ: incommensurate components) need the Fourier components cos and sin.
- Satellites of different arms of the star of k (k transformed by the
  operations of the paramagnetic group) are magnetic domains.
- The reciprocal lattice of a centred cell contains only (h, k, l) allowed
  by the centring (C: h + k even; I: h + k + l even; F: all even or all
  odd), so (0, 1, 0) of a C-centred cell is a point of the zone boundary
  and not Γ.

Coordinates are in the conventional reciprocal basis (as the indices of the
reflections); the labels follow Bradley and Cracknell (The Mathematical
Theory of Symmetry in Solids, 1972) and the k-vector tables of the Bilbao
Crystallographic Server (KVEC). For the A and B settings of a centred
orthorhombic cell permute the axes of the C setting. Other lattices
(I- and F-centred orthorhombic, rhombohedral, monoclinic, triclinic) are not
listed: there the k ≡ −k points are those with 2k in the reciprocal lattice
of the centred cell.

| Lattice | Point | k | k ≡ −k | Notes |
|---|---|---|---|---|
| cubic P | X | 0 1/2 0 | yes | Three arms. |
| cubic P | M | 1/2 1/2 0 | yes | Three arms. |
| cubic P | R | 1/2 1/2 1/2 | yes | |
| cubic F | X | 0 1 0 | yes | Type I antiferromagnet of the fcc lattice; three arms. |
| cubic F | L | 1/2 1/2 1/2 | yes | Type II antiferromagnet (MnO, NiO); four arms. |
| cubic F | W | 1/2 1 0 | yes | Type III antiferromagnet. |
| cubic F | K | 3/4 3/4 0 | no | |
| cubic I | H | 0 0 1 | yes | Antiferromagnet of the bcc lattice (body centre opposite to the corner). |
| cubic I | N | 1/2 1/2 0 | yes | Six arms. |
| cubic I | P | 1/2 1/2 1/2 | no | P and −P are different points. |
| tetragonal P | Z | 0 0 1/2 | yes | |
| tetragonal P | X | 0 1/2 0 | yes | Two arms. |
| tetragonal P | M | 1/2 1/2 0 | yes | |
| tetragonal P | R | 0 1/2 1/2 | yes | Two arms. |
| tetragonal P | A | 1/2 1/2 1/2 | yes | |
| tetragonal I | Z | 0 0 1 | yes | Called M when c < a. The body centre opposite to the corner (ThCr2Si2-type antiferromagnets). |
| tetragonal I | X | 1/2 1/2 0 | yes | Two arms. |
| tetragonal I | N | 0 1/2 1/2 | yes | Four arms. |
| tetragonal I | P | 1/2 1/2 1/2 | no | |
| orthorhombic P | X | 1/2 0 0 | yes | |
| orthorhombic P | Y | 0 1/2 0 | yes | |
| orthorhombic P | Z | 0 0 1/2 | yes | |
| orthorhombic P | S | 1/2 1/2 0 | yes | |
| orthorhombic P | U | 1/2 0 1/2 | yes | |
| orthorhombic P | T | 0 1/2 1/2 | yes | |
| orthorhombic P | R | 1/2 1/2 1/2 | yes | |
| orthorhombic C | Y | 0 1 0 | yes | The same point as 1 0 0. |
| orthorhombic C | Z | 0 0 1/2 | yes | TbTe3 (Cmcm): k1 = (0, 0, 1/2). |
| orthorhombic C | T | 0 1 1/2 | yes | |
| orthorhombic C | S | 1/2 1/2 0 | yes | Two arms (1/2 1/2 0 and 1/2 −1/2 0); TbTe3: k2 = (1/2, 1/2, 0). |
| orthorhombic C | R | 1/2 1/2 1/2 | yes | Two arms. |
| hexagonal P | A | 0 0 1/2 | yes | |
| hexagonal P | M | 1/2 0 0 | yes | Three arms. |
| hexagonal P | L | 1/2 0 1/2 | yes | Three arms. |
| hexagonal P | K | 1/3 1/3 0 | no | 120° structures of triangular lattices. |
| hexagonal P | H | 1/3 1/3 1/2 | no | |
