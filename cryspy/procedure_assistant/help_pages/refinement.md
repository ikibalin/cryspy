---
title: Refinement with RhoChi
keywords: refine, refinement, fit, chi square, chi2, rhochi, least squares, correlation, hessian, fix, free, уточнение, подгонка, хи-квадрат, корреляции, зафиксировать, освободить
categories: chi2, refine_ls, inversed_hessian, exclude, range
---
RhoChi refines the parameters marked for refinement against all experiments
of the object at once.

1. Mark parameters: write `()` after the value in rcif (`_cell_length_a 5.43()`),
   or in Python `item.<attribute>_refinement = True`, or use
   `refine_all_variables()` / `fix_variables()` of an item, loop or data block.
2. Run `cryspy.rhochi_rietveld_refinement(obj)`; it returns a dictionary with
   `chi_sq`, `parameter_name` and, if available, `correlation_matrix`. The
   object is updated in place, the inverse Hessian is stored in the
   `inversed_hessian` item. `rhochi_rietveld_refinement_with_parameters` takes
   another scipy minimizer (`optimization_method="basinhopping"`, ...).
3. `cryspy.rhochi_no_refinement(obj)` only calculates the model and chi-square.

Parameters constrained by symmetry (for example `fract_y = 1/4` on a mirror
plane, or `chi_22 = chi_11` in a cubic group) cannot be refined; their flag
`<attribute>_constraint` is True and they follow the free parameters.

Choice of the data used in chi-square is set per experiment:

- `chi2` category: `sum`, `diff`, `up`, `down` for polarized powder data,
  `asymmetry` for flip ratios;
- `range` category limits the angular (or time) range of powder data;
- `exclude` loop removes regions of a powder pattern;
- `_diffrn_refln_excluded True` removes a single crystal reflection.

After refinement check chi-square per point, the standard deviations and the
correlations above 0.7 in absolute value: strongly correlated parameters
should not be refined together.
