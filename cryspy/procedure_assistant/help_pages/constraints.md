---
title: Constraints between refined parameters
keywords: constraint, constraints, restraint, soft constraint, sigma, linear, equation, sum, occupancy, mark, label, x+y, punishment, связь, ограничение, сумма заселённостей
categories: constraint
---
A mark is a name written after the parentheses of a refined value:
`0.7()x`. Marks connect refined parameters:

- equal marks make values equal: `2.1()a` and `2.1()a`;
- a mark ending with `-` makes the value opposite: `-2.1()a-`;
- a linear constraint between different marks is written in the loop
  `constraint` of the global object:

```text
global_

loop_
_constraint_id
_constraint_expression
1 "x+y=1"
2 "2*a-b=0"

data_fe_co
...
_atom_site_occupancy
Fe1 ... 0.7()x
Co1 ... 0.3()y
```

An expression is a linear equation of the marks: numbers, marks, `+`, `-`,
`*` and `/` by numbers, parentheses. During the refinement the last
parameters of each equation are calculated from the others, so the
equation holds exactly.

A soft constraint (restraint) has a sigma: the equation `f = g` adds
`((f-g)/sigma)**2` to chi-square, so it holds within about sigma. Any formula
of the marks is allowed: `**`, `sqrt`, `exp`, `log`, `abs`, `sin`, `cos`,
`tan` (radians), `sind`, `cosd` (degrees). Hard and soft constraints share
one loop; `.` is written for the sigma of a hard one:

```text
loop_
_constraint_id
_constraint_expression
_constraint_sigma
1 "x+y=1"                 .
2 "sqrt(a**2+b**2)=2.5"   0.05
```

The assistant adds a constraint with the tool `add_constraint`, e.g.
expression `x+y=1` and parameters
`x=fe_co.atom_site[Fe1].occupancy, y=fe_co.atom_site[Co1].occupancy`
(plus `sigma` for a soft one). In CrysPy Editor: right-click the object
tree, "Add constraint..." or "Add soft constraint...".
