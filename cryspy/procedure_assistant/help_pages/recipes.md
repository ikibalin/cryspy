---
title: Recipes of the AI assistant (typical algorithms of analysis)
keywords: recipe, recipes, cookbook, algorithm, step by step, write a recipe, new recipe, edit recipe, compose a recipe, рецепт, рецепты, алгоритм, по шагам, написать рецепт, составить рецепт
categories:
---
A recipe is a text file with the steps of an analysis. The assistant reads
it (tool `get_recipe`) and calls the tools of each step itself, judging the
results; `cryspy.run_recipe(obj, name)` runs it without a model. The recipes
of CrysPy are in `cryspy/procedure_assistant/recipe_files/`; the recipes of
the user are in the folder given by `cryspy.set_user_recipe_folder(path)`
(CrysPy Editor: AI assistant → CookBook...). A recipe of the user replaces
the one of CrysPy with the same name.

```
# lines starting with '#' are comments
Recipe refine_scale: Refinement of the scale factor
Requests: refine the scale, fit the intensities
For: single_crystal, unpolarized

1. start: Fix all parameters and calculate the start model.
   Actions: fix_parameters; calculate
   Success: chi_sq/n of the start model is a number.
   Otherwise: Stop and explain the error.
2. scale: Scale factor.
   When: single crystal experiments
   Actions: refine_parameters {"pattern": "scale"}; refine

Final answer: chi_sq/n at the start and at the end.
```

- `For:` the kinds of experiments the recipe is for, all must apply:
  `single_crystal`, `powder`, `powder_1d`, `powder_2d`, `tof`, `polarized`,
  `unpolarized`, `flip_ratio`, `nuclear` (the crystal has a nuclear space
  group), `magnetic` (a magnetic space group); `any` for a recipe that is
  never suggested for an experiment. The tool `suggest_recipes` finds the
  recipes for every experiment of the object (recipe `choose_recipes`).
- `Actions:` tool calls separated by `;`: the name of the tool and its
  arguments in JSON.
- `When:` the condition of a step. When it starts with kinds of
  experiments (`powder_1d`, `single_crystal, unpolarized` = both,
  `tof or powder_1d` = either; a remark may follow after `;`),
  `run_recipe` skips the step when no experiment has these kinds; other
  conditions are judged by the assistant.
- `Success:` and `Otherwise:` are optional; by default a step is good when
  chi_sq/n (and R_F2 for single crystals) is better or unchanged without
  warnings, otherwise the refinement is undone (`undo_refinement`).

The recipes of CrysPy for the request "fit my data": `refine_structure`
chooses by the experiments `refine_powder` (background, scale, zero shift,
unit cell, peak widths, atoms), `refine_tof` or `refine_single_crystal`
(scale, extinction, λ/2, atoms); `refine_atoms` (positions, ADP,
occupancies) and `check_result` are their bricks. `refine_flip_ratio_nuclear`,
`refine_polarized_powder` and `refine_magnetic_structure` follow a good
nuclear refinement. Two tools help a small model: `diagnose_fit` reads the
residuals and names what is wrong (shifted or too narrow peaks, background,
scale, a trend with sin θ/λ, extinction, unexplained peaks), and
`scan_parameter` calculates chi_sq/n on a grid of one parameter and keeps
the best value (a start value far from the minimum). `refine` names the
parameters that stayed exactly zero (a zero extinction or moment has no
gradient) and `set_parameter` with `if_zero` gives them start values.

When the user asks to refine parameters one by one and to keep each one only
if chi_sq gets better (and no high correlation appears), `refine_one_by_one`
does the whole loop (forward selection): each candidate matching the pattern
is refined together with the parameters kept before it; it stays refined when
chi_sq/n decreases by more than `min_improvement` (0.1 % by default) and no
new pair of refined parameters reaches |r| >= `correlation_limit` (0.9),
otherwise the values return and it is fixed. With `rejected_value` (0 for
magnetic moments or their Fourier components) every candidate is compared
with the model without it; `start_value` moves a zero candidate away from 0
when its gradient vanishes there. The result gives the decision and the
reason for every parameter; `undo_refinement` returns the whole selection.

The assistant writes a recipe from an algorithm described in words: it
reads the format (tool `get_recipe_guide`), asks the user what is not clear
(tool `ask_user`, the questions come as `missing` like the other questions
of the tools), checks its text (tool `check_recipe_text`), shows it and,
when the user agrees, saves it to the folder of the user (tool
`save_recipe`, which asks before replacing a recipe of the user). A step of
a recipe may ask the user too: `Actions: ask_user {"question": "..."}`.
