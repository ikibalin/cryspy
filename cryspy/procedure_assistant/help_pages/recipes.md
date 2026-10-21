---
title: Recipes of the AI assistant (typical algorithms of analysis)
keywords: recipe, recipes, algorithm, step by step, write a recipe, new recipe, edit recipe, рецепт, рецепты, алгоритм, по шагам
categories:
---
A recipe is a text file with the steps of an analysis. The assistant reads
it (tool `get_recipe`) and calls the tools of each step itself, judging the
results; `cryspy.run_recipe(obj, name)` runs it without a model. The recipes
of CrysPy are in `cryspy/procedure_assistant/recipe_files/`; the recipes of
the user are in the folder given by `cryspy.set_user_recipe_folder(path)`
(CrysPy Editor: AI assistant → Recipes...). A recipe of the user replaces
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
- `Success:` and `Otherwise:` are optional; by default a step is good when
  chi_sq/n (and R_F2 for single crystals) is better or unchanged without
  warnings, otherwise the refinement is undone (`undo_refinement`).
