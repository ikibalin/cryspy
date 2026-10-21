"""Recipes: typical analysis algorithms of CrysPy for AI assistants.

A recipe is data: the steps of the algorithm, each with its goal, the
tool calls that do it, what counts as success and what to do otherwise.
An assistant (a language model) reads the recipe (tool get_recipe) and
calls the tools itself, step by step, judging each result. CrysPy gives it
the facts: chi-square before and after, warnings about unphysical values
(history.py) and the tool undo_refinement.

run_recipe follows a recipe without a model, taking the tool calls of the
steps as they are and undoing a refinement with warnings; it checks that
the recipes work and runs them from scripts:

    >>> import cryspy
    >>> obj = cryspy.load_file("main.rcif")
    >>> print(cryspy.run_recipe(obj, "refine_structure")["text"])
"""
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

S_UNDO = "Undo the refinement (undo_refinement) and go to the next step."


@dataclass
class RecipeStep:
    """One step of a recipe.

    actions: tool calls (name, arguments) doing the step, in their order;
    success: what the results of the actions should show;
    on_fail: what to do when they do not;
    condition: when the step applies ("" for always)."""

    name: str
    goal: str
    actions: Tuple[Tuple[str, dict], ...]
    success: str = "chi_sq is better or unchanged and there are no warnings."
    on_fail: str = S_UNDO
    condition: str = ""

    def to_text(self, number: int) -> str:
        s_action = "; ".join(f"{name:} {format_arguments(arguments):}"
                             for name, arguments in self.actions)
        ls_out = [f"{number:}. {self.name:}: {self.goal:}"]
        if self.condition:
            ls_out.append(f"   When: {self.condition:}")
        ls_out.append(f"   Actions: {s_action:}")
        ls_out.append(f"   Success: {self.success:}")
        ls_out.append(f"   Otherwise: {self.on_fail:}")
        return "\n".join(ls_out)


@dataclass
class Recipe:
    """A typical algorithm of analysis.

    when: requests of the user the recipe is for;
    report: what the final answer to the user contains."""

    name: str
    title: str
    when: str
    steps: List[RecipeStep] = field(default_factory=list)
    report: str = ""

    def to_text(self) -> str:
        ls_out = [f"Recipe {self.name:}: {self.title:}.",
                  "Do the steps in their order, one action at a time, and read the result of each "
                  "action before the next one. A step whose parameters are not found or cannot be "
                  "refined (refine_parameters gives an Error) is skipped."]
        ls_out.extend(step.to_text(number) for number, step in enumerate(self.steps, start=1))
        if self.report:
            ls_out.append(f"Final answer: {self.report:}")
        return "\n".join(ls_out)


def format_arguments(arguments: dict) -> str:
    return "{" + ", ".join(f'"{key:}": {value!r:}'.replace("'", '"')
                           for key, value in arguments.items()) + "}"


def refine_step(name: str, goal: str, pattern: str, **kwargs) -> RecipeStep:
    """A step that adds the parameters matching the pattern to the refined
    ones and refines."""
    return RecipeStep(name, goal, (("refine_parameters", {"pattern": pattern}), ("refine", {})),
                      **kwargs)


L_RECIPES: List[Recipe] = [
    Recipe(
        "refine_structure",
        "Step-by-step refinement of a crystal structure",
        "refine the structure, refine from scratch, improve the fit of the model to the data",
        [
            RecipeStep(
                "start", "Fix all parameters and calculate the start model.",
                (("fix_parameters", {}), ("calculate", {})),
                success="chi_sq of the start model is a number.",
                on_fail="Stop and explain the error: the model or the data must be corrected first."),
            refine_step("scale", "Scale factor.", "scale",
                        on_fail="Stop: without a scale factor the next steps make no sense."),
            refine_step("background", "Background and zero shift of the powder pattern.",
                        "pd_background, offset_ttheta", condition="powder experiments (pd, tof)"),
            refine_step("positions", "Free atomic positions.", "atomic positions"),
            refine_step("adp", "Isotropic thermal parameters.", "b_iso_or_equiv, u_iso_or_equiv"),
            refine_step("extinction", "Extinction.", "extinction",
                        condition="single crystal experiments"),
        ],
        report="chi_sq of the start model and at the end, the steps done, skipped and undone "
               "with the reason, and the warnings that remain."),
]


def take_recipe(name: str) -> Optional[Recipe]:
    for recipe in L_RECIPES:
        if recipe.name == name.strip():
            return recipe
    return None


def recipes_as_text() -> str:
    """One line for each recipe, for the prompt of a model."""
    return "\n".join(f"- {recipe.name:}: {recipe.title:} (for: {recipe.when:})"
                     for recipe in L_RECIPES)


def run_recipe(cryspy_object, name: str) -> dict:
    """Run a recipe without a language model: the actions of every step as
    they are written; a step is skipped when its parameters are not found,
    undone when its refinement gives warnings, and the recipe stops when
    a step without undo fails.

    Returns dict with keys 'text' (the report) and 'steps': a list of
    (name of the step, status: 'done', 'skipped', 'undone' or 'failed',
    text of the last action).
    """
    from cryspy.procedure_assistant.tools import call_assistant_tool
    recipe = take_recipe(name)
    if recipe is None:
        raise ValueError(f"Recipe '{name:}' is not known; known: "
                         + ", ".join(recipe.name for recipe in L_RECIPES))
    l_step, l_chi_sq = [], []
    for step in recipe.steps:
        status, text = "done", ""
        for action, arguments in step.actions:
            result = call_assistant_tool(cryspy_object, action, dict(arguments))
            text = str(result["text"])
            if result["flag_error"] or result["missing"]:
                status = "skipped" if action == "refine_parameters" else "failed"
                break
            if "chi_sq" in result["metrics"]:
                l_chi_sq.append(result["metrics"]["chi_sq"])
            if result["metrics"].get("warnings") and action == "refine":
                status = "undone" if step.on_fail == S_UNDO else "failed"
                if status == "undone":
                    result = call_assistant_tool(cryspy_object, "undo_refinement")
                    text = f"{text:}\n{result['text']:}"
                    if "chi_sq" in result["metrics"]:
                        l_chi_sq.append(result["metrics"]["chi_sq"])
                break
        l_step.append((step.name, status, text))
        if status == "failed":
            break
    ls_out = [f"Recipe {recipe.name:}:"]
    ls_out.extend(f"- {name:}: {status:}" for name, status, _ in l_step)
    if len(l_chi_sq) > 0:
        ls_out.append(f"chi_sq: {l_chi_sq[0]:.5g} at the start, {l_chi_sq[-1]:.5g} at the end.")
    return {"text": "\n".join(ls_out), "steps": l_step}
