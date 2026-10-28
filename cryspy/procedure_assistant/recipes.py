"""Recipes: typical analysis algorithms of CrysPy for AI assistants.

A recipe is a text file: the steps of the algorithm, each with its goal,
the tool calls that do it, what counts as success and what to do
otherwise. An assistant (a language model) reads the recipe (tool
get_recipe) and calls the tools itself, step by step, judging each
result. CrysPy gives it the facts: chi-square before and after, warnings
about unphysical values (history.py) and the tool undo_refinement.

The recipes of CrysPy are in recipe_files/; the user writes his own in a
folder given by set_user_recipe_folder (CrysPy Editor gives its settings
folder). A recipe of the user replaces the one of CrysPy with the same
name. The format (lines starting with '#' are comments):

    Recipe refine_scale: Refinement of the scale factor
    Requests: refine the scale, fit the intensities
    For: single_crystal, unpolarized

    1. start: Fix all parameters and calculate the start model.
       Actions: fix_parameters; calculate
       Success: chi_sq/n of the start model is a number.
       Otherwise: Stop and explain the error.
    2. scale: Scale factor.
       Actions: refine_parameters {"pattern": "scale"}; refine

    Final answer: chi_sq/n at the start and at the end.

'For:' names the kinds of experiments the recipe is for (experiments.py:
single_crystal, powder, flip_ratio, nuclear, magnetic, ...), so that
suggest_recipes finds it. The names of the recipe and of the steps are
words of letters, digits and '_'; a name written in several words
('1. Define background: ...') is joined ('define_background'). 'Actions:'
are tool calls separated by ';', each a tool name with its arguments in
JSON; methods of the object are called with call_method (methods.py). A step may have 'When:' (the
condition), 'Success:' and 'Otherwise:' (by default: chi_sq/n, and R_F2 for
single crystals, is better or unchanged without warnings; otherwise undo).

An assistant writes a recipe from the words of the user with the tools
get_recipe_guide (form_recipe_guide: the format and the way of work),
ask_user (its questions to the user), check_recipe_text and save_recipe.

run_recipe follows a recipe without a model, taking the tool calls of the
steps as they are and undoing a refinement with warnings; it checks that
the recipes work and runs them from scripts:

    >>> import cryspy
    >>> obj = cryspy.load_file("main.rcif")
    >>> print(cryspy.run_recipe(obj, "refine_structure")["text"])
"""
import json
import logging
import os
import re
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

from cryspy.procedure_assistant.experiments import KINDS

logger = logging.getLogger(__name__)

DIR_RECIPE_FILES = os.path.join(os.path.dirname(__file__), "recipe_files")
RECIPE_EXTENSION = ".txt"
KIND_ANY = "any"  # a recipe for any object, never suggested for an experiment
S_UNDO = "Undo the refinement (undo_refinement) and go to the next step."
S_SUCCESS = ("chi_sq/n (and R_F2 for single crystals) is better or unchanged and there "
             "are no warnings.")
S_TEMPLATE = """# A recipe: lines starting with '#' are comments.
Recipe my_recipe: Short title of the recipe
Requests: words of the requests of the user the recipe is for
For: single_crystal, unpolarized

1. start: Fix all parameters and calculate the start model.
   Actions: fix_parameters; calculate
   Success: chi_sq/n of the start model is a number.
   Otherwise: Stop and explain the error.
2. scale: Scale factor.
   Actions: refine_parameters {"pattern": "scale"}; refine

Final answer: chi_sq/n at the start and at the end, the steps done, skipped and undone.
"""

# names of recipes and steps are words of letters, digits and '_'; a name
# written in several words ('Define background') is joined (form_name)
RE_HEADER = re.compile(r"^Recipe\s+([A-Za-z][\w\s\-]*?)\s*:\s*(.+)$")
RE_STEP = re.compile(r"^\s*(\d+)[.)]\s*([A-Za-z][\w\s\-]*?)\s*:\s*(.*)$")
RE_STEP_FIELD = re.compile(r"^\s+(When|Actions|Success|Otherwise)\s*:\s*(.*)$")
RE_RECIPE_FIELD = re.compile(r"^(Requests|For|Final answer)\s*:\s*(.*)$")
RE_TOOL_NAME = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")

D_USER_FOLDER = {"path": None}


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
    success: str = S_SUCCESS
    on_fail: str = S_UNDO
    condition: str = ""

    def to_text(self, number: int) -> str:
        ls_out = [f"{number:}. {self.name:}: {self.goal:}"]
        if self.condition:
            ls_out.append(f"   When: {self.condition:}")
        ls_out.append(f"   Actions: {format_actions(self.actions):}")
        ls_out.append(f"   Success: {self.success:}")
        ls_out.append(f"   Otherwise: {self.on_fail:}")
        return "\n".join(ls_out)


@dataclass
class Recipe:
    """A typical algorithm of analysis.

    when: requests of the user the recipe is for;
    kinds: kinds of the experiments the recipe is for (experiments.KINDS);
    report: what the final answer to the user contains;
    path: the file of the recipe."""

    name: str
    title: str
    when: str
    steps: List[RecipeStep] = field(default_factory=list)
    report: str = ""
    kinds: Tuple[str, ...] = ()
    path: str = ""

    def to_text(self) -> str:
        ls_out = [f"Recipe {self.name:}: {self.title:}."]
        if self.kinds:
            ls_out.append(f"For experiments: {', '.join(self.kinds):}.")
        ls_out.append(
            "Do the steps in their order, one action at a time, and read the result of each "
            "action before the next one. A step whose parameters are not found or cannot be "
            "refined (refine_parameters gives an Error) is skipped.")
        ls_out.extend(step.to_text(number) for number, step in enumerate(self.steps, start=1))
        if self.report:
            ls_out.append(f"Final answer: {self.report:}")
        return "\n".join(ls_out)

    def flag_user(self) -> bool:
        return self.path != "" and not self.path.startswith(DIR_RECIPE_FILES)


def format_arguments(arguments: dict) -> str:
    return json.dumps(arguments, ensure_ascii=False)


def format_actions(actions) -> str:
    return "; ".join(name if len(arguments) == 0 else f"{name:} {format_arguments(arguments):}"
                     for name, arguments in actions)


def parse_actions(text: str) -> Tuple[Tuple[str, dict], ...]:
    """'fix_parameters; refine_parameters {"pattern": "scale"}; refine' ->
    (("fix_parameters", {}), ("refine_parameters", {"pattern": "scale"}), ("refine", {}))."""
    decoder = json.JSONDecoder()
    l_action, position = [], 0
    while True:
        while position < len(text) and text[position] in " \t;":
            position += 1
        if position >= len(text):
            break
        match = RE_TOOL_NAME.match(text, position)
        if match is None:
            raise ValueError(f"a tool name is expected at '{text[position:]:}'")
        position = match.end()
        while position < len(text) and text[position] in " \t":
            position += 1
        arguments = {}
        if position < len(text) and text[position] == "{":
            try:
                arguments, position = decoder.raw_decode(text, position)
            except json.JSONDecodeError as error:
                raise ValueError(f"the arguments of {match.group(0):} are not JSON: {error.msg:}") from None
        elif position < len(text) and text[position] != ";":
            raise ValueError(f"';' or JSON arguments are expected after {match.group(0):}")
        l_action.append((match.group(0), arguments))
    return tuple(l_action)


def form_name(words: str) -> str:
    """Name of a recipe or a step from its words: 'Define background
    points' -> 'define_background_points'; a name already in one word
    stays as it is."""
    words = words.strip()
    if re.fullmatch(r"[A-Za-z]\w*", words):
        return words
    return re.sub(r"[\s\-]+", "_", words).strip("_").lower()


def parse_recipe(text: str, path: str = "") -> Recipe:
    """The recipe of the text (the format is in the module docstring);
    raises ValueError with the number of the wrong line."""
    recipe, d_step, s_field, d_field = None, None, "", None
    l_step = []
    for number, line in enumerate(text.splitlines(), start=1):
        if line.strip() == "" or line.strip().startswith("#"):
            continue
        try:
            if recipe is None:
                match = RE_HEADER.match(line.strip())
                if match is None:
                    raise ValueError("the recipe starts with 'Recipe <name>: <title>'")
                recipe = {"name": form_name(match.group(1)), "title": match.group(2).strip().rstrip("."),
                          "Requests": "", "For": "", "Final answer": ""}
                d_field, s_field = recipe, "title"
                continue
            match = RE_STEP.match(line)
            if match is not None:
                d_step = {"name": form_name(match.group(2)), "goal": match.group(3).strip(), "line": number}
                l_step.append(d_step)
                d_field, s_field = d_step, "goal"
                continue
            match = RE_RECIPE_FIELD.match(line)
            if match is not None:
                d_field, s_field = recipe, match.group(1)
                d_field[s_field] = match.group(2).strip()
                continue
            match = RE_STEP_FIELD.match(line)
            if match is not None:
                if d_step is None:
                    raise ValueError(f"'{match.group(1):}' outside of a step")
                d_field, s_field = d_step, match.group(1)
                d_field[s_field] = match.group(2).strip()
                continue
            if line[0] in " \t":  # the continuation of the field above
                d_field[s_field] = f"{d_field[s_field]:} {line.strip():}".strip()
                continue
            raise ValueError(
                f"unexpected line '{line.strip():}': a line is 'Requests:', 'For:', "
                "'Final answer:', a step '<number>. <step_name>: <goal>' (e.g. '1. background: "
                "Define the background points.') or an indented field of the step ('   Actions: "
                "...'); a long text goes on indented lines")
        except ValueError as error:
            raise ValueError(f"line {number:}: {error:}") from None
    if recipe is None:
        raise ValueError("the text is empty: the recipe starts with 'Recipe <name>: <title>'")
    if len(l_step) == 0:
        raise ValueError(f"recipe {recipe['name']:} has no steps")
    l_recipe_step = []
    for d_step in l_step:
        try:
            actions = parse_actions(d_step.get("Actions", ""))
        except ValueError as error:
            raise ValueError(f"line {d_step['line']:}: step {d_step['name']:}: {error:}") from None
        l_recipe_step.append(RecipeStep(
            d_step["name"], d_step["goal"], actions, success=d_step.get("Success") or S_SUCCESS,
            on_fail=d_step.get("Otherwise") or S_UNDO, condition=d_step.get("When", "")))
    kinds = tuple(word.lower() for word in re.split(r"[\s,;]+", recipe["For"]) if word != "")
    return Recipe(recipe["name"], recipe["title"], recipe["Requests"], l_recipe_step,
                  report=recipe["Final answer"], kinds=kinds, path=path)


def check_recipe(recipe: Recipe) -> List[str]:
    """Problems of a recipe that parses: unknown tools, arguments or kinds
    of experiments, repeated names of steps."""
    from cryspy.procedure_assistant.tools import take_assistant_tool
    l_problem = []
    l_kind = [kind for kind in recipe.kinds if kind not in KINDS + (KIND_ANY, )]
    if len(l_kind) > 0:
        l_problem.append(f"unknown kinds in 'For': {', '.join(l_kind):}; known: "
                         f"{', '.join(KINDS + (KIND_ANY, )):}")
    l_name = []
    for step in recipe.steps:
        if step.name in l_name:
            l_problem.append(f"step {step.name:} is repeated")
        l_name.append(step.name)
        for name, arguments in step.actions:
            tool = take_assistant_tool(name)
            if tool is None:
                l_problem.append(f"step {step.name:}: unknown tool {name:}{explain_method_name(name):}")
                continue
            if name == "call_method" and str(arguments.get("method", "")) not in take_method_names():
                l_problem.append(f"step {step.name:}: call_method: no CrysPy class has the method "
                                 f"'{arguments.get('method', ''):}' (list_methods gives them)")
            l_argument = [key for key in arguments if key not in tool.parameters]
            if len(l_argument) > 0 and not tool.flag_extra_arguments:
                l_problem.append(f"step {step.name:}: {name:} has no arguments "
                                 f"{', '.join(l_argument):}")
    return l_problem


def take_method_names() -> dict:
    """{name of a method: names of the classes with it} of the methods an
    assistant may call (methods.py)."""
    from cryspy.procedure_assistant.methods import take_classes_with_methods, take_methods
    d_name = {}
    for cls in take_classes_with_methods():
        for name in take_methods(cls):
            d_name.setdefault(name, []).append(cls.__name__)
    return d_name


def explain_method_name(name: str) -> str:
    """'; it is a method of Pd: call_method {...}' for the name of a method
    written as a tool, else ''."""
    l_class = take_method_names().get(name)
    if l_class is None:
        return ""
    return (f"; it is a method of {', '.join(l_class[:5]):}: write call_method "
            f"{{\"method\": \"{name:}\"}} (with its arguments; \"path\" to call it for one block only)")


def set_user_recipe_folder(path: Optional[str]):
    """Folder of the recipes of the user (None: only the recipes of CrysPy)."""
    D_USER_FOLDER["path"] = path


def take_user_recipe_folder() -> Optional[str]:
    return D_USER_FOLDER["path"]


def take_recipe_paths() -> List[str]:
    """Files of the recipes: those of CrysPy, then those of the user."""
    l_path = []
    for folder in (DIR_RECIPE_FILES, D_USER_FOLDER["path"]):
        if folder is None or not os.path.isdir(folder):
            continue
        l_path.extend(os.path.join(folder, name) for name in sorted(os.listdir(folder))
                      if name.endswith(RECIPE_EXTENSION))
    return l_path


def read_recipe_file(path: str) -> Recipe:
    with open(path, "r", encoding="utf-8") as fid:
        return parse_recipe(fid.read(), path=path)


def list_recipe_files() -> List[dict]:
    """Every recipe file with its 'path', 'name', 'title', 'flag_user'
    (a recipe of the user), 'flag_replaced' (a recipe of CrysPy replaced
    by one of the user) and 'error' ("" when the recipe is correct)."""
    l_out = []
    for path in take_recipe_paths():
        d_out = {"path": path, "name": os.path.splitext(os.path.basename(path))[0], "title": "",
                 "flag_user": not path.startswith(DIR_RECIPE_FILES), "flag_replaced": False,
                 "error": ""}
        try:
            recipe = read_recipe_file(path)
            d_out.update(name=recipe.name, title=recipe.title)
            d_out["error"] = "; ".join(check_recipe(recipe))
        except (OSError, UnicodeDecodeError, ValueError) as error:
            d_out["error"] = str(error)
        l_out.append(d_out)
    l_user = [d_out["name"] for d_out in l_out if d_out["flag_user"] and d_out["error"] == ""]
    for d_out in l_out:
        d_out["flag_replaced"] = not d_out["flag_user"] and d_out["name"] in l_user
    return l_out


def take_recipes() -> List[Recipe]:
    """The correct recipes; a recipe of the user replaces the one of CrysPy
    with the same name. Wrong files are left out (list_recipe_files gives
    their errors)."""
    d_recipe = {}
    for path in take_recipe_paths():
        try:
            recipe = read_recipe_file(path)
        except (OSError, UnicodeDecodeError, ValueError) as error:
            logger.warning("Recipe %s is not read: %s", path, error)
            continue
        l_problem = check_recipe(recipe)
        if len(l_problem) > 0:
            logger.warning("Recipe %s is not used: %s", path, "; ".join(l_problem))
            continue
        d_recipe[recipe.name] = recipe
    return list(d_recipe.values())


def save_user_recipe(text: str) -> str:
    """Check the recipe and write it to the folder of the user as
    <name>.txt; returns the path. Raises ValueError for a wrong recipe or
    when the folder of the user is not set."""
    folder = D_USER_FOLDER["path"]
    if folder is None:
        raise ValueError("the folder of the recipes of the user is not set")
    recipe = parse_recipe(text)
    l_problem = check_recipe(recipe)
    if len(l_problem) > 0:
        raise ValueError("; ".join(l_problem))
    os.makedirs(folder, exist_ok=True)
    path = os.path.join(folder, f"{recipe.name:}{RECIPE_EXTENSION:}")
    with open(path, "w", encoding="utf-8") as fid:
        fid.write(text if text.endswith("\n") else f"{text:}\n")
    return path


def delete_user_recipe(name: str) -> bool:
    """Delete the recipe of the user; the recipe of CrysPy with the same
    name, if any, is used again. False when there is no such file."""
    folder = D_USER_FOLDER["path"]
    if folder is None:
        return False
    path = os.path.join(folder, f"{name:}{RECIPE_EXTENSION:}")
    if not os.path.isfile(path):
        return False
    os.remove(path)
    return True


def find_recipe_advice(recipe: Recipe) -> List[str]:
    """Remarks on a correct recipe that are not errors: what a recipe
    written from the words of the user often misses."""
    l_advice = []
    for ind, step in enumerate(recipe.steps):
        l_name = [name for name, _ in step.actions]
        if "refine_parameters" in l_name and "refine" not in l_name[l_name.index("refine_parameters"):]:
            l_later = [name for later in recipe.steps[ind + 1:] for name, _ in later.actions]
            if "refine" not in l_later:
                l_advice.append(f"step {step.name:} marks parameters for refinement, but no step "
                                "refines them (refine)")
        if len(step.actions) == 0 and step.success == S_SUCCESS:
            l_advice.append(f"step {step.name:} has no actions: the assistant does it by its goal; "
                            "write in 'Success:' what it should give")
    l_cryspy = [os.path.splitext(name)[0] for name in os.listdir(DIR_RECIPE_FILES)
                if name.endswith(RECIPE_EXTENSION)]
    if recipe.name in l_cryspy:
        l_advice.append(f"the name {recipe.name:} is the one of a recipe of CrysPy: yours replaces it")
    return l_advice


def check_recipe_text(text: str) -> str:
    """Check the text of a recipe: 'Error: ...' (where and what is wrong)
    or the short contents of the correct recipe with remarks."""
    try:
        recipe = parse_recipe(text)
    except ValueError as error:
        return f"Error: {error:}. The format: get_recipe_guide."
    l_problem = check_recipe(recipe)
    if len(l_problem) > 0:
        return f"Error in recipe {recipe.name:}: {'; '.join(l_problem):}."
    ls_out = [f"Recipe {recipe.name:} is correct: {len(recipe.steps):} steps "
              f"({', '.join(step.name for step in recipe.steps):})."]
    if not recipe.kinds:
        ls_out.append("It has no 'For:', so it is not suggested for any experiment.")
    l_advice = find_recipe_advice(recipe)
    if l_advice:
        ls_out.append("Remarks: " + "; ".join(l_advice) + ".")
    return "\n".join(ls_out)


S_RECIPE_FORMAT = """The format of a recipe (lines starting with '#' are comments):
Recipe <name>: <title>              name: letters, digits and '_'
Requests: <words of the requests of the user the recipe is for>
For: <kinds of experiments, all must apply>
<number>. <step_name>: <goal of the step>   step_name: one word, e.g. background
   When: <condition of the step, optional>
   Actions: <tool> {<arguments in JSON>}; <tool>; ...   (tools of the assistant)
   Success: <what the results should show, optional>
   Otherwise: <what to do when they do not, optional>
Final answer: <what the answer to the user contains>
By default a step is good when chi_sq/n (and R_F2 for single crystals) is better or
unchanged without warnings; otherwise the refinement is undone (undo_refinement).
Parameters are chosen by refine_parameters {"pattern": ...}: words ('scale',
'atomic positions', 'adp', 'extinction', 'unit cell', 'polarization') or paths
('atom_site.fract_x', '*.pd_background[*].intensity'); then refine.
Methods of the data blocks, loops and items (list_methods), e.g. the method
define_background_points of powder experiments, are tools too:
call_method {"method": "define_background_points", "step_ttheta": 2}
(for every block that has it; "path": "<block>" for one)."""

S_RECIPE_WRITING = """How to write a recipe from the words of the user:
1. Think over the algorithm: its steps in their order, the parameters of every step,
   the conditions ('When:'), what counts as success and what to do otherwise.
2. When something is not clear (which parameters, a condition, the kinds of
   experiments, what to do on a failure), ask the user with ask_user, all the
   questions of the moment together; never ask what a tool can find.
3. Write every step with the tools of the assistant only; an operation the user
   names that is not a tool is often a method of the object: list_methods
   {"pattern": "<a word of it>"} finds it, call_method calls it. get_recipe <name>
   shows a recipe of the CookBook as an example.
4. Check the text with check_recipe_text and correct it until it is correct.
5. Show the recipe to the user and ask with ask_user whether to save it or what to
   change; save it with save_recipe only after the user agrees.
6. Final answer: the name of the recipe, the file where it is saved and its text."""


def form_recipe_guide() -> str:
    """Everything a language model needs to write a recipe: the format, the
    kinds of experiments, the recipes of the CookBook, a template and the
    way of work."""
    # the most needed first: a program may give a model only the start of a long result
    return "\n\n".join([
        S_RECIPE_FORMAT,
        f"A template:\n{S_TEMPLATE.strip():}",
        S_RECIPE_WRITING,
        f"Kinds of experiments for 'For:': {', '.join(KINDS + (KIND_ANY, )):} "
        f"('{KIND_ANY:}': never suggested for an experiment).",
        f"Recipes of the CookBook (examples, get_recipe gives the steps):\n{recipes_as_text():}"])


def take_recipe(name: str) -> Optional[Recipe]:
    for recipe in take_recipes():
        if recipe.name == name.strip():
            return recipe
    return None


def recipes_as_text() -> str:
    """One line for each recipe, for the prompt of a model."""
    ls_out = []
    for recipe in take_recipes():
        s_kind = f"; experiments: {', '.join(recipe.kinds):}" if recipe.kinds else ""
        ls_out.append(f"- {recipe.name:}: {recipe.title:} (for: {recipe.when:}{s_kind:})")
    return "\n".join(ls_out)


def suggest_recipes(l_kind: List[str]) -> List[Recipe]:
    """Recipes for an experiment of these kinds: every kind of 'For' of the
    recipe is one of them; the most specific recipes first."""
    l_recipe = [recipe for recipe in take_recipes()
                if len(recipe.kinds) > 0 and KIND_ANY not in recipe.kinds
                and set(recipe.kinds) <= set(l_kind)]
    return sorted(l_recipe, key=lambda recipe: -len(recipe.kinds))


def run_recipe(cryspy_object, name: str, *, callback=None) -> dict:
    """Run a recipe without a language model: the actions of every step as
    they are written; a step is skipped when its parameters are not found
    or it has no actions, undone when its refinement gives warnings, and
    the recipe stops when a step without undo fails.

    callback (keyword-only, optional) is called after every step as
    callback(name of the step, status, text of the last action), e.g. to
    show the course of the work.

    Returns dict with keys 'text' (the report) and 'steps': a list of
    (name of the step, status: 'done', 'skipped', 'undone' or 'failed',
    text of the last action).
    """
    from cryspy.procedure_assistant.tools import call_assistant_tool
    recipe = take_recipe(name)
    if recipe is None:
        raise ValueError(f"Recipe '{name:}' is not known; known: "
                         + ", ".join(recipe.name for recipe in take_recipes()))
    l_step, l_metrics = [], []
    for step in recipe.steps:
        status, text = ("done", "") if len(step.actions) > 0 else ("skipped", "No actions.")
        for action, arguments in step.actions:
            result = call_assistant_tool(cryspy_object, action, dict(arguments))
            text = str(result["text"])
            if result["flag_error"] or result["missing"]:
                status = "skipped" if action == "refine_parameters" else "failed"
                break
            if "chi_sq_per_n" in result["metrics"]:
                l_metrics.append(result["metrics"])
            if result["metrics"].get("warnings") and action == "refine":
                status = "undone" if "undo_refinement" in step.on_fail else "failed"
                if status == "undone":
                    result = call_assistant_tool(cryspy_object, "undo_refinement")
                    text = f"{text:}\n{result['text']:}"
                    if "chi_sq_per_n" in result["metrics"]:
                        l_metrics.append(result["metrics"])
                break
        l_step.append((step.name, status, text))
        if callback is not None:
            callback(step.name, status, text)
        if status == "failed":
            break
    ls_out = [f"Recipe {recipe.name:}:"]
    ls_out.extend(f"- {name:}: {status:}" for name, status, _ in l_step)
    if len(l_metrics) > 0:
        d_start, d_end = l_metrics[0], l_metrics[-1]
        ls_out.append(f"chi_sq/n: {d_start['chi_sq_per_n']:.5g} at the start, "
                      f"{d_end['chi_sq_per_n']:.5g} at the end.")
        for name, r_f2 in d_end.get("r_f2", {}).items():
            if name in d_start.get("r_f2", {}):
                ls_out.append(f"R_F2 of {name:}: {d_start['r_f2'][name]:.4f} at the start, "
                              f"{r_f2:.4f} at the end.")
    return {"text": "\n".join(ls_out), "steps": l_step}
