"""MkDocs hooks: pages generated from the package at build time.

- guides/: the help pages of cryspy/procedure_assistant/help_pages
  (the same texts that cryspy.get_help shows);
- examples/: one page per folder of examples/;
- rcif-categories.md: all rcif categories with their tags, from the classes;
- api/: API reference, one page per subpackage, rendered by mkdocstrings
  from the numpy-style docstrings.
"""
import os
import re

from mkdocs.structure.files import File

DIR_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DIR_PACKAGE = os.path.join(DIR_ROOT, "cryspy")
DIR_HELP = os.path.join(DIR_PACKAGE, "procedure_assistant", "help_pages")
DIR_EXAMPLES = os.path.join(DIR_ROOT, "examples")
URL_TREE = "https://github.com/ikibalin/cryspy/tree/develop/examples"

EXAMPLES = {
    "rhochi_unpolarized_powder_1d_PbSO4":
        "Unpolarized neutron powder diffraction, PbSO4: Rietveld refinement.",
    "rhochi_unpolarized_powder_1d_xrays_PbSO4":
        "X-ray powder diffraction, PbSO4: Rietveld refinement.",
    "rhochi_unpolarized_powder_1d_mcif_URu0.96Rh0.04Si2":
        "Unpolarized neutron powder diffraction with a magnetic space group (mcif), URu0.96Rh0.04Si2.",
    "rhochi_unpolarized_tof_powder_CeCuAl":
        "Time-of-flight neutron powder diffraction, CeCuAl.",
    "rhochi_polarized_powder_1d_Ho2Ti2O7_5K_1T":
        "Polarized neutron powder diffraction (1D), Ho2Ti2O7 at 5 K and 1 T: local susceptibility.",
    "rhochi_polarized_powder_2d_DyAl_5K_5T":
        "Polarized neutron powder diffraction (2D), DyAl at 5 K and 5 T.",
    "rhochi_polozed_powder_2d_example_HB-3A":
        "Polarized neutron powder diffraction (2D), example of the HB-3A instrument.",
    "rhochi_single_susceptibility_Ho2Ti2O7":
        "Polarized neutron single-crystal diffraction (flipping ratios), Ho2Ti2O7: local susceptibility.",
    "rhochi_single_susceptibility_Ho2Ti2O7_simulation_und":
        "Simulation of flipping ratios of Ho2Ti2O7 with a script.",
    "mempy_spin_density_YTiO3":
        "Spin density of YTiO3 by the maximum entropy method.",
    "mempy_magnetization_density_Yb2Ti2O7_2K_1T":
        "Magnetization density of Yb2Ti2O7 at 2 K and 1 T by the maximum entropy method.",
}

API_PAGES = [
    ("procedure_rhochi", "Refinement (RhoChi)"),
    ("procedure_mempy", "MEM reconstruction (MEMPy)"),
    ("procedure_simulation", "Simulation"),
    ("procedure_assistant", "Help and assistant tools"),
    ("E_data_classes", "Data blocks"),
    ("C_item_loop_classes", "Categories (items and loops)"),
    ("B_parent_classes", "Parent classes"),
    ("D_functions_item_loop", "Functions on items and loops"),
    ("F_functions_data", "Functions on data blocks"),
    ("H_functions_global", "Reading and writing rcif"),
    ("A_functions_base", "Numerical functions"),
    ("exceptions", "Exceptions"),
]
# pages of modules that exist in this checkout
API_PAGES = [(package, title) for package, title in API_PAGES
             if os.path.exists(os.path.join(DIR_PACKAGE, package))
             or os.path.exists(os.path.join(DIR_PACKAGE, package + ".py"))]

EXAMPLE_TITLES = {
    "rhochi_unpolarized_powder_1d_PbSO4": "PbSO4, neutron powder",
    "rhochi_unpolarized_powder_1d_xrays_PbSO4": "PbSO4, X-ray powder",
    "rhochi_unpolarized_powder_1d_mcif_URu0.96Rh0.04Si2": "URu0.96Rh0.04Si2, mcif powder",
    "rhochi_unpolarized_tof_powder_CeCuAl": "CeCuAl, TOF powder",
    "rhochi_polarized_powder_1d_Ho2Ti2O7_5K_1T": "Ho2Ti2O7, polarized powder 1D",
    "rhochi_polarized_powder_2d_DyAl_5K_5T": "DyAl, polarized powder 2D",
    "rhochi_polozed_powder_2d_example_HB-3A": "HB-3A, polarized powder 2D",
    "rhochi_single_susceptibility_Ho2Ti2O7": "Ho2Ti2O7, flipping ratios",
    "rhochi_single_susceptibility_Ho2Ti2O7_simulation_und": "Ho2Ti2O7, simulation script",
    "mempy_spin_density_YTiO3": "YTiO3, MEM spin density",
    "mempy_magnetization_density_Yb2Ti2O7_2K_1T": "Yb2Ti2O7, MEM magnetization density",
}

N_LINE_MAX = 120  # lines of an example rcif shown on its page


def _help_pages():
    l_page = []
    for f_name in sorted(os.listdir(DIR_HELP)):
        if not f_name.endswith(".md"):
            continue
        with open(os.path.join(DIR_HELP, f_name), encoding="utf-8") as fid:
            text = fid.read()
        title = f_name[:-3].replace("_", " ").capitalize()
        match = re.match(r"---\n(.*?)\n---\n", text, flags=re.S)
        if match is not None:
            for line in match.group(1).splitlines():
                if line.startswith("title:"):
                    title = line[len("title:"):].strip()
            text = text[match.end():]
        l_page.append((f"guides/{f_name:}", title, f"# {title:}\n\n{text:}"))
    return l_page


def _example_dirs():
    return [name for name in sorted(os.listdir(DIR_EXAMPLES))
            if name in EXAMPLES and os.path.isdir(os.path.join(DIR_EXAMPLES, name))]


def _example_page(name):
    dir_example = os.path.join(DIR_EXAMPLES, name)
    f_names = sorted(os.listdir(dir_example))
    ls_out = [f"# {EXAMPLE_TITLES[name]:}", "", EXAMPLES[name], "",
              f"Folder: `examples/{name:}`", ""]
    l_image = [f for f in f_names if f.lower().endswith(".png")]
    for f_image in l_image:
        ls_out += [f"![{name:}]({name:}/{f_image:})", ""]
    l_rcif = [f for f in f_names if f.endswith(".rcif")]
    l_script = [f for f in f_names if f.endswith(".py")]
    ls_out += ["Files: " + ", ".join(
        f"[{f:}]({URL_TREE:}/{name:}/{f:})" for f in f_names), ""]
    if l_rcif and not l_script:
        if name.startswith("rhochi"):
            procedure = "rhochi_rietveld_refinement"
        elif "spin_density" in name:
            procedure = "mempy_spin_density_reconstruction"
        else:
            procedure = "mempy_magnetization_density_reconstruction"
        ls_out += ["## Run it", "", "```python", "import cryspy", "",
                   f'rcif = cryspy.load_file("{l_rcif[0]:}")',
                   f"result = cryspy.{procedure:}(rcif)",
                   'print(result["chi_sq"])', "```", ""]
    for f_script in l_script:
        with open(os.path.join(dir_example, f_script), encoding="utf-8") as fid:
            ls_out += [f"## {f_script:}", "", "```python", fid.read().rstrip(), "```", ""]
    for f_rcif in l_rcif:
        with open(os.path.join(dir_example, f_rcif), encoding="utf-8", errors="replace") as fid:
            lines = fid.read().splitlines()
        text = "\n".join(lines[:N_LINE_MAX])
        if len(lines) > N_LINE_MAX:
            text += f"\n# ... {len(lines) - N_LINE_MAX:} more lines in the file"
        ls_out += [f'??? example "{f_rcif:}"', "", "    ```text"]
        ls_out += ["    " + line for line in text.splitlines()]
        ls_out += ["    ```", ""]
    return l_image, "\n".join(ls_out)


def _examples_index():
    ls_out = ["# Examples", "",
              "Each example is a folder of [examples/](" + URL_TREE + ") with an rcif file "
              "that holds the model, the data and the refinement flags. Download "
              "it, run the refinement, change the model.", "",
              "| Example | What it shows |", "| --- | --- |"]
    for name in _example_dirs():
        ls_out.append(f"| [{EXAMPLE_TITLES[name]:}]({name:}.md) | {EXAMPLES[name]:} |")
    return "\n".join(ls_out) + "\n"


def _type_name(val_type):
    return getattr(val_type, "__name__", str(val_type))


def _categories_page():
    import cryspy  # noqa: F401  (registers all classes)
    from cryspy.H_functions_global.function_1_cryspy_objects import (
        L_DATA_CLASS, L_ITEM_CLASS, L_LOOP_CLASS)

    ls_out = [
        "# rcif categories", "",
        "Generated from the CrysPy classes. A category is written either as "
        "single items (`_cell_length_a 5.0`) or, when it has a loop class, as "
        "a `loop_` table. **Refinable** tags take refinement flags "
        "(`5.0()`, `5.0(1)`).", "",
        "## Data blocks", "",
        "| Data block | Class | Mandatory categories | Optional categories |",
        "| --- | --- | --- | --- |"]
    for cls in L_DATA_CLASS:
        mand = ", ".join(f"`{c.ITEM_CLASS.PREFIX if hasattr(c, 'ITEM_CLASS') else c.PREFIX}`"
                         for c in cls.CLASSES_MANDATORY)
        opt = ", ".join(f"`{c.ITEM_CLASS.PREFIX if hasattr(c, 'ITEM_CLASS') else c.PREFIX}`"
                        for c in cls.CLASSES_OPTIONAL)
        ls_out.append(f"| `{cls.PREFIX}` | `{cls.__name__}` | {mand} | {opt} |")

    loop_by_item = {loop_cls.ITEM_CLASS: loop_cls for loop_cls in L_LOOP_CLASS}
    item_classes = list(L_ITEM_CLASS)
    for item_cls in loop_by_item:
        if item_cls not in item_classes:
            item_classes.append(item_cls)

    ls_out += ["", "## Categories", ""]
    for item_cls in sorted(item_classes, key=lambda c: c.PREFIX):
        loop_cls = loop_by_item.get(item_cls)
        kind = f"item `{item_cls.__name__}`"
        if loop_cls is not None:
            kind += f", loop `{loop_cls.__name__}`"
        ls_out += [f"### {item_cls.PREFIX}", "", kind, ""]
        doc = (item_cls.__doc__ or "").strip()
        if doc:
            ls_out += [re.sub(r"\n\s*", "\n", doc), ""]
        ls_out += ["| Tag | Type | Mandatory | Refinable | Allowed values |",
                   "| --- | --- | --- | --- | --- |"]
        mandatory = set(item_cls.ATTR_MANDATORY_NAMES)
        constraints = getattr(item_cls, "D_CONSTRAINTS", {})
        for name, name_cif, val_type in zip(
                item_cls.ATTR_NAMES, item_cls.ATTR_CIF, item_cls.ATTR_TYPES):
            allowed = ", ".join(f"`{v}`" for v in constraints.get(name, []))
            ls_out.append(
                f"| `_{item_cls.PREFIX}_{name_cif}` | {_type_name(val_type)} | "
                f"{'yes' if name in mandatory else ''} | "
                f"{'yes' if name in item_cls.ATTR_REF else ''} | {allowed} |")
        ls_out.append("")
    return "\n".join(ls_out)


def _api_page(package, title):
    if package == "exceptions":
        return f"# {title}\n\n::: cryspy.exceptions\n    options:\n      show_root_heading: false\n"
    dir_package = os.path.join(DIR_PACKAGE, package)
    ls_out = [f"# {title}", "", f"`cryspy.{package}`", ""]
    for f_name in sorted(os.listdir(dir_package)):
        if f_name.endswith(".py") and not f_name.startswith("__"):
            ls_out += [f"::: cryspy.{package}.{f_name[:-3]}", ""]
    return "\n".join(ls_out)


def _generated():
    """(path, content or source file) of every generated page."""
    l_out = [(path, text) for path, _, text in _help_pages()]
    l_out.append(("examples/index.md", _examples_index()))
    for name in _example_dirs():
        l_image, text = _example_page(name)
        l_out.append((f"examples/{name}.md", text))
        for f_image in l_image:
            l_out.append((f"examples/{name}/{f_image}",
                          File  # marker: copy a binary file
                          ))
    l_out.append(("rcif-categories.md", _categories_page()))
    for package, title in API_PAGES:
        l_out.append((f"api/{package}.md", _api_page(package, title)))
    return l_out


def on_config(config):
    nav = [
        {"Home": "index.md"},
        {"rcif format": ["rcif-format.md", "rcif-categories.md"]},
        {"Guides": [{title: path} for path, title, _ in _help_pages()]},
        {"Examples": ["examples/index.md"] + [
            {EXAMPLE_TITLES[name]: f"examples/{name}.md"} for name in _example_dirs()]},
        {"API reference": [{title: f"api/{package}.md"} for package, title in API_PAGES]},
    ]
    config["nav"] = nav
    return config


def on_files(files, config):
    for path, content in _generated():
        if content is File:
            name, f_image = path[len("examples/"):].split("/")
            files.append(File.generated(
                config, path,
                abs_src_path=os.path.join(DIR_EXAMPLES, name, f_image)))
        else:
            files.append(File.generated(config, path, content=content))
    return files
