"""Interactive 3D view of a crystal structure with py3Dmol (3Dmol.js).

py3Dmol is an optional extra: ``pip install cryspy[view3d]``. The scene
(atoms in the chosen block of unit cells, moments, cell edges) is computed
by ``cryspy.A_functions_base.structure_scene``; 3Dmol.js draws the atoms
with Jmol colours and van der Waals radii and finds the bonds itself.

The view works in Jupyter notebooks (the returned ``py3Dmol.view`` shows
itself) and in any web view, e.g. CrysPy Editor, through
:func:`make_structure_view_html`.
"""
import dataclasses
import re

import numpy

from cryspy.A_functions_base.optional_import import optional_module

py3Dmol = optional_module("py3Dmol", extra="view3d")
ipywidgets = optional_module("ipywidgets", extra="view3d")

AXIS_COLORS = ("red", "green", "blue")


def _setting(default, description: str, **limits):
    return dataclasses.field(
        default=default, metadata=dict(description=description, **limits))


@dataclasses.dataclass
class StructureViewSettings:
    """Settings of the 3D view of a crystal structure.

    Every field carries ``metadata["description"]`` and, for numbers,
    ``"min"``, ``"max"`` and ``"step"``, so that Jupyter widgets and the
    panel of CrysPy Editor are built from this class.
    """

    n_cells_a: int = _setting(1, "Number of unit cells along a.", min=1, max=6, step=1)
    n_cells_b: int = _setting(1, "Number of unit cells along b.", min=1, max=6, step=1)
    n_cells_c: int = _setting(1, "Number of unit cells along c.", min=1, max=6, step=1)
    atom_scale: float = _setting(
        0.3, "Atom size as a fraction of the van der Waals radius.",
        min=0.05, max=1., step=0.05)
    show_bonds: bool = _setting(True, "Show bonds (found by 3Dmol.js from the distances).")
    show_labels: bool = _setting(False, "Show atom labels.")
    show_cell: bool = _setting(True, "Show the edges of the block of unit cells.")
    show_inner_cells: bool = _setting(False, "Show the edges of every unit cell.")
    show_axes: bool = _setting(True, "Show the a, b, c axes.")
    show_boundary_atoms: bool = _setting(
        True, "Show atoms on the far faces of the block (x, y or z = n).")
    show_moments: bool = _setting(True, "Show magnetic moments as arrows.")
    moment_length: float = _setting(
        2.5, "Arrow length of the largest magnetic moment, Å (the others in proportion).",
        min=0.5, max=10., step=0.1)
    moment_color: str = _setting("red", "Colour of the moment arrows.")
    background: str = _setting("white", "Background colour.")
    width: int = _setting(640, "Width of the view in pixels.", min=200, max=2000, step=10)
    height: int = _setting(480, "Height of the view in pixels.", min=200, max=2000, step=10)

    @property
    def n_cells(self):
        return (self.n_cells_a, self.n_cells_b, self.n_cells_c)


def take_element(type_symbol: str) -> str:
    """Chemical element of an atom type: 'Fe3+' -> 'Fe', 'O2-' -> 'O'."""
    found = re.search(r"[A-Z][a-z]?", str(type_symbol))
    return found.group(0) if found is not None else "X"


def _point(xyz) -> dict:
    return {"x": float(xyz[0]), "y": float(xyz[1]), "z": float(xyz[2])}


def make_structure_view(
        scene: dict, atom_element, atom_label,
        settings: StructureViewSettings = None, js: str = None):
    """Draw a structure scene with py3Dmol.

    Parameters
    ----------
    scene : dict
        Output of ``structure_scene.calc_structure_scene``.
    atom_element : sequence of str
        Chemical element of every atom of the asymmetric unit.
    atom_label : sequence of str
        Label of every atom of the asymmetric unit.
    settings : StructureViewSettings
        Settings of the view; the defaults if None.
    js : str
        URL of 3Dmol-min.js (e.g. a local file for an offline program);
        the py3Dmol default (CDN) if None.

    Returns
    -------
    py3Dmol.view
    """
    if settings is None:
        settings = StructureViewSettings()
    kwargs = {"width": settings.width, "height": settings.height}
    if js is not None:
        kwargs["js"] = js
    view = py3Dmol.view(**kwargs)

    index = scene["atom_index"]
    atom_xyz = scene["atom_xyz"]
    lines = [str(index.size), "cryspy"]
    for i_atom, xyz in zip(index, atom_xyz.transpose()):
        lines.append(f"{atom_element[i_atom]} {xyz[0]:.5f} {xyz[1]:.5f} {xyz[2]:.5f}")
    view.addModel("\n".join(lines), "xyz")
    style = {"sphere": {"scale": settings.atom_scale, "colorscheme": "Jmol"}}
    if settings.show_bonds:
        style["stick"] = {"radius": 0.1, "colorscheme": "Jmol"}
    view.setStyle({}, style)

    if settings.show_cell:
        for start, end in zip(scene["edges_xyz"][0].transpose(),
                              scene["edges_xyz"][1].transpose()):
            view.addCylinder({"start": _point(start), "end": _point(end),
                              "radius": 0.02, "color": "gray",
                              "fromCap": 1, "toCap": 1})
    if settings.show_axes:
        origin = numpy.zeros(3)
        for name, vector, color in zip(
                "abc", scene["axes_xyz"].transpose(), AXIS_COLORS):
            end = 0.3 * vector
            view.addArrow({"start": _point(origin), "end": _point(end),
                           "radius": 0.06, "color": color})
            view.addLabel(name, {"position": _point(end * 1.1),
                                 "fontColor": color, "fontSize": 14,
                                 "showBackground": False, "inFront": True})
    moment_xyz = scene["moment_xyz"]
    if settings.show_moments and moment_xyz is not None:
        norm = numpy.sqrt(numpy.square(moment_xyz).sum(axis=0))
        norm_max = numpy.nanmax(norm, initial=0.)
        for xyz, moment in zip(atom_xyz.transpose(), moment_xyz.transpose()):
            if numpy.any(numpy.isnan(moment)) or not numpy.any(moment):
                continue
            half = 0.5 * settings.moment_length / norm_max * moment
            view.addArrow({"start": _point(xyz - half), "end": _point(xyz + half),
                           "radius": 0.08, "radiusRatio": 2.,
                           "color": settings.moment_color})
    if settings.show_labels:
        for i_atom, xyz in zip(index, atom_xyz.transpose()):
            view.addLabel(str(atom_label[i_atom]), {
                "position": _point(xyz), "fontSize": 10, "fontColor": "black",
                "backgroundColor": "white", "backgroundOpacity": 0.5})
    view.setBackgroundColor(settings.background)
    # c up (z || c), looking slightly from above and from the side
    view.rotate(-80, "x")
    view.rotate(-20, "z")
    view.zoomTo()
    return view


def make_structure_view_html(view, fullpage: bool = True) -> str:
    """Self-contained HTML page of a py3Dmol view (for a web view widget).

    Parameters
    ----------
    view : py3Dmol.view
        The view made by :func:`make_structure_view`.
    fullpage : bool
        Fill the whole page instead of the width and height of the view.
    """
    html = view._make_html()
    if fullpage:
        html = re.sub(r"width: (\S+);", "width: 100%;", html, count=1)
        html = re.sub(r"height: (\S+);", "height: 100vh;", html, count=1)
    return ('<!DOCTYPE html>\n<html>\n<head><meta charset="utf-8"></head>\n'
            '<body style="margin: 0; padding: 0; overflow: hidden;">\n'
            f"{html}\n</body>\n</html>\n")


def make_settings_widgets(settings: StructureViewSettings) -> dict:
    """ipywidgets for every field of the settings (for Jupyter)."""
    widgets = {}
    for field in dataclasses.fields(settings):
        value = getattr(settings, field.name)
        meta = field.metadata
        description = field.name.replace("_", " ")
        kwargs = {"value": value, "description": description,
                  "tooltip": meta.get("description", ""),
                  "style": {"description_width": "initial"}}
        if isinstance(value, bool):
            widget = ipywidgets.Checkbox(**kwargs)
        elif isinstance(value, int):
            widget = ipywidgets.IntSlider(
                min=meta.get("min", 0), max=meta.get("max", 100),
                step=meta.get("step", 1), **kwargs)
        elif isinstance(value, float):
            widget = ipywidgets.FloatSlider(
                min=meta.get("min", 0.), max=meta.get("max", 1.),
                step=meta.get("step", 0.01), **kwargs)
        else:
            widget = ipywidgets.Text(**kwargs)
        widgets[field.name] = widget
    return widgets
