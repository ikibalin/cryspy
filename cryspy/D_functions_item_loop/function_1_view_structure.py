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
import warnings

import numpy

from cryspy.A_functions_base.function_3_den_file import read_den_file_grid
from cryspy.A_functions_base.optional_import import optional_module
from cryspy.A_functions_base.structure_scene import (
    calc_box_mesh, calc_ellipsoid_matrix_by_adp,
    calc_ellipsoid_matrix_by_susceptibility, calc_ellipsoid_mesh,
    calc_flags_non_positive_definite, calc_principal_axes_by_tensor,
    form_cube_text)

na = numpy.newaxis

py3Dmol = optional_module("py3Dmol", extra="view3d")
ipywidgets = optional_module("ipywidgets", extra="view3d")

AXIS_COLORS = ("red", "green", "blue")

ATOM_DISPLAY_SPHERES = "spheres"
ATOM_DISPLAY_ADP = "ADP ellipsoids"
ATOM_DISPLAY_SUSCEPTIBILITY = "susceptibility ellipsoids"
ATOM_DISPLAYS = (ATOM_DISPLAY_SPHERES, ATOM_DISPLAY_ADP,
                 ATOM_DISPLAY_SUSCEPTIBILITY)

# atoms with non-positive definite displacement parameters (NPD) are drawn
# with |U|: as boxes (the ORTEP/SHELX/Olex2 convention) or as half
# transparent ellipsoids crossed by bars along the axes with U_i <= 0
NPD_DISPLAY_CUBES = "cubes"
NPD_DISPLAY_MARKED = "ellipsoids with bars"
NPD_DISPLAYS = (NPD_DISPLAY_CUBES, NPD_DISPLAY_MARKED)
NPD_BAR_COLOR = "#FF00FF"
NPD_OPACITY = 0.7
# a small negative U gives a mark hidden inside the bonds: the largest half
# axis of the mark is at least this, Å (the shape keeps its proportions)
NPD_MIN_SEMI_AXIS = 0.4

# Jmol colours (as 3Dmol.js uses them) for the ellipsoids drawn as meshes
JMOL_COLORS = {
    "H": "#FFFFFF", "He": "#D9FFFF", "Li": "#CC80FF", "Be": "#C2FF00",
    "B": "#FFB5B5", "C": "#909090", "N": "#3050F8", "O": "#FF0D0D",
    "F": "#90E050", "Ne": "#B3E3F5", "Na": "#AB5CF2", "Mg": "#8AFF00",
    "Al": "#BFA6A6", "Si": "#F0C8A0", "P": "#FF8000", "S": "#FFFF30",
    "Cl": "#1FF01F", "Ar": "#80D1E3", "K": "#8F40D4", "Ca": "#3DFF00",
    "Sc": "#E6E6E6", "Ti": "#BFC2C7", "V": "#A6A6AB", "Cr": "#8A99C7",
    "Mn": "#9C7AC7", "Fe": "#E06633", "Co": "#F090A0", "Ni": "#50D050",
    "Cu": "#C88033", "Zn": "#7D80B0", "Ga": "#C28F8F", "Ge": "#668F8F",
    "As": "#BD80E3", "Se": "#FFA100", "Br": "#A62929", "Kr": "#5CB8D1",
    "Rb": "#702EB0", "Sr": "#00FF00", "Y": "#94FFFF", "Zr": "#94E0E0",
    "Nb": "#73C2C9", "Mo": "#54B5B5", "Tc": "#3B9E9E", "Ru": "#248F8F",
    "Rh": "#0A7D8C", "Pd": "#006985", "Ag": "#C0C0C0", "Cd": "#FFD98F",
    "In": "#A67573", "Sn": "#668080", "Sb": "#9E63B5", "Te": "#D47A00",
    "I": "#940094", "Xe": "#429EB0", "Cs": "#57178F", "Ba": "#00C900",
    "La": "#70D4FF", "Ce": "#FFFFC7", "Pr": "#D9FFC7", "Nd": "#C7FFC7",
    "Pm": "#A3FFC7", "Sm": "#8FFFC7", "Eu": "#61FFC7", "Gd": "#45FFC7",
    "Tb": "#30FFC7", "Dy": "#1FFFC7", "Ho": "#00FF9C", "Er": "#00E675",
    "Tm": "#00D452", "Yb": "#00BF38", "Lu": "#00AB24", "Hf": "#4DC2FF",
    "Ta": "#4DA6FF", "W": "#2194D6", "Re": "#267DAB", "Os": "#266696",
    "Ir": "#175487", "Pt": "#D0D0E0", "Au": "#FFD123", "Hg": "#B8B8D0",
    "Tl": "#A6544D", "Pb": "#575961", "Bi": "#9E4FB5", "Po": "#AB5C00",
    "At": "#754F45", "Rn": "#428296", "Fr": "#420066", "Ra": "#007D00",
    "Ac": "#70ABFA", "Th": "#00BAFF", "Pa": "#00A1FF", "U": "#008FFF",
    "Np": "#0080FF", "Pu": "#006BFF", "Am": "#545CF2", "Cm": "#785CE3",
}


def _setting(default, description: str, **limits):
    return dataclasses.field(
        default=default, metadata=dict(description=description, **limits))


@dataclasses.dataclass
class StructureViewSettings:
    """Settings of the 3D view of a crystal structure.

    Every field carries ``metadata["description"]``; numbers also
    ``"min"``, ``"max"``, ``"step"`` and maybe ``"unit"``, a choice among
    fixed strings ``"choices"``; Jupyter widgets and the panel of CrysPy
    Editor are built from this class.
    """

    n_cells_a: int = _setting(1, "Number of unit cells along a.", min=1, max=6, step=1)
    n_cells_b: int = _setting(1, "Number of unit cells along b.", min=1, max=6, step=1)
    n_cells_c: int = _setting(1, "Number of unit cells along c.", min=1, max=6, step=1)
    atom_display: str = _setting(
        ATOM_DISPLAY_SPHERES,
        "How atoms are drawn: spheres; ellipsoids of the atomic displacement "
        "parameters (isotropic or anisotropic); magnetization ellipsoids of "
        "the susceptibility tensor for magnetic atoms (the other atoms stay "
        "spheres).",
        choices=ATOM_DISPLAYS)
    atom_scale: float = _setting(
        0.3, "Sphere size as a fraction of the van der Waals radius.",
        min=0.05, max=1., step=0.05)
    adp_probability: float = _setting(
        0.5, "Probability to find the atom inside its displacement ellipsoid.",
        min=0.1, max=0.99, step=0.05)
    npd_display: str = _setting(
        NPD_DISPLAY_CUBES,
        "How atoms with non-positive definite displacement parameters "
        "(negative Biso, Uiso or a negative eigenvalue of Uani) are drawn "
        "with |U|: cubes (boxes along the principal axes), or half "
        "transparent ellipsoids crossed by magenta bars along the axes "
        "with negative U_i; a mark is at least 0.4 Å in size to be seen.",
        choices=NPD_DISPLAYS)
    show_bonds: bool = _setting(True, "Show bonds (found by 3Dmol.js from the distances).")
    show_labels: bool = _setting(False, "Show atom labels.")
    show_cell: bool = _setting(True, "Show the edges of the block of unit cells.")
    show_inner_cells: bool = _setting(False, "Show the edges of every unit cell.")
    show_axes: bool = _setting(True, "Show the a, b, c axes.")
    show_boundary_atoms: bool = _setting(
        True, "Show atoms on the far faces of the block (x, y or z = n).")
    hide_dummy_atoms: bool = _setting(
        True, "Hide dummy atoms (type_symbol 'Dummy', e.g. points of local axes).")
    show_moments: bool = _setting(True, "Show magnetic moments as arrows.")
    moment_scale: float = _setting(
        1., "Length per 1 μB, Å/μB: scale of the moment arrows and of the "
        "susceptibility ellipsoids.",
        min=0.05, max=10., step=0.05, unit="Å/μB")
    susceptibility_field: float = _setting(
        1., "Magnetic field, T: a susceptibility ellipsoid shows the moments "
        "this field induces in every direction (semi-axes |χ|·B).",
        min=0.1, max=50., step=0.1, unit="T")
    moment_color: str = _setting("red", "Colour of the moment arrows.")
    density_file: str = _setting(
        "", "Path of a .den file (e.g. the spin or magnetization density "
        "of MEM) drawn as isosurfaces; empty: no density.")
    density_file_subtracted: str = _setting(
        "", "Path of a second .den file subtracted from the first, e.g. "
        "to see a deformation density; empty: none. Both files need the "
        "same grid.")
    density_level: float = _setting(
        0.1, "Isosurface level as a fraction of the largest |density|: "
        "surfaces at +level and -level.",
        min=0.01, max=1., step=0.01)
    density_opacity: float = _setting(
        0.6, "Opacity of the density isosurfaces.", min=0.1, max=1., step=0.05)
    density_color_positive: str = _setting(
        "blue", "Colour of the isosurface of positive density.")
    density_color_negative: str = _setting(
        "red", "Colour of the isosurface of negative density.")
    background: str = _setting("white", "Background colour.")
    width: int = _setting(
        640, "Width of the view in pixels.", min=200, max=2000, step=10, unit="px")
    height: int = _setting(
        480, "Height of the view in pixels.", min=200, max=2000, step=10, unit="px")

    def __post_init__(self):
        if self.atom_display not in ATOM_DISPLAYS:
            raise ValueError(
                f"atom_display is {self.atom_display!r}, not one of {ATOM_DISPLAYS}.")
        if self.npd_display not in NPD_DISPLAYS:
            raise ValueError(
                f"npd_display is {self.npd_display!r}, not one of {NPD_DISPLAYS}.")

    @property
    def n_cells(self):
        return (self.n_cells_a, self.n_cells_b, self.n_cells_c)


def take_element(type_symbol: str) -> str:
    """Chemical element of an atom type: 'Fe3+' -> 'Fe', 'O2-' -> 'O'."""
    found = re.search(r"[A-Z][a-z]?", str(type_symbol))
    return found.group(0) if found is not None else "X"


def _point(xyz, decimals: int = None) -> dict:
    if decimals is not None:  # short numbers: smaller page for many meshes
        return {"x": round(float(xyz[0]), decimals),
                "y": round(float(xyz[1]), decimals),
                "z": round(float(xyz[2]), decimals)}
    return {"x": float(xyz[0]), "y": float(xyz[1]), "z": float(xyz[2])}


def calc_display_ellipsoids(scene: dict, settings: StructureViewSettings):
    """Matrices (Å) turning a unit sphere into the drawn ellipsoids.

    Parameters
    ----------
    scene : dict
        Output of ``structure_scene.calc_structure_scene``; its
        ``tensor_xyz`` holds U (Å²) or χ (μB/T) as ``settings.atom_display``
        asks.
    settings : StructureViewSettings

    Returns
    -------
    numpy.ndarray, shape (3, 3, n), or None for spheres
        NaN for atoms drawn as spheres. A non-positive definite U gives
        the ellipsoid of |U| (see :func:`calc_display_non_positive_adp`).
    """
    tensor_xyz = scene.get("tensor_xyz")
    if settings.atom_display == ATOM_DISPLAY_SPHERES or tensor_xyz is None:
        return None
    flags = numpy.all(numpy.isfinite(tensor_xyz), axis=(0, 1))
    res = numpy.full(tensor_xyz.shape, numpy.nan)
    if not numpy.any(flags):
        return res
    if settings.atom_display == ATOM_DISPLAY_ADP:
        matrix = calc_ellipsoid_matrix_by_adp(
            tensor_xyz[:, :, flags], settings.adp_probability,
            flag_absolute=True)
    else:
        matrix = calc_ellipsoid_matrix_by_susceptibility(
            tensor_xyz[:, :, flags], settings.susceptibility_field
            ) * settings.moment_scale
    # a zero tensor (e.g. no displacement parameters) stays a sphere
    zero = ~numpy.any(numpy.abs(matrix) > 1e-6, axis=(0, 1))
    matrix[:, :, zero] = numpy.nan
    res[:, :, flags] = matrix
    return res


def calc_display_non_positive_adp(scene: dict, settings: StructureViewSettings):
    """Atoms of the scene with non-positive definite displacement tensors.

    Parameters
    ----------
    scene : dict
        Output of ``structure_scene.calc_structure_scene``.
    settings : StructureViewSettings

    Returns
    -------
    numpy.ndarray, shape (n,), bool
        True for an atom with NPD displacement parameters; all False
        unless ADP ellipsoids are shown.
    """
    tensor_xyz = scene.get("tensor_xyz")
    n_atom = scene["atom_index"].size
    flags = numpy.zeros(n_atom, dtype=bool)
    if settings.atom_display != ATOM_DISPLAY_ADP or tensor_xyz is None:
        return flags
    finite = numpy.all(numpy.isfinite(tensor_xyz), axis=(0, 1))
    if numpy.any(finite):
        flags[finite] = calc_flags_non_positive_definite(tensor_xyz[:, :, finite])
    return flags


def calc_display_density_grid(axes_xyz, settings: StructureViewSettings):
    """Density of the .den files of the settings on the grid of the cell.

    Parameters
    ----------
    axes_xyz : numpy.ndarray, shape (3, 3)
        Cartesian vectors a, b, c of the shown crystal as columns, Å.
    settings : StructureViewSettings
        ``density_file`` minus ``density_file_subtracted`` (if given).

    Returns
    -------
    numpy.ndarray, shape (n_a, n_b, n_c), or None without a density file
        A warning is given when a file has another unit cell than the
        crystal (the density is drawn in the cell of the crystal).
    """
    if not settings.density_file:
        return None
    lengths = numpy.sqrt(numpy.square(numpy.asarray(axes_xyz)).sum(axis=0))
    grid = None
    for file_name, sign in ((settings.density_file, 1.),
                            (settings.density_file_subtracted, -1.)):
        if not file_name:
            continue
        cell_parameters, grid_file = read_den_file_grid(file_name)
        if not numpy.allclose(cell_parameters[:3], lengths, rtol=1e-2):
            warnings.warn(
                f"The cell of {file_name!r} ({cell_parameters[:3]} Å) differs "
                f"from the cell of the crystal ({lengths} Å).", stacklevel=2)
        if grid is None:
            grid = sign * grid_file
        elif grid.shape != grid_file.shape:
            raise ValueError(
                f"The grid of {file_name!r} {grid_file.shape} differs from "
                f"the grid of {settings.density_file!r} {grid.shape}: "
                "the densities cannot be subtracted.")
        else:
            grid = grid + sign * grid_file
    return grid


def _add_density(view, scene: dict, settings: StructureViewSettings):
    """Isosurfaces at +-density_level * max|density| of the .den files."""
    grid = calc_display_density_grid(scene["axes_xyz"], settings)
    if grid is None:
        return
    density_max = numpy.nanmax(numpy.abs(grid))
    if not density_max > 0.:
        return
    level = settings.density_level * density_max
    cube = form_cube_text(scene["axes_xyz"], grid, settings.n_cells)
    for sign, color in ((1., settings.density_color_positive),
                        (-1., settings.density_color_negative)):
        if numpy.nanmax(sign * grid) >= level:
            view.addVolumetricData(cube, "cube", {
                "isoval": float(sign * level), "color": color,
                "opacity": settings.density_opacity})


def _add_mesh(view, vertices, normals, faces, color: str, opacity: float = 1.):
    shape = {
        "vertexArr": [_point(v, 3) for v in vertices.transpose()],
        "normalArr": [_point(n, 3) for n in normals.transpose()],
        "faceArr": [int(i) for i in faces.flatten()],
        "color": color}
    if opacity < 1.:
        shape["opacity"] = opacity
    view.addCustom(shape)


def _add_non_positive_atom(view, xyz, matrix, u_xyz, color: str,
                           settings: StructureViewSettings):
    """A cube, or an ellipsoid with bars along the negative axes, of |U|."""
    # the drawn surface c |U|^(1/2) along the principal axes of U
    val, vec = calc_principal_axes_by_tensor(u_xyz[:, :, na])
    val, vec = val[:, 0], vec[:, :, 0]
    semi = numpy.sqrt(numpy.square(matrix @ vec).sum(axis=0))
    # a zero eigenvalue would give a flat box: keep a tenth of the largest axis
    semi = numpy.maximum(semi, 0.1 * semi.max())
    semi = semi * max(1., NPD_MIN_SEMI_AXIS / semi.max())
    axes = vec * semi[na, :]
    if settings.npd_display == NPD_DISPLAY_CUBES:
        _add_mesh(view, *calc_box_mesh(xyz, axes), color)
        return
    _add_mesh(view, *calc_ellipsoid_mesh(xyz, axes, n_theta=10, n_phi=16),
              color, opacity=NPD_OPACITY)
    for i_axis in numpy.flatnonzero(val <= 1e-9 * numpy.abs(val).max()):
        half = 1.3 * axes[:, i_axis]
        view.addCylinder({"start": _point(xyz - half), "end": _point(xyz + half),
                          "radius": 0.04, "color": NPD_BAR_COLOR,
                          "fromCap": 1, "toCap": 1})


def make_structure_view(  # noqa: C901, PLR0912, PLR0915
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
    ellipsoid_matrix = calc_display_ellipsoids(scene, settings)
    if ellipsoid_matrix is None:
        flags_ellipsoid = numpy.zeros(index.size, dtype=bool)
    else:
        flags_ellipsoid = ~numpy.any(numpy.isnan(ellipsoid_matrix), axis=(0, 1))
    flags_npd = calc_display_non_positive_adp(scene, settings) & flags_ellipsoid
    stick = {"radius": 0.1, "colorscheme": "Jmol"}
    style = {"sphere": {"scale": settings.atom_scale, "colorscheme": "Jmol"}}
    if settings.show_bonds:
        style["stick"] = stick
    view.setStyle({}, style)
    if numpy.any(flags_ellipsoid):
        # 3Dmol.js has no ellipsoids: the atom keeps only its bonds and the
        # ellipsoid is a mesh
        style_ellipsoid = ({"stick": stick} if settings.show_bonds else
                           {"sphere": {"hidden": True}})
        view.setStyle({"index": [int(i) for i in numpy.flatnonzero(flags_ellipsoid)]},
                      style_ellipsoid)
        for i_scene in numpy.flatnonzero(flags_ellipsoid):
            xyz = atom_xyz[:, i_scene]
            matrix = ellipsoid_matrix[:, :, i_scene]
            color = JMOL_COLORS.get(atom_element[index[i_scene]], "#FF1493")
            if flags_npd[i_scene]:
                _add_non_positive_atom(view, xyz, matrix,
                                       scene["tensor_xyz"][:, :, i_scene],
                                       color, settings)
                continue
            _add_mesh(view, *calc_ellipsoid_mesh(xyz, matrix, n_theta=10, n_phi=16),
                      color)

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
        for xyz, moment in zip(atom_xyz.transpose(), moment_xyz.transpose()):
            if numpy.any(numpy.isnan(moment)) or not numpy.any(moment):
                continue
            half = 0.5 * settings.moment_scale * moment
            view.addArrow({"start": _point(xyz - half), "end": _point(xyz + half),
                           "radius": 0.08, "radiusRatio": 2.,
                           "color": settings.moment_color})
    _add_density(view, scene, settings)
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
        if "unit" in meta:
            kwargs["description"] = f"{description}, {meta['unit']}"
        if "choices" in meta:
            widget = ipywidgets.Dropdown(options=list(meta["choices"]), **kwargs)
        elif isinstance(value, bool):
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
