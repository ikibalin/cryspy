"""Lazy import of optional dependencies.

Visualization libraries (matplotlib, ...) are optional extras of cryspy:
``pip install cryspy`` installs only the calculation part and
``pip install cryspy[plot]`` adds matplotlib for the ``plot_*`` methods.
Modules import them through :func:`optional_module`, so ``import cryspy``
works without them and the error appears only when a plot is requested.
"""
import importlib


class MissingOptionalDependency(ImportError):
    """An optional dependency of cryspy is not installed."""


class OptionalModule:
    """Module proxy that imports the module on first attribute access.

    Submodules are imported on access as well, so ``matplotlib.colors``
    works on a proxy of ``matplotlib``.
    """

    def __init__(self, name: str, extra: str):
        self._name = name
        self._extra = extra
        self._module = None

    def _load(self):
        if self._module is None:
            try:
                self._module = importlib.import_module(self._name)
            except ImportError as error:
                package = self._name.split(".")[0]
                raise MissingOptionalDependency(
                    f"'{package}' is required for this function but it is "
                    "not installed. Install it with "
                    f"'pip install cryspy[{self._extra}]'."
                ) from error
        return self._module

    def __getattr__(self, attr: str):
        if attr.startswith("_"):
            raise AttributeError(attr)
        module = self._load()
        try:
            return getattr(module, attr)
        except AttributeError:
            return importlib.import_module(f"{self._name}.{attr}")

    def __repr__(self) -> str:
        return f"<optional module '{self._name}' (cryspy[{self._extra}])>"


def optional_module(name: str, extra: str = "plot") -> OptionalModule:
    """Return a lazy proxy of an optional module.

    Parameters
    ----------
    name : str
        Full module name, e.g. ``"matplotlib.pyplot"``.
    extra : str
        Name of the cryspy extra that installs the module.
    """
    return OptionalModule(name, extra)


def is_available(name: str) -> bool:
    """Check whether an optional module can be imported."""
    try:
        importlib.import_module(name)
    except ImportError:
        return False
    return True
