"""Descriptors of the attributes of CIF items.

A class derived from ItemN describes its CIF attributes by the tuples
ATTR_NAMES, ATTR_TYPES, ATTR_REF, ... (they are also read by CrysPy Editor
and by the AI assistants). ItemN.__init_subclass__ turns every name of them
into a descriptor of the class (install_cif_descriptors):

- CifValue for ATTR_NAMES: converts the value to its type (ATTR_TYPES),
  checks D_CONSTRAINTS, clips to D_MIN / D_MAX and drops the internal
  parameters calculated from the old value;
- CifSigma (ATTR_SIGMA), CifFlag (ATTR_CONSTR_FLAG, ATTR_REF_FLAG),
  CifMark (ATTR_CONSTR_MARK);
- CifPlain (ATTR_INT_PROTECTED_NAMES): '.' and '?' only become None.

The value is kept in the __dict__ of the object under the same name, so
the code writing ``self.__dict__[name]`` directly keeps working. An
attribute that is not given raises AttributeError, as before
(``is_attribute`` relies on it).
"""
import warnings

import numpy

from cryspy.exceptions import CrysPyWarning, RcifFormatError


def normalize_cif_value(value):
    """None for '.', '?' and None; the element of a numpy array of size 1."""
    if isinstance(value, str) and value in (".", "?"):
        return None
    if isinstance(value, numpy.ndarray) and value.size == 1:
        return value[0]
    return value


class CifPlain:
    """Attribute stored in the __dict__ of the object: '.' and '?' are None."""

    __slots__ = ("name",)

    def __init__(self, name: str):
        self.name = name

    def __get__(self, instance, owner=None):
        if instance is None:
            return self
        try:
            return instance.__dict__[self.name]
        except KeyError:
            raise AttributeError(
                f"Attribute '{self.name:}' is not defined in "
                f"'{type(instance).__name__:}'") from None

    def __set__(self, instance, value):
        instance.__dict__[self.name] = self.convert(
            instance, normalize_cif_value(value))

    def __delete__(self, instance):
        try:
            del instance.__dict__[self.name]
        except KeyError:
            raise AttributeError(self.name) from None

    def convert(self, instance, value):
        return value

    def __repr__(self):
        return f"{type(self).__name__:}('{self.name:}')"


class CifSigma(CifPlain):
    """Standard uncertainty of a refinable attribute (float)."""

    __slots__ = ()

    def convert(self, instance, value):
        return None if value is None else float(value)


class CifFlag(CifPlain):
    """Refinement or constraint flag of a refinable attribute (bool)."""

    __slots__ = ()

    def convert(self, instance, value):
        return None if value is None else bool(value)


class CifMark(CifPlain):
    """Constraint mark of a refinable attribute (str)."""

    __slots__ = ()

    def convert(self, instance, value):
        return None if value is None else str(value).strip()


class CifValue(CifPlain):
    """CIF attribute with the type ``value_type``."""

    __slots__ = ("value_type",)

    def __init__(self, name: str, value_type=str):
        super().__init__(name)
        self.value_type = value_type

    def __set__(self, instance, value):
        value = normalize_cif_value(value)
        if value is None:
            instance.__dict__[self.name] = None
            return
        name, value_type = self.name, self.value_type
        if value_type is bool and isinstance(value, str):
            value_new = value == "True"
        else:
            try:
                value_new = value_type(value)
            except (TypeError, ValueError) as error_value:
                raise RcifFormatError(
                    f"Value '{value:}' of '{name:}' in {type(instance).__name__:} "
                    f"is not {value_type.__name__:}.") from error_value
        flag_write = True
        d_constraints = take_class_or_object_dict(instance, "D_CONSTRAINTS")
        if name in d_constraints and value_new not in d_constraints[name]:
            flag_write = False
            warnings.warn(
                f"Value '{value_new:}' of '{name:}' in {type(instance).__name__:} "
                f"is ignored: allowed values are "
                f"{', '.join(str(_) for _ in d_constraints[name]):}.",
                CrysPyWarning, stacklevel=2)
        d_min = take_class_or_object_dict(instance, "D_MIN")
        if name in d_min and value_new < d_min[name]:
            value_new = d_min[name]
        d_max = take_class_or_object_dict(instance, "D_MAX")
        if name in d_max and value_new > d_max[name]:
            value_new = d_max[name]
        # the internal parameters were calculated from the old value
        instance.delete_internal_parameters()
        if flag_write:
            instance.__dict__[name] = value_new

    def __repr__(self):
        return f"CifValue('{self.name:}', {getattr(self.value_type, '__name__', self.value_type)})"


def take_class_or_object_dict(instance, name: str) -> dict:
    """D_MIN, D_MAX, D_CONSTRAINTS of the object (set in __init__) or of
    its class; {} when neither has it."""
    value = instance.__dict__.get(name)
    if value is None:
        value = getattr(type(instance), name, None)
    return value if isinstance(value, dict) else {}


def install_cif_descriptors(cls):
    """Descriptors for every name of the ATTR_* tuples of the class cls.

    A name the class (or a parent other than ItemN) defines as a method
    or a property is left as it is.
    """
    l_type = list(getattr(cls, "ATTR_TYPES", ()))
    d_descriptor = {}
    for name in getattr(cls, "ATTR_INT_PROTECTED_NAMES", ()):
        d_descriptor[name] = CifPlain(name)
    for ind, name in enumerate(getattr(cls, "ATTR_NAMES", ())):
        value_type = l_type[ind] if ind < len(l_type) else str
        d_descriptor[name] = CifValue(name, value_type)
    for name in getattr(cls, "ATTR_SIGMA", ()):
        d_descriptor[name] = CifSigma(name)
    for name in (tuple(getattr(cls, "ATTR_CONSTR_FLAG", ())) +
                 tuple(getattr(cls, "ATTR_REF_FLAG", ()))):
        d_descriptor[name] = CifFlag(name)
    for name in getattr(cls, "ATTR_CONSTR_MARK", ()):
        d_descriptor[name] = CifMark(name)
    for name, descriptor in d_descriptor.items():
        existing = None
        for parent in cls.__mro__:
            if name in parent.__dict__:
                existing = parent.__dict__[name]
                break
        if existing is not None and not isinstance(existing, CifPlain):
            continue  # a method, property or class constant of that name
        setattr(cls, name, descriptor)
