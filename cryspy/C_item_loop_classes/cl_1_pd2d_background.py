"""Pd2dBackground class."""
import numpy
import scipy
import scipy.interpolate
from typing import NoReturn, Union
from cryspy.A_functions_base.function_1_strings import \
    string_to_value_error, value_error_to_string
from cryspy.B_parent_classes.cl_1_item import ItemN


class Pd2dBackground(ItemN):
    """
    Pd2dBackground class.

    PdInstrReflexAsymmetry describes asymmetry of Bragg reflections for
    1d powder diffractometer.

    Attributes
    ----------
        - ttheta_phi_intensity

    Internal Attributes
    -------------------
        - ttheta, phi, intensity, intensity_sigma
        - intensity_refinement, intensity_constraint

    Methods
    -------
        - get_variable_names
        - get_variable_by_name, set_variable_by_name
        - form_ttheta_phi_intensity
        - interpolate_by_points
    """

    ATTR_MANDATORY_NAMES = ("ttheta_phi_intensity", )
    ATTR_MANDATORY_TYPES = (str, )
    # ("matrix", "matrix", "matrix", "matrix")
    ATTR_MANDATORY_CIF = ("2theta_phi_intensity", )

    ATTR_OPTIONAL_NAMES = ()
    ATTR_OPTIONAL_TYPES = ()
    ATTR_OPTIONAL_CIF = ()

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ("ttheta", "phi", "intensity", "intensity_sigma",
                      "intensity_refinement", "intensity_constraint")
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered are refined parameters
    ATTR_REF = ()
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])

    # constraints on the parameters
    D_CONSTRAINTS = {}

    # default values for the parameters
    D_DEFAULT = {}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False

    PREFIX = "pd2d_background"

    def __init__(self, **kwargs) -> NoReturn:
        super(Pd2dBackground, self).__init__()

        # defined for any integer and float parameters
        D_MIN = {}

        # defined for ani integer and float parameters
        D_MAX = {}

        self.__dict__["D_MIN"] = D_MIN
        self.__dict__["D_MAX"] = D_MAX
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)

    def form_object(self):
        """Form object."""
        l_1 = (self.ttheta_phi_intensity).strip().split("\n")

        l_ttheta = [float(_) for _ in l_1[0].strip().split()[1:]]
        l_phi, ll_intensity = [], []
        for line in l_1[1:]:
            l_1 = line.strip().split()
            l_phi.append(float(l_1[0]))
            ll_intensity.append(l_1[1:])
        ll_intensity = [[string_to_value_error(ll_intensity[_2][_1])
                         for _2 in range(len(ll_intensity))]
                        for _1 in range(len(ll_intensity[0]))]
        np_int_sigma = numpy.array(ll_intensity, dtype=float)
        self.__dict__["ttheta"] = l_ttheta
        self.__dict__["phi"] = l_phi
        self.__dict__["intensity"] = np_int_sigma[:, :, 0]
        self.__dict__["intensity_sigma"] = numpy.where(
            numpy.isnan(np_int_sigma[:, :, 1]), 0., np_int_sigma[:, :, 1])
        self.__dict__["intensity_refinement"] = numpy.where(
            numpy.isnan(np_int_sigma[:, :, 1]), False, True)
        self.__dict__["intensity_constraint"] = numpy.zeros(
            shape=np_int_sigma[:, :, 0].shape, dtype=bool)

    def form_ttheta_phi_intensity(self) -> NoReturn:
        """Form 2theta_phi_intensity from internal attributes."""
        if ((self.phi is not None) & (self.ttheta is not None) &
                (self.intensity is not None)):
            ls_out = []
            ls_out.append(f"{len(self.phi):12} " + " ".join(
                [f"{_:6.2f}      " for _ in self.ttheta]))
            ll_intensity = self.intensity
            ll_intensity_sigma = self.intensity_sigma
            ll_intensity_refinement = self.intensity_refinement

            ll_intensity = [[ll_intensity[_2][_1] for _2 in
                             range(len(ll_intensity))] for _1 in
                            range(len(ll_intensity[0]))]
            ll_intensity_sigma = [[ll_intensity_sigma[_2][_1] for _2 in
                                   range(len(ll_intensity_sigma))] for _1 in
                                  range(len(ll_intensity_sigma[0]))]
            ll_intensity_refinement = [[ll_intensity_refinement[_2][_1] for _2
                                        in range(len(ll_intensity_refinement))]
                                       for _1 in range(len(
                                               ll_intensity_refinement[0]))]
            for phi, l_int, l_int_sig, l_int_ref in \
                zip(self.phi, ll_intensity, ll_intensity_sigma,
                    ll_intensity_refinement):
                ls_out.append("{:12.2f} ".format(phi) +
                              " ".join(
                                  [f"{value_error_to_string(int, int_sig):12}"
                                   if int_ref else f"{int:12}"
                                   for int, int_sig, int_ref in
                                   zip(l_int, l_int_sig, l_int_ref)]))

            self.__dict__["ttheta_phi_intensity"] = "\n".join(ls_out)

    def get_variable_names(self) -> list:
        """
        Get names of variable as a list.

        (((#prefix, #NAME), (#attribute, (#index_1, #index_2)))

        Returns
        -------
        list
            List of names of variable.

        """
        prefix = self.PREFIX
        np_i, np_j = numpy.where(self.intensity_refinement == True)
        return [((prefix, None), ("intensity", (i, j)))
                for i, j in zip(np_i, np_j)]

    def get_variable_by_name(self, name: tuple) -> Union[float, int, str]:
        """
        Get variable given by name.

        Parameters
        ----------
        name : tuple
            (((#prefix, ), (#attribute, (#index_1, #index_2)))

        Returns
        -------
        Union[float, int, str]
            Value.

        """
        prefix_t = name[0]
        if prefix_t[0] != self.PREFIX:
            return None
        if len(name) == 1:
            return self
        attr_t = name[1]
        attr_name, ind_ij = attr_t
        return getattr(self, attr_name)[ind_ij]

    def set_variable_by_name(self, name: tuple, value) -> NoReturn:
        """
        Set value to variable given by name.

        Parameters
        ----------
        name : tuple
            (((#prefix, ), (#attribute, ))
        value : TYPE
            DESCRIPTION.

        Returns
        -------
        NoReturn

        """
        prefix_t, attr_t = name
        if prefix_t[0] != self.PREFIX:
            return
        attr_name, ind_ij = attr_t
        np_val = getattr(self, attr_name)
        np_val[ind_ij] = value

    def interpolate_by_points(self, tth, phi):
        """Interpolate by points."""
        l_phi_b = self.phi
        l_tth_b = self.ttheta
        ll_int_b = self.intensity
        ll_int_b = [[float(ll_int_b[_2][_1]) for _2 in range(len(ll_int_b))]
                    for _1 in range(len(ll_int_b[0]))]
        if len(l_tth_b) == 0:
            int_2d = numpy.zeros((tth.size, phi.size), dtype=float)
        else:
            phi_b = numpy.array(l_phi_b, dtype=float)
            tth_b = numpy.array(l_tth_b, dtype=float)
            int_b = numpy.array(ll_int_b, dtype=float)
            func = scipy.interpolate.interp2d(tth_b, phi_b, int_b)
            # tth_2d, phi_2d = numpy.meshgrid(tth, phi, indexing="ij")
            int_2d = func(tth, phi)
            int_2d = int_2d.transpose()
        return int_2d

    def is_variables(self):
        """
        Redefine function
        """
        return numpy.any(self.intensity_refinement)

# s_cont = """
#   _pd2d_background_2theta_phi_intensity
#   ;
#       2    4.5     40.0     80.0
#   -3.000 -350    -350.0(15)   -400.0
#   41.000 -350.0   -347(15)   -400.0
#   ;

# """

# obj = Pd2dBackground.from_cif(s_cont)
# print(obj, end="\n\n")
# print(obj.ttheta, end="\n\n")
# print(obj.phi, end="\n\n")
# print(obj.intensity, type(obj.intensity), obj.intensity.dtype, end="\n\n")
# print(obj.intensity_sigma, type(obj.intensity_sigma),
#       obj.intensity_sigma.dtype, end="\n\n")
# print(obj.intensity_refinement, type(obj.intensity_refinement),
#       obj.intensity_refinement.dtype, end="\n\n")
# print(obj.intensity_constraint, type(obj.intensity_constraint),
#       obj.intensity_constraint.dtype, end="\n\n")
# for name in obj.get_variable_names():
#     obj.set_variable_by_name(name, 17)

# obj.form_ttheta_phi_intensity()
# print(obj, end="\n\n")
