"""Pd2dBackground class."""
import numpy
import scipy
import scipy.interpolate
from typing import NoReturn, Union

from cryspy.A_functions_base.function_1_strings import \
    string_to_value_error_mark, value_error_mark_to_string

from cryspy.B_parent_classes.cl_1_item import ItemN

from cryspy.C_item_loop_classes.cl_1_pd2d_meas import Pd2dMeas

na = numpy.newaxis

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

    ATTR_MANDATORY_NAMES = ()
    ATTR_MANDATORY_TYPES = ()
    # ("matrix", "matrix", "matrix", "matrix")
    ATTR_MANDATORY_CIF = ()

    ATTR_OPTIONAL_NAMES = ("gamma_nu_intensity", )
    ATTR_OPTIONAL_TYPES = (str, str)
    ATTR_OPTIONAL_CIF = ("gamma_nu_intensity", )

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ("intensity", "intensity_sigma",
                      "intensity_refinement", "intensity_constraint",
                      "gamma", "nu")
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered are refined parameters
    ATTR_REF = ()
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    # constraints on the parameters
    D_CONSTRAINTS = {}

    # default values for the parameters
    D_DEFAULT = {}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False
    for key in ATTR_CONSTR_MARK:
        D_DEFAULT[key] = ""

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
        self_keys = self.__dict__.keys()

        l_1 = (self.gamma_nu_intensity).strip().split("\n")
        l_gamma = [float(_) for _ in l_1[0].strip().split()[1:]]
        l_nu, ll_intensity = [], []
        for line in l_1[1:]:
            l_1 = line.strip().split()
            l_nu.append(float(l_1[0]))
            ll_intensity.append(l_1[1:])
        ll_intensity_2 = [[string_to_value_error_mark(ll_intensity[_2][_1])[:2]
                         for _2 in range(len(ll_intensity))]
                        for _1 in range(len(ll_intensity[0]))]
        np_int_sigma = numpy.array(ll_intensity_2, dtype=float)
        ll_mark = [[string_to_value_error_mark(ll_intensity[_2][_1])[2]
                         for _2 in range(len(ll_intensity))]
                        for _1 in range(len(ll_intensity[0]))]
        np_int_mark = numpy.array(ll_mark, dtype=str)
        self.__dict__["gamma"] = l_gamma
        self.__dict__["nu"] = l_nu
        self.__dict__["intensity"] = np_int_sigma[:, :, 0]
        self.__dict__["intensity_sigma"] = numpy.where(
            numpy.isnan(np_int_sigma[:, :, 1]), 0., np_int_sigma[:, :, 1])
        self.__dict__["intensity_refinement"] = numpy.where(
            numpy.isnan(np_int_sigma[:, :, 1]), False, True)
        self.__dict__["intensity_constraint"] = numpy.zeros(
            shape=np_int_sigma[:, :, 0].shape, dtype=bool)
        self.__dict__["intensity_mark"] = np_int_mark

    def form_gamma_nu_intensity(self) -> NoReturn:
        """Form 2theta_phi_intensity from internal attributes."""
        if ((self.nu is not None) & (self.gamma is not None) &
                (self.intensity is not None)):
            ls_out = []
            ls_out.append(f"{len(self.nu):12} " + " ".join(
                [f"{_:6.2f}      " for _ in self.gamma]))
            ll_intensity = self.intensity
            ll_intensity_sigma = self.intensity_sigma
            ll_intensity_refinement = self.intensity_refinement
            ll_intensity_mark = self.intensity_mark

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
            ll_intensity_mark = [[ll_intensity_mark[_2][_1] for _2
                                        in range(len(ll_intensity_mark))]
                                       for _1 in range(len(
                                               ll_intensity_mark[0]))]
            for nu, l_int, l_int_sig, l_int_ref, l_int_mark in \
                zip(self.nu, ll_intensity, ll_intensity_sigma,
                    ll_intensity_refinement, ll_intensity_mark):
                ls_out.append("{:12.2f} ".format(nu) +
                              " ".join(
                                  [f"{value_error_mark_to_string(int, int_sig, int_mark):12}"
                                   if int_ref else f"{int:12}"
                                   for int, int_sig, int_ref, int_mark in
                                   zip(l_int, l_int_sig, l_int_ref, l_int_mark)]))

            self.__dict__["gamma_nu_intensity"] = "\n".join(ls_out)

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
        if self.is_attribute("intensity"):
            np_i, np_j = numpy.where(self.intensity_refinement == True)
            return [((prefix, None), ("intensity", (i, j)))
                    for i, j in zip(np_i, np_j)]
        else:
            return []

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
        l_phi_b = self.nu
        l_tth_b = self.gamma
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


    def define_points(self, pd2d_meas: Pd2dMeas, step_gamma: float = 10., step_nu: float = 10.):
        gamma = pd2d_meas.gamma
        nu = pd2d_meas.nu
        if pd2d_meas.is_attribute("gamma_nu_intensity_plus"):
            intensity = (pd2d_meas.intensity_plus + pd2d_meas.intensity_minus)
        else:
            intensity = pd2d_meas.intensity

        points_tth = int((gamma.max()-gamma.min())/step_gamma + 2)
        points_phi = int((nu.max()-nu.min())/step_nu + 2)

        ttheta_bkgr = numpy.linspace(gamma.min(), gamma.max(), points_tth, endpoint=True)
        phi_bkgr = numpy.linspace(nu.min(), nu.max(), points_phi, endpoint=True)

        flags_tth = numpy.abs(gamma[na, :]-ttheta_bkgr[:, na])<step_gamma
        flags_phi = numpy.abs(nu[na, :]-phi_bkgr[:, na])<step_nu
        flags = flags_tth[:, :, na, na] * flags_phi[na, na, :, :]
        int_bkrg = numpy.zeros((points_tth, points_phi), dtype=float)
        for i_tth in range(points_tth):
            for i_phi in range(points_phi):
                flag = flags[i_tth, :, i_phi, :]
                flag_nan = numpy.isnan(intensity[flag])
                if numpy.all(flag_nan):
                    int_bkrg[i_tth, i_phi] = 0.
                else:
                    int_bkrg[i_tth, i_phi] = numpy.nanmin(intensity[flag][numpy.logical_not(flag_nan)])
        self.gamma = ttheta_bkgr
        self.nu = phi_bkgr
        self.intensity = numpy.round(int_bkrg, decimals=5)
        self.intensity_refinement = numpy.zeros(int_bkrg.shape, dtype=bool)
        self.intensity_sigma = numpy.zeros(int_bkrg.shape, dtype=float)
        self.intensity_mark = numpy.zeros(int_bkrg.shape, dtype=str)
        self.form_gamma_nu_intensity()
        return

    def get_dictionary(self):
        res = {}
        res["background_gamma"] = numpy.array(self.gamma, dtype=float)*numpy.pi/180.
        res["background_nu"] = numpy.array(self.nu, dtype=float)*numpy.pi/180.
        res["background_intensity"] = numpy.array(self.intensity, dtype=float)
        res["flags_background_intensity"] = numpy.array(self.intensity_refinement, dtype=bool)
        return res