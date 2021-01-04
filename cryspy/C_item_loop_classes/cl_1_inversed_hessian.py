"""Pd2dBackground class."""
import numpy
import scipy
import scipy.interpolate
from typing import NoReturn, Union

from cryspy.B_parent_classes.cl_1_item import ItemN

def hessian_to_correlation(np_hessian):
    """Calculate correlation matrix and sigmas for hessian matrix."""
    np_sigma_sq = numpy.diagonal(np_hessian)
    np_na = numpy.newaxis
    np_sigma = numpy.sqrt(numpy.abs(np_sigma_sq))
    np_sigma_sq_matrix = np_sigma[:, np_na] * np_sigma[np_na, :]
    np_correlation_matrix = np_hessian / np_sigma_sq_matrix
    return np_correlation_matrix, np_sigma

class InversedHessian(ItemN):
    """Inversed Hessian matrix.

    In mathematics, the Hessian matrix or Hessian is a square matrix of
    second-order partial derivatives of a scalar-valued function, or scalar
    field. It describes the local curvature of a function of many variables.

    Attributes
    ----------
        - with_labels

    Internal Attributes
    -------------------
        - label, matrix, correlation_matrix

    Methods
    -------
        - get_variable_names
        - get_variable_by_name, set_variable_by_name
    """

    ATTR_MANDATORY_NAMES = ("with_labels", )
    ATTR_MANDATORY_TYPES = (str, )
    # ("matrix", "matrix", "matrix", "matrix")
    ATTR_MANDATORY_CIF = ("with_labels", )

    ATTR_OPTIONAL_NAMES = ()
    ATTR_OPTIONAL_TYPES = ()
    ATTR_OPTIONAL_CIF = ()

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ("label", "matrix", "correlation_matrix", "sigma")
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

    PREFIX = "inversed_hessian"

    def __init__(self, **kwargs) -> NoReturn:
        super(InversedHessian, self).__init__()

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
        l_1 = (self.with_labels).strip().split("\n")

        l_label, l_matrix = [], []
        for line in l_1[0:]:
            l_1 = line.strip().split()
            l_label.append(str(l_1[0]))
            l_matrix.append(l_1[1:])

        np_label = numpy.array(l_label, dtype=str)
        np_matrix = numpy.array(l_matrix, dtype=float)
        np_correlation_matrix, np_sigma = hessian_to_correlation(np_matrix)
        
        self.set_labels(np_label)
        self.set_inversed_hessian(np_matrix)
        self.set_correlation_matrix(np_correlation_matrix)
        self.set_sigmas(np_sigma)

    def set_labels(self, labels) -> NoReturn:
        """Form 2theta_phi_intensity from internal attributes."""
        np_label = numpy.array(labels, dtype=str)
        self.__dict__["label"] = np_label

    def set_inversed_hessian(self, matrix) -> NoReturn:
        """Form 2theta_phi_intensity from internal attributes."""
        self.__dict__["matrix"] = matrix

    def set_correlation_matrix(self, matrix) -> NoReturn:
        """Form 2theta_phi_intensity from internal attributes."""
        self.__dict__["correlation_matrix"] = matrix

    def set_sigmas(self, sigmas) -> NoReturn:
        """Form 2theta_phi_intensity from internal attributes."""
        self.__dict__["sigma"] = sigmas

    def form_inversed_hessian(self):
        f_label = self.is_attribute("label")
        f_hessian = self.is_attribute("matrix")
        f_correlation = self.is_attribute("correlation_matrix")
        f_sigma = self.is_attribute("sigma")
        if f_label & f_hessian:
            ls_out = [f"{lab:7} "+" ".join([f"{val:.10f}" for val in hess])
                      for lab, hess in zip(self.label, self.matrix)]
            self.with_labels = "\n".join(ls_out)
        elif f_label & f_correlation & f_sigma:
            np_correlation = self.correlation_matrix
            np_sigma = self.sigma
            np_na = numpy.newaxis
            np_hessian = np_correlation*np_sigma[:, np_na]*np_sigma[np_na, :]
            self.set_inversed_hessian(np_hessian)
            ls_out = [f"{lab:7} "+" ".join([f"{val:.10f}" for val in hess])
                      for lab, hess in zip(self.label, self.matrix)]
            self.with_labels = "\n".join(ls_out)

    def report(self):
        if self.is_defined():
            self.form_object()
            np_correlation = self.correlation_matrix
            np_sigma = self.sigma
            np_label = self.label
            ls_out = ["# Sigmas and correlation matrix"]
            ls_out.append("|SIGMAS | Cor.Matrix:|"+"|".join(np_label) + "|")
            ls_out.append("|-----|-----|"+"|".join(["-----" for hh in
                                                    np_label]) + "|")
            
            for i1, lab, mat, sig in zip(range(len(np_label)), np_label,
                                         np_correlation, np_sigma):
                ls_out.append(f"|{sig:.5f}|{lab:}|"+"|".join([
                    f"{val:5.2f}" for i2, val in enumerate(mat)]) + "|")
            # ls_out.extend([f"{i_n+1:2} -- {lab:7}: {sig:.5f}" for i_n, lab, sig in zip(
            #     range(len(np_sigma)), np_label, np_sigma)])
            # ls_out.append("\nCorrelation matrix:")
            # ls_out.extend([f" {i_n+1:2}: "+" ".join([
            #     f"{val:5.2f}" for val in mat]) for i_n, mat in enumerate(np_correlation)])
            return "\n".join(ls_out)

    def report_html(self):
        ls_html = []
        if self.is_defined():
            self.form_object()
            np_correlation = self.correlation_matrix
            np_sigma = self.sigma
            np_label = self.label
            ls_html.append("<table>")
            ls_html.append("<tr><th>SIGMAS</th><th>CORRELATION MATRIX:</th>"+"".join(
                [f"<th>{lab:}</th>" for lab in np_label]) + "</tr>")
            for i1, lab, mat, sig in zip(range(len(np_label)), np_label, np_correlation, np_sigma):
                ls_html.append(f"<tr><td>{sig:.5f}</td><th>{lab:}</th>"+"".join([
                    f"<td bgcolor='#ff2200'>{val:5.2f}</td>" if ((abs(val)>0.7) & (i1 != i2 ))
                    else f"<td>{val:5.2f}</td>" for i2, val in enumerate(mat)]) + "</tr>")
            ls_html.append("</table>")
        return "".join(ls_html)
            

# s_cont = """
#   _hessian_matrix_with_labels
#   ;
#   length_a 8 3 0
#   length_b 3 9 2
#   length_c 2 8 3
#   ;

# """

# obj = InversedHessian.from_cif(s_cont)
# obj.form_object()
# print(obj, end="\n\n")
# print(obj.matrix, end="\n\n")
# print(obj.label, end="\n\n")
# print(obj.correlation_matrix, end="\n\n")
# print(obj.sigma, end="\n\n")
# print(obj.with_labels, end="\n\n")
# obj.form_inversed_hessian()

# print(obj, end="\n\n")
