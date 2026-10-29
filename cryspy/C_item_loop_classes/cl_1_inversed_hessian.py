"""InversedHessian class."""
import html
import numpy
from typing import NoReturn

from cryspy.A_functions_base.function_1_markdown import HTML_TABLE_OPEN
from cryspy.B_parent_classes.cl_1_item import ItemN

def inversed_hessian_to_correlation(inv_hessian):
    """Calculate correlation matrix and sigmas for inversed hessian matrix."""
    np_sigma_sq = numpy.diagonal(inv_hessian)
    np_sigma = numpy.sqrt(numpy.abs(np_sigma_sq))
    # np_sigma_sq_matrix = np_sigma[:, np_na] * np_sigma[np_na, :]
    np_sigma_sq_matrix = numpy.outer(np_sigma, np_sigma)
    eps = 1e-12 
    np_safe_sigma = numpy.where(np_sigma_sq_matrix == 0, eps, np_sigma_sq_matrix)
    np_correlation_matrix = inv_hessian / np_safe_sigma
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

        # the order of the parameters is kept: it is the order of the log
        # of the refinement
        np_label = numpy.array(l_label, dtype=str)
        np_matrix = numpy.array(l_matrix, dtype=float)
        np_correlation_matrix, np_sigma = inversed_hessian_to_correlation(np_matrix)
        
        self.set_labels(np_label)
        self.set_inversed_hessian(np_matrix)
        self.set_correlation_matrix(np_correlation_matrix)
        self.set_sigmas(np_sigma)
        # self.form_inversed_hessian()

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
            ls_out = [f"{lab:7} "+" ".join([f"{val:.10e}" for val in hess])
                      for lab, hess in zip(self.label, self.matrix)]
            self.with_labels = "\n".join(ls_out)

        elif f_label & f_correlation & f_sigma:
            np_correlation = self.correlation_matrix
            np_sigma = self.sigma
            
            np_na = numpy.newaxis
            np_hessian = np_correlation*np_sigma[:, np_na]*np_sigma[np_na, :]
            self.set_inversed_hessian(np_hessian)
            ls_out = [f"{lab:7} "+" ".join([f"{val:.10e}" for val in hess])
                      for lab, hess in zip(self.label, self.matrix)]
            self.with_labels = "\n".join(ls_out)

    def report(self):
        """Sigmas and correlation matrix in markdown."""
        if not self.is_defined():
            return ""
        self.form_object()
        np_correlation = self.correlation_matrix
        np_sigma = self.sigma
        np_label = self.label
        l_number = [str(i_n+1) for i_n in range(np_label.size)]
        ls_out = ["# Sigmas and correlation matrix", ""]
        ls_out.append("| # | Parameter | Sigma |")
        ls_out.append("|---|---|---|")
        for numb, lab, sig in zip(l_number, np_label, np_sigma):
            ls_out.append(f"| {numb:} | {lab:} | {sig:.5g} |")
        ls_out.append("")
        ls_out.append("**Correlation matrix:**")
        ls_out.append("")
        ls_out.append("| # | " + " | ".join(l_number) + " |")
        ls_out.append("|---|" + "|".join(["---" for hh in l_number]) + "|")
        for numb, mat in zip(l_number, np_correlation):
            ls_out.append(f"| {numb:} | " + " | ".join(
                [f"{val:.2f}" for val in mat]) + " |")
        return "\n".join(ls_out)

    def report_html(self):
        """Sigmas and correlation matrix in html (Qt rich text).

        Strong correlations (|r| > 0.7, off the diagonal) are on a red
        background with white text, readable in light and dark themes.
        """
        if not self.is_defined():
            return ""
        self.form_object()
        np_correlation = self.correlation_matrix
        np_sigma = self.sigma
        np_label = self.label
        l_number = [str(i_n+1) for i_n in range(np_label.size)]
        s_table = HTML_TABLE_OPEN

        ls_html = ["<h3>Sigmas and correlation matrix</h3>", s_table]
        ls_html.append("<tr><th>#</th><th align='left'>Parameter</th><th>Sigma</th></tr>")
        for numb, lab, sig in zip(l_number, np_label, np_sigma):
            ls_html.append(
                f"<tr><td align='right'>{numb:}</td><td>{html.escape(str(lab)):}</td>"
                f"<td align='right'>{sig:.5g}</td></tr>")
        ls_html.append("</table>")

        ls_html.append("<p><b>Correlation matrix</b> (|r| &gt; 0.7 is marked):</p>")
        ls_html.append(s_table)
        ls_html.append("<tr><th>#</th>" + "".join(
            [f"<th>{numb:}</th>" for numb in l_number]) + "</tr>")
        for i1, numb, mat in zip(range(len(l_number)), l_number, np_correlation):
            l_td = []
            for i2, val in enumerate(mat):
                if i1 == i2:
                    l_td.append(f"<td align='right'>{val:.2f}</td>")
                elif abs(val) > 0.7:
                    l_td.append(
                        "<td align='right' bgcolor='#c62828'>"
                        f"<font color='#ffffff'><b>{val:.2f}</b></font></td>")
                else:
                    l_td.append(f"<td align='right'>{val:.2f}</td>")
            ls_html.append(f"<tr><th>{numb:}</th>" + "".join(l_td) + "</tr>")
        ls_html.append("</table>")
        return "".join(ls_html)
