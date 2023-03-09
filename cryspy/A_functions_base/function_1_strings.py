"""
Functions to work with strings.
"""

from typing import Tuple, Union
import numpy
from fractions import Fraction


def ciftext_to_html(text: str) -> str:
    """Convert cif text to html format."""
    l_text_new = []
    for line in text.split("\n"):
        if line.strip().startswith("###"):
            line_new = f"<h3>{line:}</h3>"
        elif line.strip().startswith("##"):
            line_new = f"<h2>{line:}</h2>"
        elif line.strip().startswith("#"):
            line_new = f"<h1>{line:}</h1>"
        else:
            line_new = line
        l_text_new.append(line_new)
    text_html = "<br>".join(l_text_new)
    # text_html = "<b>"+text.replace("\n", "<br>")+"<\b>"
    return text_html


def find_prefix(*argv) -> str:
    """Find prefix
    
    Example 1 (trivial)
    -------------------
    
    loop_
    _atom_site_label
    _atom_site_fract_x
    _atom_site_fract_y
    
    Output: "_atom_site"

    Example 2 (not trivial)
    -----------------------
    loop_
        _geom_torsion
        _geom_torsion_atom_site_label_1
        _geom_torsion_atom_site_label_2
    
    Output: "_geom_torsion"
    (name for item "geom_torsion" is label)
    
    """
    c_string = common_string(*argv)
    ind_point = c_string.rfind(".")
    if ind_point != -1:
        prefix = c_string[:ind_point]
    elif ((len(argv) > 1) & (c_string in argv)):
        prefix = c_string
    else:
        ind_und = c_string.rfind("_")
        prefix = c_string[:ind_und]
    return prefix


def common_string(*argv) -> str:
    """Return common begining of input strings."""
    ls_out = []
    for args in zip(*argv):
        arg_1 = args[0]
        flag = all([arg == arg_1 for arg in args])
        if flag:
            ls_out.append(arg_1)
        else:
            break
    return "".join(ls_out)



def string_to_value_error_mark(string: str) -> Tuple[float, Union[float, None], str]:
    """
    Convert string to float and error.

    Parameters
    ----------
    string : str
        DESCRIPTION.

    Returns
    -------
    value : float
        Value.
    error : float
        Error.

    """
    value, error, mark = None, None, ""
    ind_1 = string.find("(")
    s_sigma = ""
    if string == ".":
        pass 
    elif ind_1 != -1:
        ind_2 = string.find(")")
        if ind_2 > ind_1:
            s_sigma = string[(ind_1+1):ind_2]
            if not(s_sigma.isdigit()):
                s_sigma = ""
        str_1 = string.split("(")[0]
        value = float(str_1)
        mark = string[ind_2+1:].strip()
        if s_sigma != "":
            s_h = "".join(["0" if _h.isdigit() else _h for _h in
                           str_1[:-len(s_sigma)]])
            error = abs(float(s_h+s_sigma))
        else:
            error = 0.
    else:
        try:
            value = float(string)
        except ValueError:
            value = None
    return value, error, mark

def value_error_mark_to_string(value: float, error: float, mark: str) -> str:
    """
    Convert value and error to string

    Parameters
    ----------
    value : float
        DESCRIPTION.
    error : float
        DESCRIPTION.

    Returns
    -------
    str
        DESCRIPTION.

    """
    if ((value is None) | (value is numpy.nan)):
        string = "."
    if not((error == 0.) | (error is None)):
        val_hh = numpy.log10(error)
        n_power = int(val_hh)
        if val_hh <= 0:
            val_1 = numpy.round(value, decimals = -1*int(n_power)+2)
            val_2 = numpy.round(error, decimals = -1*int(n_power)+2)
            i_val_11 = int(val_1)
            s_sign = ""
            if val_1 < 0.:
                s_sign = "-"
            s_val_12 = ("{:}".format(int(abs(val_1) % 1.*10**(-n_power+2)))
                        ).rjust(int(-n_power+2), "0")
            val_2 = int(val_2*10**(-n_power+2))
            string = f"{s_sign:}{abs(i_val_11):}.{s_val_12:}({int(val_2):}){mark.strip():}"
        else:
            val_1 = numpy.round(value)
            val_2 = numpy.round(error)
            string = "{:}({:}){:}".format(int(val_1), int(val_2),mark.strip())
    elif (error == 0.):
        string = "{:}(){:}".format(value, mark.strip())
    else:
        string = "{:}".format(value)
    return string

##FIXME to del
#def string_to_value_error(string: str) -> Tuple[float, Union[float, None]]:
#    """
#    Convert string to float and error.
#
#    Parameters
#    ----------
#    string : str
#        DESCRIPTION.
#
#    Returns
#    -------
#    value : float
#        Value.
#    error : float
#        Error.
#
#    """
#    value, error = None, None
#    ind_1 = string.find("(")
#    s_sigma = ""
#    if value == ".":
#       pass 
#    elif ind_1 != -1:
#        ind_2 = string.find(")")
#        if ind_2 > ind_1:
#            s_sigma = string[(ind_1+1):ind_2]
#            if not(s_sigma.isdigit()):
#                s_sigma = ""
#        str_1 = string.split("(")[0]
#        value = float(str_1)
#        if s_sigma != "":
#            s_h = "".join(["0" if _h.isdigit() else _h for _h in
#                           str_1[:-len(s_sigma)]])
#            error = abs(float(s_h+s_sigma))
#        else:
#            error = 0.
#    else:
#        try:
#            value = float(string)
#        except ValueError:
#            value = None
#    return value, error
#
# #FIXME to del
# def value_error_to_string(value: float, error: float) -> str:
#     """
#     Convert value and error to string
# 
#     Parameters
#     ----------
#     value : float
#         DESCRIPTION.
#     error : float
#         DESCRIPTION.
# 
#     Returns
#     -------
#     str
#         DESCRIPTION.
# 
#     """
#     if ((value is None) | (value is numpy.nan)):
#         string = "."
#     if not((error == 0.) | (error is None)):
#         val_hh = numpy.log10(error)
#         n_power = int(val_hh)
#         if val_hh <= 0:
#             val_1 = numpy.round(value, decimals = -1*int(n_power)+2)
#             val_2 = numpy.round(error, decimals = -1*int(n_power)+2)
#             i_val_11 = int(val_1)
#             s_sign = ""
#             if val_1 < 0.:
#                 s_sign = "-"
#             s_val_12 = ("{:}".format(int(abs(val_1) % 1.*10**(-n_power+2)))
#                         ).rjust(int(-n_power+2), "0")
#             val_2 = int(val_2*10**(-n_power+2))
#             string = f"{s_sign:}{abs(i_val_11):}.{s_val_12:}({int(val_2):})"
#         else:
#             val_1 = numpy.round(value)
#             val_2 = numpy.round(error)
#             string = "{:}({:})".format(int(val_1), int(val_2))
#     elif (error == 0.):
#         string = "{:}()".format(value)
#     else:
#         string = "{:}".format(value)
#     return string
# 

def transform_string_to_r_b(name: str, labels=("x", "y", "z")) -> Tuple:
    """
    transform string to rotation part and offset: 
    x,y,-z   ->   0.0 1 0 0  0.0 0 1 0  0.0 0 0 -1
    """
    l_name = "".join(name.strip().split()).lstrip("(").rstrip(")").split(",")
    rij, bi = [], []
    for _name in l_name:
        coefficients, offset = transform_string_to_digits(_name, labels)
        rij.append(coefficients)
        bi.append(offset)
    res_rij = numpy.array(rij, dtype=object)
    res_bi = numpy.array(bi, dtype=object)
    return res_rij, res_bi



def transform_string_to_digits(name: str, labels: Tuple[str]):
    """
    Multiplication has to be implicit, division must be explicit.
    White space within the string is optional 

    Example: 
    transform_string_to_digits('-a+2x/7-0.7t+4', ('a', 'x', 't')) = 
    = (Fraction(-1,1),fraction(2,7),Fraction(-7,10)), Fraction(4,1)
    """
    coefficients = []
    offset = Fraction(0, 1).limit_denominator(10)
    l_name = name.strip().replace(" ", "").replace("+", " +").replace("-", " -").split()
    for _label in labels:
        res = Fraction(0, 1).limit_denominator(10)
        for _name in l_name:
            if _label in _name:
                s_1 = _name.replace(_label, "").replace("+/", "+1/").replace("-/", "-1/")
                if s_1 == "": s_1 = "1"
                if s_1.startswith("/"): s_1 = "1" + s_1
                if s_1.endswith("+"): s_1 = s_1 + "1"
                if s_1.endswith("-"): s_1 = s_1 + "1"
                res += Fraction(s_1).limit_denominator(10)
        coefficients.append(res)
    res = Fraction(0, 1).limit_denominator(10)
    for _name in l_name:
        flag = all([not (_label in _name) for _label in labels])
        if flag:
            res += Fraction(_name).limit_denominator(10)
    offset = res
    return coefficients, offset


def transform_fraction_with_label_to_string(number: Fraction, label: str) -> str:
    if isinstance(number, Fraction):
        val = number
    else:
        val = Fraction(number).limit_denominator(10)
    if val == Fraction(0, 1):
        res = ""
    elif val == Fraction(1, 1):
        res = f"{label:}"
    elif val == Fraction(-1, 1):
        res = f"-{label:}"
    elif val.denominator == 1:
        res = f"{val.numerator}{label:}"
    elif val.numerator == 1:
        res = f"{label:}/{val.denominator}"
    else:
        res = f"{val.numerator}{label:}/{val.denominator}"
    return res


def transform_digits_to_string(labels: Tuple[str], coefficients, offset: Fraction) -> str:
    """
Form a string from digits.

Keyword arguments: 

    labels: the tuple of lablels (ex.: ('x', 'y', 'z') or ('a', 'b', 'c')))
    coefficients: the parameters in front of label (ex.: (1.0, 0.5, 0.0))
    offset: the number (ex.: 2/3)

Output arguments is string.

Example:
>>> transform_digits_to_string(('x', 'y', 'z'), (1.0, 0.5, 0.0), 0.6666667)
x+1/2y+2/3
    """
    l_res = []
    for _coefficient, _label in zip(coefficients, labels):
        _name = transform_fraction_with_label_to_string(_coefficient, _label)
        if _name == "":
            pass
        elif _name.startswith("-"):
            l_res.append(_name)
        elif l_res == []:
            l_res.append(_name)
        else:
            l_res.append(f"+{_name:}")
    _name = str(Fraction(offset).limit_denominator(10))
    if _name == "0":
        if l_res == []:
            l_res.append(_name)
    elif ((l_res == []) | (_name.startswith("-"))):
        l_res.append(_name)
    else:
        l_res.append(f"+{_name:}")
    return "".join(l_res)


def transform_r_b_to_string(r, b, labels=("x", "y", "z")) -> str:
    l_res = [transform_digits_to_string(labels, _ri, _bi) for _ri, _bi in zip(r, b)]
    return ",".join(l_res)
