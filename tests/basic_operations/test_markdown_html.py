"""Markdown of the reports into html for the View tab of CrysPy Editor."""
import numpy

from cryspy.A_functions_base.function_1_markdown import md_to_html
from cryspy.A_functions_base.function_1_agreement_factors import \
    calc_agreement_factors, report_powder_agreement


def test_md_to_html_closes_all_tags():
    s_md = "# Title\n| a | 1 |\n| b | 2 |\n## Next\n- item 1\n- item 2\ntext"
    s_html = md_to_html(s_md)
    for tag in ("table", "tr", "td", "ul", "li", "p", "h3", "h4"):
        assert s_html.count(f"<{tag:}>") + s_html.count(f"<{tag:} ") == \
            s_html.count(f"</{tag:}>"), tag
    # the heading after the table is not inside the table
    assert s_html.index("</table>") < s_html.index("<h4>")


def test_md_to_html_table_header_only_with_separator():
    s_html = md_to_html("| name | value |\n|:---|---:|\n| a | 1.5 |")
    assert "<th align='left'>name</th>" in s_html
    assert "<td align='right'>1.5</td>" in s_html
    s_html = md_to_html("| Number of reflections | 95 |\n| Range of h | -1 … 2 |")
    assert "<th" not in s_html
    assert "<td align='left'>Number of reflections</td>" in s_html
    assert "<td align='right'>95</td>" in s_html
    assert "border='1'" in s_html


def test_md_to_html_escape_and_stars():
    s_html = md_to_html("(N-1)*sigma < abs(x) < N*sigma and **bold** *it*")
    assert "(N-1)*sigma &lt; abs(x) &lt; N*sigma" in s_html
    assert "<b>bold</b>" in s_html
    assert "<i>it</i>" in s_html


def test_md_to_html_rows_of_different_length():
    s_html = md_to_html("| a | b | c |\n|---|---|---|\n| 1 |")
    assert s_html.count("<td") == 3


def test_md_to_html_empty():
    assert md_to_html("  \n ") == ""


def test_calc_agreement_factors():
    y_exp = numpy.array([10., 20., 30., numpy.nan], dtype=float)
    y_sigma = numpy.array([1., 2., 3., 1.], dtype=float)
    y_calc = numpy.array([11., 20., 27., 5.], dtype=float)
    res = calc_agreement_factors(y_exp, y_sigma, y_calc)
    assert res["n_point"] == 3
    assert numpy.isclose(res["chi_sq"], 2.)
    assert numpy.isclose(res["chi_sq_per_n"], 2./3.)
    assert numpy.isclose(res["r_p"], 4./60.)
    assert numpy.isclose(res["r_wp"], numpy.sqrt(2./300.))
    res = calc_agreement_factors(
        y_exp, y_sigma, y_calc, flag_in=numpy.array([False, True, True, True]))
    assert res["n_point"] == 2


def test_report_powder_agreement_polarized():
    n = 5
    dict_signal = {
        "intensity_plus_net": numpy.full(n, 4.), "intensity_minus_net": numpy.full(n, 2.),
        "intensity_bkg_calc": numpy.full(n, 1.),
        "intensity_plus": numpy.full(n, 5.), "intensity_plus_sigma": numpy.ones(n),
        "intensity_minus": numpy.full(n, 2.), "intensity_minus_sigma": numpy.ones(n),
        "excluded": numpy.array([True, False, False, False, False])}
    s_md = report_powder_agreement(dict_signal, flag_sum=True, flag_diff=False)
    assert "| Excluded points | 1 |" in s_md
    assert "| up + down | yes | 4 |" in s_md
    assert "| up − down | no | 5 |" in s_md
    s_html = md_to_html(s_md)
    assert s_html.count("<table") == 2
