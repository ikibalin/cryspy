"""Unknown CIF categories are written in the order of the file."""
import cryspy

RCIF_UNKNOWN = """data_a

_zeta_one 1
_alpha_two 2
_mu_three 3

loop_
_omega_x
_omega_y
1 2

loop_
_beta_x
_beta_y
3 4
"""


def test_unknown_categories_in_order_of_file():
    text = cryspy.str_to_globaln(RCIF_UNKNOWN).to_cif()
    l_position = [text.find(name) for name in
                  ("_zeta_one", "_alpha_two", "_mu_three", "_omega_x", "_beta_x")]
    assert -1 not in l_position
    assert l_position == sorted(l_position)
