"""Reading the multi-pattern format of FullProf (NPATT): the patterns, the
phases that contribute to each of them and the refusals."""
import warnings

import pytest

from cryspy.H_functions_global.fullprof import FullProfFormatError
from cryspy.H_functions_global.fullprof_multi_pattern import (
    is_multi_pattern, read_fullprof_multi_pcr)

# two patterns (neutrons, X-rays), two phases: 'a' in both patterns, 'b' in the second only
PHASE_PROFILE = """  {scale:}  0.0 0.0 0.0 0.0 0.0 0
  {code:}  0.0 0.0 0.0 0.0 0.0
  0.1 -0.2 0.3 0.0 0.05 0.0 0.0 0
  0.0 0.0 0.0 0.0 0.0 0.0 0.0
  {cell:}
  31.0 31.0 31.0 0.0 0.0 0.0
  0.0 0.0 0.1 0.0 0.0 0.0 0.0 0.0
  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
"""
TEXT_PCR = """COMM two patterns
NPATT 2 1 1 <- Flags
W_PAT 0.5 0.5
!Nph Dum Ias Nre Cry Opt Aut
   2   0   0   0   0   0   1
!Job Npr Nba Nex Nsc Nor Iwg Ilo Res Ste Uni Cor Anm Int
   1   7   0   1   0   0   0   0   0   0   0   0   0   0
   0   7   2   0   0   0   0   0   0   0   0   0   0   0
neutrons.dat
xrays.dat
   0   1   0   0   0   0
   0   0   1   0   4   0   1   6   1   1   0
   0   0   1   0   0   0   1  10   1   1   0
 1.9 1.9 1.0 70.0 6.0 0.0 0.0 160.0 0.0 0.0
 1.54 1.54 0.0 50.0 6.0 0.8 0.0 160.0 0.0 0.0
  3  0.10  1.00  1.00  1.00  1.00
    10.0 0.05 155.0 0.0 0.0
    12.0 0.02 140.0 0.0 0.0
  0.0 10.0
  12.0 100.0 0.0
  140.0 90.0 0.0
  5 !Number of refined parameters
 -0.08 81.0 0.0 0.0 0.0 0.0 0.0 0.0 0
 200.0 40.0 0.0 0.0 0.0 0.0
 151.0 0.0 0.0 0.0 0.0 0.0
 0.02 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0
a
   1   0   0   0   0   0   0  0.0   0   1
   0   3   0   0   0   0      0
 1 1
   0   7   1   0   0
  0.0 0.0 1.0 1.0 0.0 0.0 0.0
   0   7   0   0   0
  0.0 0.0 1.0 1.0 0.0 0.0 0.0
  3.5 0.0 BVS
  2 1 0.0
  PB+2 S+6
  O-2
F m -3 m
Ce CE 0.0 0.0 0.0 0.5 0.0208 0 0 0 0
  0.0 0.0 0.0 0.0 0.0
""" + PHASE_PROFILE.format(scale=1.5, code=11.0, cell="5.41 5.41 5.41 90 90 90") \
    + PHASE_PROFILE.format(scale=0.002, code=21.0, cell="5.41 5.41 5.41 90 90 90") + """b
   1   0   0   0   0   0   0  0.0   0   0
 0 1
   0   7   0   0   0
  0.0 0.0 1.0 1.0 0.0 0.0 0.0
P 1
Si SI 0.0 0.0 0.0 0.5 1.0 0 0 0 0
  0.0 0.0 0.0 0.0 0.0
""" + PHASE_PROFILE.format(scale=0.003, code=0.0, cell="3.0 4.0 5.0 90 90 90")


def write_pcr(tmp_path, text: str) -> str:
    file_pcr = tmp_path / "two.pcr"
    file_pcr.write_text(text)
    return str(file_pcr)


def test_is_multi_pattern(tmp_path):
    assert is_multi_pattern(write_pcr(tmp_path, TEXT_PCR))
    # the title is optional
    assert is_multi_pattern(write_pcr(tmp_path, TEXT_PCR.split("\n", 1)[1]))
    assert not is_multi_pattern(write_pcr(tmp_path, "COMM one pattern\n 1 7 1 0 0\n"))


def test_patterns_and_contributions(tmp_path):
    title, l_pattern, l_crystal = read_fullprof_multi_pcr(write_pcr(tmp_path, TEXT_PCR))
    assert title == "two patterns"
    assert [pattern.file_data for pattern in l_pattern] == ["neutrons.dat", "xrays.dat"]
    assert [pattern.job for pattern in l_pattern] == [1, 0]
    assert [pattern.wavelength for pattern in l_pattern] == [1.9, 1.54]
    assert l_pattern[0].excluded == [[0.0, 10.0]]
    assert len(l_pattern[1].background_points) == 2
    assert l_pattern[0].zero == pytest.approx(-0.08)
    assert l_pattern[0].background_polynomial[0] == pytest.approx(200.)
    # 'a' contributes to both patterns with its own scale, 'b' to the second one
    assert [[phase.label for phase in pattern.phases] for pattern in l_pattern] == [
        ["a"], ["a", "b"]]
    assert [phase.scale for phase in l_pattern[1].phases] == [0.002, 0.003]
    assert l_pattern[0].phases[0].scale_code == 11.
    # the atoms are common, the crystals are the phases of the first pattern
    assert [crystal.label for crystal in l_crystal] == ["a", "b"]
    assert l_crystal[0].atoms is l_pattern[1].phases[0].atoms
    assert l_crystal[1].cell == [3.0, 4.0, 5.0, 90., 90., 90.]


def test_weights_are_warned(tmp_path):
    text = TEXT_PCR.replace("W_PAT 0.5 0.5", "W_PAT 0.8 0.2")
    with warnings.catch_warnings(record=True) as l_warning:
        warnings.simplefilter("always")
        read_fullprof_multi_pcr(write_pcr(tmp_path, text))
    assert any("W_PAT" in str(warning.message) for warning in l_warning)


@pytest.mark.parametrize("old, new, message", [
    # integrated intensities (a single crystal) in the second pattern
    ("   0   7   2   0   0   0   0   0   0   0   0   0   0   0",
     "   0   7   2   0   0   0   0   0   0   0   0   0   0   1", "integrated intensities"),
    # TOF
    ("   1   7   0   1   0   0   0   0   0   0   0   0   0   0",
     "  -1   7   0   1   0   0   0   0   0   0   0   0   0   0", "Job = -1"),
    # the profile of phase 'b' in the second pattern
    ("   0   7   0   0   0\n  0.0 0.0 1.0 1.0 0.0 0.0 0.0\nP 1",
     "   0   5   0   0   0\n  0.0 0.0 1.0 1.0 0.0 0.0 0.0\nP 1", "Npr = 5")])
def test_refused_options(tmp_path, old, new, message):
    assert old in TEXT_PCR
    with pytest.raises(FullProfFormatError, match=message):
        read_fullprof_multi_pcr(write_pcr(tmp_path, TEXT_PCR.replace(old, new)))
