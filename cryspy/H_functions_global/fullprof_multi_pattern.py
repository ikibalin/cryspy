"""Multi-pattern .pcr files of FullProf (the line 'NPATT n'): every
pattern is one experiment of CrysPy, the crystals are common.

FullProf refines several patterns together with the common structural
parameters (atoms, cell); the profile parameters (scale, U V W ...) are
given for each phase and each pattern. CrysPy does the same with several
experiments (pd) that name the same crystals in their phase loops.
"""
import copy
import logging
import os
from typing import List, Tuple

import numpy

from cryspy.B_parent_classes.cl_4_global import GlobalN
from cryspy.H_functions_global.function_1_cryspy_objects import str_to_globaln
from cryspy.H_functions_global.fullprof import (
    JOBS_SIMULATION, CodeCounter, FullProfFormatError, PatternPcr, form_atom_lines,
    form_phase_head, take_pd_data, warn_special_positions, write_xydata, FullProfPattern, FullProfPhase, PcrLines,
    calc_occupancies, check_data_format, check_pattern_options, check_phase_options,
    form_crystal_rcif, form_pd_rcif, read_atom_lines, read_background_points_and_excluded,
    read_distance_lines, read_fullprof_data, read_more_options, read_powder_profile,
    read_wavelength_line, read_zero_and_background, take_background, take_label,
    take_space_group_symbol, warn_conversion, warn_shared_codes)

logger = logging.getLogger(__name__)


def is_multi_pattern(file_pcr: str) -> bool:
    """Whether a .pcr file has the multi-pattern format (a line 'NPATT n')."""
    with open(file_pcr, "r", errors="replace") as fid:
        lines = PcrLines(fid.read())
    # the title line 'COMM ...' is optional
    return any(line.upper().startswith("NPATT") for line in lines.lines[:2])


def read_pattern_flags(lines: PcrLines) -> Tuple[List[int], List[float]]:
    """'NPATT n flag_1 ... flag_n' and 'W_PAT w_1 ... w_n'."""
    l_word = lines.next_text().split()
    n_pattern = int(l_word[1])
    l_flag = [int(word) for word in l_word[2:2 + n_pattern]]
    l_flag += [1] * (n_pattern - len(l_flag))
    l_word = lines.next_text().split()
    if l_word[0].upper() != "W_PAT":
        raise FullProfFormatError("The line 'W_PAT' is expected after the line 'NPATT'.")
    l_weight = [float(word) for word in l_word[1:1 + n_pattern]]
    return l_flag, l_weight


def warn_flags_and_weights(l_flag: List[int], l_weight: List[float]):
    """Warn about excluded patterns and unequal weights (not converted)."""
    if any(flag == 0 for flag in l_flag):
        warn_conversion(
            "Some patterns are excluded from the refinement in FullProf (flag 0); "
            "CrysPy refines all experiments.")
    if len(l_weight) > 1 and not numpy.allclose(l_weight, l_weight[0]):
        warn_conversion(
            f"The weights of the patterns W_PAT = {l_weight:} are not converted: CrysPy "
            "adds the chi-squares of the experiments with equal weights.")


def read_job_lines(lines: PcrLines, n_pattern: int) -> Tuple[List[FullProfPattern], List[int]]:
    """The lines 'Job Npr Nba Nex Nsc Nor Iwg Ilo Res Ste Uni Cor Anm Int',
    one for each pattern: the patterns and their numbers of excluded regions."""
    l_pattern, l_nex = [], []
    for index in range(n_pattern):
        l_job = lines.next_integers(12)
        job, _, nba, nex, nsc, nor = l_job[:6]
        res = l_job[8]
        if len(l_job) > 13 and l_job[13] != 0:
            raise FullProfFormatError(
                f"Pattern {index + 1:} of the multi-pattern file has integrated intensities "
                "(Int = 1, a single crystal or a powder): only powder patterns are "
                "converted from multi-pattern files.")
        check_pattern_options(job, nba, nsc, res)
        if nba == -1:
            raise FullProfFormatError(
                f"Pattern {index + 1:}: the background Nba = -1 is not converted from "
                "multi-pattern files.")
        pattern = FullProfPattern()
        pattern.job, pattern.nba, pattern.nor = job, nba, nor
        l_pattern.append(pattern)
        l_nex.append(nex)
    return l_pattern, l_nex


def read_multi_phase(lines: PcrLines, l_pattern: List[FullProfPattern], l_used: List[str]
                     ) -> FullProfPhase:
    """A phase of a multi-pattern file; its copy with the profile parameters
    of each pattern is added to the phases of the pattern."""
    phase = FullProfPhase()
    phase.label = take_label(lines.next_text(), l_used)
    l_value = lines.next_numbers(10)
    nat, dis, ang, jbt, isy, i_str, furth = [int(value) for value in l_value[:7]]
    nvk, more = int(l_value[8]), int(l_value[9])
    d_option = {"Isy": isy, "Str": i_str, "Furth": furth, "Nvk": nvk, "Dis": dis, "Ang": ang}
    check_phase_options(phase.label, jbt, 7, d_option)
    jdi = read_more_options(lines, phase.label) if more == 1 else 0
    l_contribution = lines.next_integers(len(l_pattern))[:len(l_pattern)]
    l_index = [index for index, flag in enumerate(l_contribution) if flag != 0]
    for _ in l_index:
        # Irf Npr Jtyp [Nsp_Ref Ph_Shift] and Pr1 Pr2 Pr3 ... for each pattern
        l_irf = lines.next_integers(3)
        check_phase_options(phase.label, 0, l_irf[1], {"Irf": l_irf[0]})
        phase.pref_axis = lines.next_numbers(3)[:3]
    if more == 1:
        read_distance_lines(lines, jdi)
    phase.space_group = take_space_group_symbol(lines.next_text().split("<--")[0].strip())
    phase.atoms = read_atom_lines(lines, nat)
    for index in l_index:
        phase_pattern = copy.copy(phase)
        read_powder_profile(lines, phase_pattern, l_pattern[index].nor)
        l_pattern[index].phases.append(phase_pattern)
    return phase


def read_fullprof_multi_pcr(file_pcr: str) -> Tuple[str, List[FullProfPattern], List[FullProfPhase]]:
    """Read a multi-pattern .pcr file of FullProf: the title, the patterns
    (with their phases and profile parameters) and the crystals (the phases
    with the cell of the first pattern they contribute to)."""
    with open(file_pcr, "r", errors="replace") as fid:
        lines = PcrLines(fid.read())
    title = ""
    if not lines.lines[0].upper().startswith("NPATT"):
        # the title, with or without 'COMM'
        title = lines.next_text()
        title = title[4:].strip() if title.upper().startswith("COMM") else title.strip()
    l_flag, l_weight = read_pattern_flags(lines)
    warn_flags_and_weights(l_flag, l_weight)
    l_global = lines.next_integers(5)
    nph, nre, cry = l_global[0], l_global[3], l_global[4]
    if cry != 0:
        raise FullProfFormatError("Single crystal data (Cry != 0) of FullProf are not converted.")
    l_pattern, l_nex = read_job_lines(lines, len(l_flag))
    for pattern in l_pattern:
        pattern.title = title
        pattern.file_data = lines.next_text().strip()
    lines.next_integers(6)  # Mat Pcr NLI Rpa Sym Sho
    for pattern in l_pattern:
        pattern.ins = lines.next_integers(11)[7]
        check_data_format(pattern.ins)
    for pattern in l_pattern:
        read_wavelength_line(lines, pattern)
    lines.next_numbers(6)  # NCY Eps R_at R_an R_pr R_gl
    for pattern in l_pattern:
        pattern.ttheta_min, pattern.ttheta_step, pattern.ttheta_max = lines.next_numbers(3)[:3]
    for pattern, nex in zip(l_pattern, l_nex):
        read_background_points_and_excluded(lines, pattern, nex)
    lines.next_numbers(1)  # number of refined parameters
    for pattern in l_pattern:
        read_zero_and_background(lines, pattern)
    l_phase = []
    for _ in range(nph):
        phase = read_multi_phase(lines, l_pattern, [phase.label for phase in l_phase])
        l_phase.append(phase)
    if nre > 0:
        logger.info("The limits of %d parameters in the .pcr file are not converted.", nre)
    return title, l_pattern, take_crystal_phases(l_phase, l_pattern)


def take_crystal_phases(l_phase: List[FullProfPhase], l_pattern: List[FullProfPattern]
                        ) -> List[FullProfPhase]:
    """The phases with the cell of the first pattern they contribute to;
    other cells give a warning (CrysPy has one cell for each crystal)."""
    l_crystal = []
    for phase in l_phase:
        l_copy = [item for pattern in l_pattern for item in pattern.phases
                  if item.label == phase.label]
        if len(l_copy) == 0:
            raise FullProfFormatError(f"Phase '{phase.label:}' contributes to no pattern.")
        for item in l_copy[1:]:
            if not numpy.allclose(item.cell, l_copy[0].cell, rtol=1e-6):
                warn_conversion(
                    f"Phase '{phase.label:}' has other cell parameters in other patterns: "
                    "CrysPy keeps the cell of the first pattern.")
        l_crystal.append(l_copy[0])
    return l_crystal


def find_pattern_file(file_data: str, file_pcr: str) -> str:
    """The data file of a pattern, next to the .pcr file."""
    directory = os.path.dirname(os.path.abspath(file_pcr))
    for extension in ("", ".dat", ".xys", ".xye", ".DAT"):
        file_name = os.path.join(directory, file_data + extension)
        if os.path.isfile(file_name):
            return file_name
    raise FullProfFormatError(f"The data file {file_data:} of {file_pcr:} is not found.")


def load_fullprof_multi_pattern(file_pcr: str) -> GlobalN:
    """Read a multi-pattern .pcr file of FullProf into crystals and one
    powder experiment for each pattern."""
    _, l_pattern, l_crystal = read_fullprof_multi_pcr(file_pcr)
    ls_text = ["global_", ""]
    for phase in l_crystal:
        ls_text.extend([form_crystal_rcif(phase, calc_occupancies(phase)), ""])
    l_used = [phase.label for phase in l_crystal]
    for pattern in l_pattern:
        if pattern.job in JOBS_SIMULATION:
            ttheta = numpy.arange(pattern.ttheta_min, pattern.ttheta_max + 0.5 * pattern.ttheta_step,
                                  pattern.ttheta_step)
            intensity, sigma = numpy.zeros_like(ttheta), numpy.ones_like(ttheta)
        else:
            file_data = find_pattern_file(pattern.file_data, file_pcr)
            ttheta, intensity, sigma = read_fullprof_data(file_data, pattern.ins)
        warn_shared_codes(pattern)
        name = take_label(os.path.splitext(os.path.basename(pattern.file_data))[0], l_used)
        l_used.append(name)
        ls_text.extend([form_pd_rcif(pattern, name, ttheta, intensity, sigma,
                                     take_background(pattern, file_pcr)), ""])
    logger.info("FullProf multi-pattern file %s is converted: %d phases, %d patterns.",
                file_pcr, len(l_crystal), len(l_pattern))
    return str_to_globaln("\n".join(ls_text))


def form_multi_phase_lines(l_pattern: list, l_crystal: list, counter: CodeCounter
                           ) -> Tuple[List[str], int]:
    """The phases of a multi-pattern .pcr file (the crystals used by some
    pattern) and their number."""
    n_pattern = len(l_pattern)
    l_label = [[str(crystal.data_name).lower() for _, crystal in pattern.l_phase]
               for pattern in l_pattern]
    l_used = [crystal for crystal in l_crystal
              if any(str(crystal.data_name).lower() in labels for labels in l_label)]
    ls_phase, flag_special = [], False
    for index, crystal in enumerate(l_used):
        label = str(crystal.data_name).lower()
        l_contribution = [1 if label in labels else 0 for labels in l_label]
        ls_atom, flag = form_atom_lines(crystal, l_pattern[0].radiation, counter)
        flag_special = flag_special or flag
        ls_after_nat = [f"!Contributions (0/1) of this phase to the {n_pattern:2} patterns",
                        " " + " ".join(str(value) for value in l_contribution)]
        for i_pattern, pattern in enumerate(l_pattern):
            if l_contribution[i_pattern] == 1:
                ls_after_nat.extend([
                    f"!Irf Npr Jtyp  Nsp_Ref Ph_Shift for Pattern#  {i_pattern + 1:}",
                    f"   0   7 {pattern.job:4}      0      0",
                    f"! Pr1    Pr2    Pr3   Brind.   Rmua   Rmub   Rmuc     for Pattern#  {i_pattern + 1:}",
                    "  0.000  0.000  1.000  1.000  0.000  0.000  0.000"])
        ls_phase.extend(form_phase_head(
            crystal, index + 1, ls_atom, ls_after_nat,
            "!Nat Dis Ang Jbt Isy Str Furth        ATZ     Nvk More",
            f"{len(crystal.atom_site.items):4}   0   0   0   0   0   0          0.0000   0   0"))
        for i_pattern, pattern in enumerate(l_pattern):
            for phase, crystal_pattern in pattern.l_phase:
                if crystal_pattern is crystal:
                    ls_phase.extend(pattern.form_profile_lines(
                        phase, crystal, counter,
                        f"Pattern #  {i_pattern + 1:}  ----> Phase # {index + 1:3}"))
    warn_special_positions(flag_special)
    return ls_phase, len(l_used)


def form_fullprof_multi_pcr(l_pd: list, l_crystal: list, l_file_data_name: List[str],
                            l_ttheta: List[numpy.ndarray], title: str) -> str:
    """Text of a multi-pattern .pcr file for several powder experiments:
    the common crystals with their atoms once, the profile parameters of
    each phase for each pattern."""
    counter = CodeCounter()
    l_pattern = [PatternPcr(pd, l_crystal, ttheta, counter) for pd, ttheta in zip(l_pd, l_ttheta)]
    n_pattern = len(l_pattern)
    ls_phase, n_phase = form_multi_phase_lines(l_pattern, l_crystal, counter)
    if len({pattern.radiation for pattern in l_pattern}) > 1:
        warn_conversion(
            "The patterns have neutrons and X-rays: the scattering factors of the atoms "
            "are written for the radiation of the first pattern; check the types of the "
            "atoms in FullProf.")
    ls_out = [
        f"COMM {title:} (written by CrysPy)",
        f"NPATT {n_pattern:4}  " + " ".join(["1"] * n_pattern)
        + " <- Flags for patterns (1:refined, 0: excluded)",
        "W_PAT " + " ".join(f"{1. / n_pattern:.3f}" for _ in range(n_pattern)),
        "!Nph Dum Ias Nre Cry Opt Aut",
        f"{n_phase:4}   0   0   0   0   0   1",
        "!Job Npr Nba Nex Nsc Nor Iwg Ilo Res Ste Uni Cor Anm Int"]
    for index, pattern in enumerate(l_pattern):
        ls_out.append(
            f"{pattern.job:4}   7{len(pattern.l_bkg):4}{len(pattern.l_excluded):4}"
            f"   0   0   0   0   0   0   0   0   0   0  !-> Patt#: {index + 1:}")
    ls_out.extend(["!", "!File names of data(patterns) files"])
    ls_out.extend([f"{name:}.dat" for name in l_file_data_name])
    ls_out.extend(["!", "!Mat Pcr NLI Rpa Sym Sho", "   0   1   0   0   0   0",
                   "!Ipr Ppl Ioc Ls1 Ls2 Ls3 Prf Ins Hkl Fou Ana"])
    ls_out.extend([f"   0   0   1   0   4   0   3  10   0   0   0  !-> Patt#: {index + 1:}"
                   for index in range(n_pattern)])
    ls_out.append("!")
    for index, pattern in enumerate(l_pattern):
        ls_out.extend(pattern.form_wavelength_lines(index + 1))
    ls_out.extend(["!NCY  Eps  R_at  R_an  R_pr  R_gl", " 10  0.10  1.00  1.00  1.00  1.00"])
    for index, pattern in enumerate(l_pattern):
        ls_out.extend([f"!     Thmin       Step       Thmax    PSD    Sent0  -> Patt#: {index + 1:}",
                       f"  {pattern.form_range_numbers():}"])
    for index, pattern in enumerate(l_pattern):
        ls_out.extend(["!", *pattern.form_background_lines(index + 1)])
    ls_out.extend(["!", "NREFINED", "!"])
    for index, pattern in enumerate(l_pattern):
        ls_out.extend(pattern.form_zero_lines(index + 1))
    ls_out.extend(ls_phase)
    for index, pattern in enumerate(l_pattern):
        ls_out.extend([f"!  2Th1/TOF1    2Th2/TOF2  Pattern # {index + 1:}",
                       f"{pattern.ttheta_min:12.3f}{pattern.ttheta_max:12.3f}       1"])
    text = "\n".join(ls_out) + "\n"
    return text.replace("NREFINED", f"{counter.number:8}    !Number of refined parameters")


def save_fullprof_multi_pattern(l_pd: list, l_crystal: list, file_base: str) -> str:
    """Write several powder experiments as one multi-pattern .pcr file
    (file_base.pcr) and the data files file_base_<experiment>.dat; return
    the name of the .pcr file."""
    l_ttheta, l_name = [], []
    os.makedirs(os.path.dirname(os.path.abspath(file_base)), exist_ok=True)
    for pd in l_pd:
        ttheta, intensity, sigma = take_pd_data(pd)
        name = f"{os.path.basename(file_base):}_{pd.data_name:}"
        write_xydata(os.path.join(os.path.dirname(os.path.abspath(file_base)), name + ".dat"),
                     pd.data_name, ttheta, intensity, sigma)
        l_ttheta.append(ttheta)
        l_name.append(name)
    text = form_fullprof_multi_pcr(l_pd, l_crystal, l_name, l_ttheta,
                                   ", ".join(str(pd.data_name) for pd in l_pd))
    with open(file_base + ".pcr", "w") as fid:
        fid.write(text)
    logger.info("FullProf multi-pattern file %s.pcr is written: %d patterns.", file_base, len(l_pd))
    return file_base + ".pcr"
