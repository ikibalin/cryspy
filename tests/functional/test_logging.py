"""Refinement messages go to the "cryspy" logger, nothing to stdout."""
import logging
import os

import cryspy

F_NAME = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "examples",
    "rhochi_unpolarized_powder_1d_PbSO4", "rhochi_unpolarized_powder_1d_PbSO4.rcif")


def test_refinement_logs_and_does_not_print(capsys, caplog):
    rcif = cryspy.load_file(F_NAME)
    with caplog.at_level(logging.INFO, logger="cryspy"):
        cryspy.rhochi_rietveld_refinement(rcif)
    captured = capsys.readouterr()
    assert captured.out == ""
    messages = [record.getMessage() for record in caplog.records]
    assert all(record.name.startswith("cryspy.") for record in caplog.records)
    assert "Rietveld refinement by CrysPy (module RhoChi)" in messages
    assert "Optimal parameters:" in messages


def test_progress_lines_are_marked(caplog):
    """Iteration lines carry progress=True so that a program can replace
    the previous one instead of adding a line per iteration."""
    logger = logging.getLogger("cryspy.test")
    with caplog.at_level(logging.INFO, logger="cryspy"):
        logger.info("Iteration 1", extra=cryspy.PROGRESS)
        logger.info("Done")
    progress, done = caplog.records
    assert cryspy.is_progress_record(progress)
    assert not cryspy.is_progress_record(done)
