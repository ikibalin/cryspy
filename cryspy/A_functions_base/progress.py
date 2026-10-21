"""Progress lines of the "cryspy" logger.

An iterative procedure reports each iteration with
logger.info(message, extra=PROGRESS). Such a record carries the attribute
progress=True: a program that shows the log (CrysPy Editor, a console
handler) may replace the previous progress line with the new one instead of
adding a line per iteration, as print(..., end="\\r") did before the
library moved to logging. Without such a program the record is an ordinary
INFO message.
"""
import logging

PROGRESS = {"progress": True}


def is_progress_record(record: logging.LogRecord) -> bool:
    """True if the record is a progress line that the next one replaces."""
    return bool(getattr(record, "progress", False))
