"""Transformation of markdown in different formals

The html is written for the Qt rich text subset (QLabel, QTextBrowser):
CrysPy Editor shows ``report_html()`` of the objects in its View tab
as it is. Colors are not set, so the text keeps the palette of the
application and is readable in light and dark themes.
"""
import html
import re

HTML_TABLE_OPEN = (
    "<table border='1' cellspacing='0' cellpadding='4' "
    "style='border-collapse: collapse; border-color: #808080; "
    "margin-top: 4px; margin-bottom: 8px;'>")

HTML_HEADINGS = {1: "h3", 2: "h4", 3: "h5"}

# emphasis starts after a space (or a bracket) and ends before a space
# or a punctuation mark: "N*sigma" and "(N-1)*sigma" are not emphasis
_RE_BOLD = re.compile(
    r"(?<![^\s(\[])\*\*(?=\S)([^*]+?)(?<=\S)\*\*(?![^\s.,;:!?)\]])")
_RE_ITALIC = re.compile(
    r"(?<![^\s(\[])\*(?=[^\s*])([^*]+?)(?<=[^\s*])\*(?![^\s.,;:!?)\]])")
_RE_CODE = re.compile(r"`([^`]+)`")
_RE_SEPARATOR = re.compile(r"^:?-+:?$")
_RE_NUMBER = re.compile(r"^[\s\d.,:;=/()%+\-−–…×eE]*\d[\s\d.,:;=/()%+\-−–…×eE]*$")
_RE_OL = re.compile(r"^\d+[.)]\s+")


def md_inline_to_html(s_text: str) -> str:
    """Escape the text and convert `code`, **bold** and *italic* marks.

    A single star between words or numbers (``N*sigma``) is kept as it is.
    """
    s_out = html.escape(s_text, quote=False)
    s_out = _RE_CODE.sub(r"<code>\1</code>", s_out)
    s_out = _RE_BOLD.sub(r"<b>\1</b>", s_out)
    s_out = _RE_ITALIC.sub(r"<i>\1</i>", s_out)
    return s_out


def _split_table_row(s_line: str) -> list:
    s_line = s_line.strip()
    if s_line.startswith("|"):
        s_line = s_line[1:]
    if s_line.endswith("|"):
        s_line = s_line[:-1]
    return [cell.strip() for cell in s_line.split("|")]


def _is_separator_row(l_cell: list) -> bool:
    return all([_RE_SEPARATOR.match(cell.replace(" ", "")) is not None
                for cell in l_cell if cell != ""]) and \
        any([cell != "" for cell in l_cell])


def _separator_align(cell: str):
    cell = cell.replace(" ", "")
    if cell.startswith(":") and cell.endswith(":"):
        return "center"
    if cell.endswith(":"):
        return "right"
    if cell.startswith(":"):
        return "left"
    return None


def _is_number(cell: str) -> bool:
    return _RE_NUMBER.match(cell) is not None


def md_table_to_html(l_line: list) -> str:
    """Convert the lines of a markdown table into a bordered html table.

    The first row is a header only if the second row is a separator
    (``|---|--:|``). Without the separator all rows are data rows
    (a "name | value" table). Alignment is taken from the separator,
    otherwise numerical columns are aligned right and text columns left.
    """
    ll_cell = [_split_table_row(line) for line in l_line]
    l_header, l_align = None, None
    if len(ll_cell) > 1 and _is_separator_row(ll_cell[1]):
        l_header = ll_cell[0]
        l_align = [_separator_align(cell) for cell in ll_cell[1]]
        ll_body = ll_cell[2:]
    else:
        ll_body = ll_cell
    ll_body = [l_cell for l_cell in ll_body if not _is_separator_row(l_cell)]

    n_column = max([len(l_cell) for l_cell in ll_cell])
    if l_header is not None:
        l_header = l_header + [""] * (n_column - len(l_header))
    ll_body = [l_cell + [""] * (n_column - len(l_cell)) for l_cell in ll_body]
    if l_align is None:
        l_align = []
    l_align = l_align + [None] * (n_column - len(l_align))

    for i_column in range(n_column):
        if l_align[i_column] is not None:
            continue
        l_value = [l_cell[i_column] for l_cell in ll_body
                   if l_cell[i_column] != ""]
        n_number = sum([_is_number(val) for val in l_value])
        if len(l_value) > 0 and 2 * n_number >= len(l_value):
            l_align[i_column] = "right"
        else:
            l_align[i_column] = "left"

    ls_html = [HTML_TABLE_OPEN]
    if l_header is not None:
        ls_html.append("<tr>" + "".join([
            f"<th align='{align:}'>{md_inline_to_html(cell):}</th>"
            for cell, align in zip(l_header, l_align)]) + "</tr>")
    for l_cell in ll_body:
        ls_html.append("<tr>" + "".join([
            f"<td align='{align:}'>{md_inline_to_html(cell):}</td>"
            for cell, align in zip(l_cell, l_align)]) + "</tr>")
    ls_html.append("</table>")
    return "".join(ls_html)


def md_to_html(s_content: str) -> str:
    """Convert markdown of the reports into html (Qt rich text subset).

    Supported: headings (``#``, ``##``, ``###``), tables, bullet and
    numbered lists, paragraphs (lines are kept), code blocks (```),
    images, horizontal rules and inline ``**bold**``, ``*italic*``,
    ```code```. Every opened tag is closed.
    """
    if s_content.strip() == "":
        return ""
    l_html = []
    l_table, l_paragraph = [], []
    l_code = None
    s_list = None

    def close_blocks():
        nonlocal s_list
        if len(l_table) != 0:
            l_html.append(md_table_to_html(l_table))
            l_table.clear()
        if len(l_paragraph) != 0:
            l_html.append("<p>" + "<br>".join(l_paragraph) + "</p>")
            l_paragraph.clear()
        if s_list is not None:
            l_html.append(f"</{s_list:}>")
            s_list = None

    for line in s_content.split("\n"):
        s_line = line.strip()
        if l_code is not None:
            if s_line.startswith("```"):
                l_html.append("<pre>" + html.escape("\n".join(l_code), quote=False) + "</pre>")
                l_code = None
            else:
                l_code.append(line.rstrip())
            continue
        if s_line.startswith("|"):
            if len(l_table) == 0:
                close_blocks()
            l_table.append(s_line)
            continue
        if s_line.startswith("```"):
            close_blocks()
            l_code = []
            continue
        if s_line == "":
            close_blocks()
            continue
        s_mark, _, s_title = s_line.partition(" ")
        if s_mark != "" and set(s_mark) == {"#"}:
            close_blocks()
            level = len(s_mark)
            s_title = md_inline_to_html(s_title.strip().rstrip("#").strip())
            if level in HTML_HEADINGS:
                tag = HTML_HEADINGS[level]
                l_html.append(f"<{tag:}>{s_title:}</{tag:}>")
            else:
                l_html.append(f"<p><b>{s_title:}</b></p>")
            continue
        if len(set(s_line)) == 1 and s_line[0] in "-*_" and len(s_line) >= 3:
            close_blocks()
            l_html.append("<hr>")
            continue
        flag_ul = s_line[0] in "-+*" and s_line[1:2] == " "
        match_ol = _RE_OL.match(s_line)
        if flag_ul or match_ol is not None or s_line.startswith("*."):
            s_tag = "ul" if flag_ul else "ol"
            if s_list != s_tag:
                close_blocks()
                l_html.append(f"<{s_tag:}>")
                s_list = s_tag
            if flag_ul:
                s_item = s_line[2:]
            elif match_ol is not None:
                s_item = s_line[match_ol.end():]
            else:
                s_item = s_line[2:]
            l_html.append(f"<li>{md_inline_to_html(s_item.strip()):}</li>")
            continue
        if s_line.startswith("![") and "](" in s_line:
            close_blocks()
            l_h = s_line.split("](")
            s_alt = html.escape(l_h[0][2:].rstrip("]"))
            s_src = html.escape(l_h[1].rstrip(")"))
            l_html.append(f"<p><img alt='{s_alt:}' src='{s_src:}'></p>")
            continue
        if s_list is not None or len(l_table) != 0:
            close_blocks()
        l_paragraph.append(md_inline_to_html(s_line))
    if l_code is not None:
        l_html.append("<pre>" + html.escape("\n".join(l_code), quote=False) + "</pre>")
    close_blocks()
    return "".join(l_html)
