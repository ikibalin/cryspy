"""Transformation of markdown in different formals
"""


def md_to_html(s_content: str):
    if s_content.strip() == "":
        return ""
    l_html = []
    l_content = s_content.split("\n")
    f_table_head = False
    f_ul = False
    f_ol = False

    for line in l_content:
        s_line = line.strip()
        if s_line.find("**") != -1:
            ls_line = []
            for i_val, val in enumerate(s_line.split("**")):
                if i_val%2 == 0:
                    ls_line.append(f"{val:}<b>")
                elif i_val%2 == 1:
                    ls_line.append(f"</b>{val:}")
            if i_val%2 == 0:
                ls_line.append(f"</b>")
            s_line = "".join(ls_line)

        if s_line.find("*") != -1:
            ls_line = []
            for i_val, val in enumerate(s_line.split("*")):
                if i_val%2 == 0:
                    ls_line.append(f"{val:}<i>")
                elif i_val%2 == 1:
                    ls_line.append(f"</i>{val:}")
            if i_val%2 == 0:
                ls_line.append(f"</i>")
            s_line = "".join(ls_line)

        if s_line.startswith("# "):
            l_html.append(f"<h1>{s_line.strip('# '):}</h1>\n")
        elif s_line.startswith("## "):
            l_html.append(f"<h2>{s_line.strip('# '):}</h2>\n")
        elif s_line.startswith("### "):
            l_html.append(f"<h3>{s_line.strip('# '):}</h3>\n")
        elif s_line.startswith("#### "):
            l_html.append(f"<h4>{s_line.strip('# '):}</h4>\n")
        elif s_line.startswith("|"):
            if len(set(s_line) - set("|:- ")) == 0:
                f_table_head = True
            else:
                if not(f_table_head):
                    l_html.append("<table>")
                    l_td = [f"<th>{th.strip():}</th>" for th in s_line.split("|") if th != ""]
                    f_table_head = True
                else:
                    l_td = [f"<td>{th.strip():}</td>" for th in s_line.split("|") if th != ""]
                l_html.append(f"<tr align='right'>{''.join(l_td):}</tr>")
        elif (s_line.startswith("-") | s_line.startswith("+")):
            if not(f_ul):
                f_ul = True
                l_html.append("<ul>")
            l_html.append(f"<li>{s_line[1:].strip():}</li>")
        elif (s_line.startswith("*.")):
            if not(f_ol):
                f_ol = True
                l_html.append("<ol>")
            l_html.append(f"<li>{s_line[1:].strip():}</li>")
        elif (s_line.startswith("![")):
            l_h = s_line.split("](")
            s_alt = l_h[0].lstrip("![")
            s_src = l_h[1].rstrip(")")
            l_html.append(f"<p><img alt='{s_alt:}' src='{s_src:}'></p>")
        else:
            if f_table_head:
                f_table_head = False
                l_html.append("</table>")
            if f_ul:
                f_ul = False
                l_html.append("</ul>")
            if f_ol:
                f_ol = False
                l_html.append("</ol>")
            l_html.append(f"<p>{s_line:}</p>")

    return "".join(l_html)


