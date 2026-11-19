"""Tools of the AI assistant for the documents of the user (documents.py):
search_documents, list_documents, index_documents."""
from typing import Callable, Optional

from cryspy.procedure_assistant import documents
from cryspy.procedure_assistant.tool_base import ArgumentsNeeded

MAX_FRAGMENT_CHARACTERS = 1200
MAX_TEXT_CHARACTERS = 5000  # the fragments must fit the context of a small model with the history
N_RESULTS = 4
MAX_ERRORS_SHOWN = 5
D_FOLDER_QUESTION = {
    "name": "folder", "question": "Folder with your documents (PDF, txt, md)?", "type": "string",
    "hint": "a path such as ~/papers; the subfolders are read too"}


def take_folder(folder: str = "") -> str:
    """Set the folder if it is given; ask it when none is known."""
    if folder.strip() != "":
        documents.set_document_folder(folder.strip())
    if documents.take_document_folder() is None:
        raise ArgumentsNeeded([dict(D_FOLDER_QUESTION)], "No folder of documents is given.")
    return documents.take_document_folder()


def form_errors_text(d_error: dict) -> str:
    """Files that could not be read, with the advice to install pypdf."""
    if len(d_error) == 0:
        return ""
    l_line = [f"- {name:}: {text:}" for name, text in list(d_error.items())[:MAX_ERRORS_SHOWN]]
    if len(d_error) > MAX_ERRORS_SHOWN:
        l_line.append(f"- ... and {len(d_error) - MAX_ERRORS_SHOWN:} more")
    return "Files not read:\n" + "\n".join(l_line)


def form_method_text() -> str:
    name = documents.take_document_embedding_name()
    return f"BM25 + embedding {name:}" if name else "BM25 (words only)"


def cut_text(text: str, max_length: int) -> str:
    if len(text) <= max_length:
        return text
    return text[:max_length].rsplit(" ", 1)[0] + " [...]"


def tool_search_documents(cryspy_object, query: str, n_results: int = N_RESULTS,
                          folder: str = "") -> dict:
    """Fragments of the documents with their sources ("sources" key: the
    dicts of documents.search_documents)."""
    path = take_folder(folder)
    d_summary = documents.update_document_index(flag_embed=False)
    l_found = documents.search_documents(query, n_results=max(1, min(int(n_results), 10)))
    s_errors = form_errors_text(d_summary["errors"])
    if len(l_found) == 0:
        s_text = (f"Nothing is found for '{query:}' in the {d_summary['n_files']:} documents of "
                  f"{path:}. Try other words (the words of the documents, e.g. in English).")
        return {"text": f"{s_text:}\n{s_errors:}".strip(), "sources": []}
    ls_out = [f"Fragments of the documents of the user for '{query:}', found by {form_method_text():}. "
              "Answer from them only and cite every fact as [file, page]; say so if they do "
              "not answer the question."]
    n_left = MAX_TEXT_CHARACTERS
    l_source = []
    for ind, d_found in enumerate(l_found):
        text = cut_text(d_found["text"], min(MAX_FRAGMENT_CHARACTERS, n_left))
        ls_out.append(f"[{ind + 1:}] {d_found['file']:}, {d_found['location']:}:\n{text:}")
        l_source.append(d_found)
        n_left -= len(text)
        if n_left < 200:
            break
    if s_errors:
        ls_out.append(s_errors)
    return {"text": "\n\n".join(ls_out), "sources": l_source}


def tool_list_documents(cryspy_object, folder: str = "") -> str:
    take_folder(folder)
    d_summary = documents.update_document_index(flag_embed=False)
    d_info = documents.describe_documents()
    ls_out = [f"Folder of documents: {d_info['folder']:}",
              f"{d_info['n_files']:} files, {d_info['n_chunks']:} fragments; search: {form_method_text():}."]
    if d_info["embedding_name"]:
        ls_out.append(f"Fragments with embedding vectors: {d_info['n_embedded']:} of "
                      f"{d_info['n_chunks']:} (index_documents makes the others).")
    for name, d_file in d_info["files"].items():
        s_pages = f", {d_file['n_pages']:} pages" if d_file["n_pages"] else ""
        s_error = f" [{d_file['error']:}]" if d_file["error"] else ""
        ls_out.append(f"- {name:}: {d_file['n_chunks']:} fragments{s_pages:}{s_error:}")
    s_errors = form_errors_text({name: text for name, text in d_summary["errors"].items()
                                 if name not in d_info["files"]})
    if s_errors:
        ls_out.append(s_errors)
    return "\n".join(ls_out)


def tool_index_documents(cryspy_object, folder: str = "", *,
                         progress: Optional[Callable] = None) -> str:
    path = take_folder(folder)
    d_summary = documents.update_document_index(progress=progress)
    ls_out = [f"Documents of {path:}: {d_summary['n_files']:} files, {d_summary['n_chunks']:} fragments "
              f"({d_summary['n_new']:} new, {d_summary['n_changed']:} changed, "
              f"{d_summary['n_removed']:} removed files)."]
    if documents.take_document_embedding_name():
        ls_out.append(f"Embedding vectors ({documents.take_document_embedding_name():}): "
                      f"{d_summary['n_embedded']:} of {d_summary['n_chunks']:} fragments.")
    else:
        ls_out.append("No embedding model is given: the search is by words only (BM25).")
    if d_summary["flag_stopped"]:
        ls_out.append("Stopped: call index_documents again to go on.")
    s_errors = form_errors_text(d_summary["errors"])
    if s_errors:
        ls_out.append(s_errors)
    return "\n".join(ls_out)


def tool_search_knowledge(cryspy_object, query: str, n_results: int = N_RESULTS) -> dict:
    """Fragments of the internal knowledge base of CrysPy (knowledge.py)
    with their files ("sources" key)."""
    from cryspy.procedure_assistant.knowledge import list_knowledge_files, search_knowledge
    l_found = search_knowledge(query, n_results=max(1, min(int(n_results), 10)))
    if len(l_found) == 0:
        return {"text": f"Nothing is found for '{query:}' in the knowledge base of CrysPy "
                        f"({', '.join(list_knowledge_files()):}). Try other words in English.",
                "sources": []}
    ls_out = [f"Fragments of the knowledge base of CrysPy for '{query:}', found by "
              f"{form_method_text():}. Answer from them and name the file."]
    n_left = MAX_TEXT_CHARACTERS
    for ind, d_found in enumerate(l_found):
        text = cut_text(d_found["text"], min(MAX_FRAGMENT_CHARACTERS, n_left))
        ls_out.append(f"[{ind + 1:}] knowledge/{d_found['file']:}, {d_found['location']:}:\n{text:}")
        n_left -= len(text)
        if n_left < 200:
            break
    return {"text": "\n\n".join(ls_out), "sources": l_found}
