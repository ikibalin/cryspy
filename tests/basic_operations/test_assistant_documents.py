"""Search in the documents of the user (procedure_assistant/documents.py)."""
import math
import os

import numpy
import pytest

import cryspy
from cryspy.procedure_assistant import documents
from cryspy.procedure_assistant.documents import (
    Bm25Index, DocumentIndex, split_index_words, split_text_chunks)

TEXT_EXTINCTION = """Extinction correction

The extinction of strong reflections reduces their intensity. The model of
Becker and Coppens describes the mosaic spread and the size of domains.
"""
TEXT_MAGNETIC = """Magnetization density

The magnetization density of the crystal is reconstructed by the maximum
entropy method from the flipping ratios of polarized neutrons.
"""
TEXT_RUSSIAN = """Магнитная структура

Магнитные моменты атомов тербия упорядочены вдоль оси c.
"""


@pytest.fixture(autouse=True)
def no_documents():
    yield
    cryspy.set_document_folder(None)
    cryspy.set_document_embedding(None)


@pytest.fixture
def folder(tmp_path):
    path = tmp_path / "papers"
    (path / "sub").mkdir(parents=True)
    (path / "extinction.txt").write_text(TEXT_EXTINCTION, encoding="utf-8")
    (path / "sub" / "magnetic.md").write_text(TEXT_MAGNETIC, encoding="utf-8")
    (path / "russian.txt").write_text(TEXT_RUSSIAN, encoding="utf-8")
    (path / "data.rcif").write_text("data_x\n_cell_length_a 5\n", encoding="utf-8")
    (path / ".hidden.txt").write_text("extinction", encoding="utf-8")
    return path


def test_split_text_chunks_overlap_and_lines():
    text = "\n".join(f"w{ind:}" for ind in range(25))  # a word per line
    l_chunk = split_text_chunks(text, chunk_words=10, overlap=3)
    assert [line for line, _ in l_chunk] == [1, 8, 15, 22]
    assert l_chunk[0][1].split()[-3:] == l_chunk[1][1].split()[:3]
    assert l_chunk[-1][1].split()[-1] == "w24"
    assert split_text_chunks("  \n ") == []


def test_split_text_chunks_joins_hyphen_break():
    assert split_text_chunks("magne-\ntization density")[0][1] == "magnetization density"


def test_split_index_words_stems_and_stop_words():
    assert split_index_words("The magnetic magnetization") == ["magneti", "magneti"]
    assert split_index_words("магнитного магнитный") == ["магнитн", "магнитн"]


def calc_bm25_reference(l_text, query):
    """BM25 written term by term, as in Robertson and Zaragoza (2009)."""
    l_doc = [split_index_words(text) for text in l_text]
    mean_length = sum(len(doc) for doc in l_doc) / len(l_doc)
    l_score = []
    for doc in l_doc:
        score = 0.
        for word in set(split_index_words(query)):
            n_with = sum(word in doc_2 for doc_2 in l_doc)
            if n_with == 0:
                continue
            idf = math.log(1. + (len(l_doc) - n_with + 0.5) / (n_with + 0.5))
            count = doc.count(word)
            score += idf * count * 2.5 / (count + 1.5 * (0.25 + 0.75 * len(doc) / mean_length))
        l_score.append(score)
    return l_score


def test_bm25_matches_reference():
    l_text = [TEXT_EXTINCTION, TEXT_MAGNETIC, TEXT_RUSSIAN, "extinction extinction model"]
    query = "extinction model of magnetization"
    score = Bm25Index(l_text).score(query)
    assert numpy.allclose(score, calc_bm25_reference(l_text, query), rtol=1e-12)
    assert numpy.allclose(Bm25Index([]).score(query), [])


def test_search_documents_by_words(folder):
    cryspy.set_document_folder(str(folder), cache_folder=str(folder.parent / "cache"))
    l_found = cryspy.search_documents("mosaic spread of domains")
    assert l_found[0]["file"] == "extinction.txt"
    assert l_found[0]["location"] == "line 1"
    assert l_found[0]["rank_bm25"] == 1 and l_found[0]["rank_embedding"] is None
    # stems: another ending of the Russian word
    assert cryspy.search_documents("магнитного момента")[0]["file"] == "russian.txt"
    assert cryspy.search_documents("density")[0]["file"] == "sub/magnetic.md"
    d_info = documents.describe_documents()
    assert sorted(d_info["files"]) == ["extinction.txt", "russian.txt", "sub/magnetic.md"]
    assert cryspy.search_documents("zzzz") == []


def test_index_is_cached_and_updated(folder, monkeypatch):
    cache = str(folder.parent / "cache")
    index = DocumentIndex(str(folder), cache)
    assert index.update()["n_new"] == 3

    def fail(*args):
        raise AssertionError("an unchanged file is read again")
    monkeypatch.setattr(documents, "read_document_chunks", fail)
    index_2 = DocumentIndex(str(folder), cache)
    assert len(index_2.l_chunk) == len(index.l_chunk)
    assert index_2.update() == {"n_new": 0, "n_changed": 0, "n_removed": 0, "errors": {}}
    monkeypatch.undo()

    (folder / "russian.txt").unlink()
    path = folder / "extinction.txt"
    path.write_text(TEXT_EXTINCTION + "\nPrimary and secondary extinction.\n", encoding="utf-8")
    os.utime(path, (1, 1))  # a new modification time on any file system
    d_summary = DocumentIndex(str(folder), cache).update()
    assert (d_summary["n_changed"], d_summary["n_removed"]) == (1, 1)
    l_found = DocumentIndex(str(folder), cache).search("secondary")
    assert l_found[0]["file"] == "extinction.txt"


L_CONCEPTS = (("extinction", "mosaic", "погасан"), ("magnet", "магнит"), ("entropy", "энтроп"))


def embed_concepts(l_text):
    """A fake embedding model: a vector of the concepts in the text, the
    same for English and Russian words."""
    l_vector = []
    for text in l_text:
        text = text.lower()
        l_vector.append([sum(text.count(word) for word in concept) + 0.01 for concept in L_CONCEPTS])
    return numpy.array(l_vector)


def test_embedding_finds_text_in_another_language(folder):
    cache = folder.parent / "cache"
    cryspy.set_document_folder(str(folder), cache_folder=str(cache))
    query = "погасание"  # 'extinction' in Russian: no common word with the documents
    assert cryspy.search_documents(query) == []
    n_call = []

    def embed(l_text):
        n_call.append(len(l_text))
        return embed_concepts(l_text)
    cryspy.set_document_embedding(embed, "concepts")
    d_summary = cryspy.update_document_index()
    assert d_summary["n_embedded"] == d_summary["n_chunks"] == 3
    l_found = cryspy.search_documents(query)
    assert l_found[0]["file"] == "extinction.txt"
    assert l_found[0]["rank_embedding"] == 1 and l_found[0]["rank_bm25"] is None
    # the vectors are read from the cache: only the query is embedded
    cryspy.set_document_folder(None)
    cryspy.set_document_folder(str(folder), cache_folder=str(cache))
    n_call.clear()
    assert cryspy.search_documents(query)[0]["file"] == "extinction.txt"
    assert n_call == [1]
    # another model does not use these vectors
    cryspy.set_document_embedding(embed, "other model")
    assert documents.describe_documents()["n_embedded"] == 0


def test_hybrid_ranking_merges_both(folder):
    cryspy.set_document_folder(str(folder), cache_folder=str(folder.parent / "cache"))
    cryspy.set_document_embedding(embed_concepts, "concepts")
    # BM25 finds the candidates' vectors itself, without index_documents
    l_found = cryspy.search_documents("maximum entropy magnetization")
    assert l_found[0]["file"] == "sub/magnetic.md"
    assert l_found[0]["rank_bm25"] == 1 and l_found[0]["rank_embedding"] == 1


def test_tool_search_documents_asks_folder_and_cites(folder):
    result = cryspy.call_assistant_tool(None, "search_documents", {"query": "extinction"})
    assert result["missing"][0]["name"] == "folder"
    result = cryspy.call_assistant_tool(
        None, "search_documents", {"query": "extinction of strong reflections", "folder": str(folder)})
    assert not result["flag_error"]
    assert "[1] extinction.txt, line 1:" in result["text"]
    assert "cite" in result["text"]
    assert result["sources"][0]["file"] == "extinction.txt"
    result = cryspy.call_assistant_tool(None, "list_documents", {})
    assert "3 files" in result["text"] and "sub/magnetic.md" in result["text"]
    result = cryspy.call_assistant_tool(None, "index_documents", {})
    assert "3 files" in result["text"] and "BM25" in result["text"]


def test_tool_index_documents_reports_progress_and_stops(folder):
    cryspy.set_document_folder(str(folder), cache_folder=str(folder.parent / "cache"))
    l_info = []

    def progress(info):
        l_info.append(info)
        return False
    result = cryspy.call_assistant_tool(None, "index_documents", {}, progress=progress)
    assert "Stopped" in result["text"]
    assert l_info[0]["label"] == "Reading documents" and l_info[0]["n_step"] == 3


def test_missing_folder_is_an_error(tmp_path):
    result = cryspy.call_assistant_tool(None, "search_documents",
                                        {"query": "x", "folder": str(tmp_path / "none")})
    assert result["flag_error"] and "not found" in result["text"]


def make_pdf(l_page_text) -> bytes:
    """A minimal PDF with one line of text on every page."""
    n_page = len(l_page_text)
    l_object = ["<< /Type /Catalog /Pages 2 0 R >>",
                "<< /Type /Pages /Kids [" + " ".join(f"{3 + 2 * ind:} 0 R" for ind in range(n_page))
                + f"] /Count {n_page:} >>"]
    for ind, text in enumerate(l_page_text):
        stream = f"BT /F1 12 Tf 72 720 Td ({text:}) Tj ET"
        l_object.append(f"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 612 792] "
                        f"/Contents {4 + 2 * ind:} 0 R /Resources << /Font << /F1 {3 + 2 * n_page:} 0 R >> >> >>")
        l_object.append(f"<< /Length {len(stream):} >>\nstream\n{stream:}\nendstream")
    l_object.append("<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>")
    out = b"%PDF-1.4\n"
    l_offset = []
    for ind, s_object in enumerate(l_object):
        l_offset.append(len(out))
        out += f"{ind + 1:} 0 obj\n{s_object:}\nendobj\n".encode("latin-1")
    xref = len(out)
    out += f"xref\n0 {len(l_object) + 1:}\n0000000000 65535 f \n".encode("latin-1")
    out += "".join(f"{offset:010d} 00000 n \n" for offset in l_offset).encode("latin-1")
    out += (f"trailer\n<< /Size {len(l_object) + 1:} /Root 1 0 R >>\n"
            f"startxref\n{xref:}\n%%EOF\n").encode("latin-1")
    return out


def test_pdf_pages_are_cited(tmp_path):
    pytest.importorskip("pypdf")
    (tmp_path / "paper.pdf").write_bytes(make_pdf(
        ["Introduction to neutron scattering", "The flipping ratio of TbTe3 was measured"]))
    cryspy.set_document_folder(str(tmp_path), cache_folder=str(tmp_path / ".cache"))
    l_found = cryspy.search_documents("flipping ratio")
    assert (l_found[0]["file"], l_found[0]["location"]) == ("paper.pdf", "p. 2")


def test_pdf_without_pypdf_is_reported(tmp_path, monkeypatch):
    (tmp_path / "paper.pdf").write_bytes(make_pdf(["text"]))

    def no_pypdf(path):
        raise ImportError("'pypdf' is required ... 'pip install cryspy[pdf]'.")
    monkeypatch.setattr(documents, "read_pdf_pages", no_pypdf)
    cryspy.set_document_folder(str(tmp_path), cache_folder=str(tmp_path / ".cache"))
    result = cryspy.call_assistant_tool(None, "search_documents", {"query": "text"})
    assert "Nothing is found" in result["text"] and "cryspy[pdf]" in result["text"]


def test_help_page_of_documents():
    assert "search_documents" in cryspy.get_help("search in my papers pdf")
