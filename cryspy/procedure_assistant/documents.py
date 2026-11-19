"""Search in the documents of the user (papers, books, reports) for AI assistants.

The user gives a folder with PDF, text and Markdown files
(set_document_folder). Their text is cut into fragments of about
CHUNK_WORDS words that overlap by CHUNK_OVERLAP words; each fragment keeps
its file and its page (PDF) or first line (text files), so that an answer
can cite the source. Search is hybrid (retrieval-augmented generation):

- BM25 (Robertson and Zaragoza, "The probabilistic relevance framework:
  BM25 and beyond", Found. Trends Inf. Retr. 3, 333 (2009)) on words cut to
  their first STEM_LENGTH letters (a crude stem that fits English, French,
  German and Russian endings), written with numpy only;
- optionally, the cosine similarity of embedding vectors given by the
  program (set_document_embedding): CrysPy holds no language model, the
  front end (CrysPy Editor, a script) gives a function that turns texts into
  vectors (a local GGUF model, Ollama, ...).

Both rankings are merged by reciprocal rank fusion (Cormack, Clarke and
Buettcher, SIGIR 2009: score = sum of 1 / (RRF_K + rank)).

The fragments and the embedding vectors are kept on disk (cache folder), so
only new and changed files are read again; the vectors are kept by the
hash of the text of a fragment and of the name of the embedding model.

Reading PDF needs the optional extra `pip install cryspy[pdf]` (pypdf);
text files need nothing.

    >>> import cryspy
    >>> cryspy.set_document_folder("~/papers")
    >>> for d_found in cryspy.search_documents("flipping ratio extinction"):
    ...     print(d_found["file"], d_found["location"], d_found["text"][:80])
"""
import hashlib
import json
import logging
import os
import re
from typing import Callable, Dict, List, Optional

import numpy

from cryspy.A_functions_base.optional_import import optional_module
from cryspy.procedure_assistant.help import SET_STOP_WORDS

logger = logging.getLogger(__name__)

pypdf = optional_module("pypdf", "pdf")

INDEX_VERSION = 1
INDEX_FILE = "documents_index.json"
EMBEDDING_FILE_PREFIX = "embeddings_"
SUFFIXES_TEXT = (".txt", ".md", ".markdown", ".rst", ".tex")
SUFFIX_PDF = ".pdf"
CHUNK_WORDS = 180  # about 1100 characters: a paragraph, fits a small context
CHUNK_OVERLAP = 30  # a sentence on the border is found in both fragments
STEM_LENGTH = 7  # "magnetic", "magnetization", "магнитный", "магнитного" -> same stem
BM25_K1 = 1.5
BM25_B = 0.75
RRF_K = 60
N_CANDIDATES = 50  # of each ranking, merged by reciprocal rank fusion
EMBEDDING_BATCH = 16
MAX_FILE_BYTES = 200 * 1024**2  # larger files are skipped (not documents)
RE_WORD = re.compile(r"\w+")
RE_HYPHEN_BREAK = re.compile(r"(\w)-\n(\w)")  # "magne-\ntic" in PDF text


class DocumentChunk:
    """A fragment of a document.

    file: path relative to the document folder; page: number of the page
    from 1 (PDF), 0 for text files; line: first line from 1 (text files),
    0 for PDF; text: the fragment."""

    __slots__ = ("file", "page", "line", "text", "hash")

    def __init__(self, file: str, page: int, line: int, text: str):
        self.file = file
        self.page = page
        self.line = line
        self.text = text
        # the vectors are kept by the hash of the text
        self.hash = hashlib.sha1(text.encode("utf-8")).hexdigest()

    @property
    def location(self) -> str:
        """'p. 12' or 'line 40': where the fragment is in its file."""
        if self.page > 0:
            return f"p. {self.page:}"
        return f"line {self.line:}"

    def to_list(self) -> list:
        return [self.page, self.line, self.text]


def split_index_words(text: str) -> List[str]:
    """Words of the text for BM25: lower case, stop words and single letters
    removed, cut to STEM_LENGTH letters."""
    return [word[:STEM_LENGTH] for word in RE_WORD.findall(text.lower())
            if len(word) > 1 and word not in SET_STOP_WORDS]


def split_text_chunks(text: str, chunk_words: int = CHUNK_WORDS,
                      overlap: int = CHUNK_OVERLAP) -> List[tuple]:
    """Fragments of a text of about chunk_words words, overlapping by
    overlap words.

    Returns [(first line of the fragment from 1, text of the fragment)]."""
    text = RE_HYPHEN_BREAK.sub(r"\1\2", text)
    l_word, l_line = [], []
    for ind_line, line in enumerate(text.splitlines()):
        for word in line.split():
            l_word.append(word)
            l_line.append(ind_line + 1)
    if len(l_word) == 0:
        return []
    step = max(1, chunk_words - overlap)
    l_chunk = []
    for start in range(0, len(l_word), step):
        l_chunk.append((l_line[start], " ".join(l_word[start:start + chunk_words])))
        if start + chunk_words >= len(l_word):
            break
    return l_chunk


def read_pdf_pages(path: str) -> List[str]:
    """Text of every page of a PDF file (pypdf, extra cryspy[pdf])."""
    reader = pypdf.PdfReader(path)
    l_page = []
    for page in reader.pages:
        try:
            l_page.append(page.extract_text() or "")
        except Exception as error:  # a broken page must not lose the others
            logger.debug("Page of %s is not read: %s", path, error)
            l_page.append("")
    return l_page


def read_document_chunks(path: str, file_name: str) -> List[DocumentChunk]:
    """Fragments of one document; file_name is kept in them."""
    if path.lower().endswith(SUFFIX_PDF):
        l_chunk = []
        for ind_page, text in enumerate(read_pdf_pages(path)):
            l_chunk.extend(DocumentChunk(file_name, ind_page + 1, 0, chunk)
                           for _, chunk in split_text_chunks(text))
        return l_chunk
    with open(path, "r", encoding="utf-8", errors="replace") as fid:
        text = fid.read()
    return [DocumentChunk(file_name, 0, line, chunk) for line, chunk in split_text_chunks(text)]


def is_document_file(file_name: str) -> bool:
    name = file_name.lower()
    return name.endswith(SUFFIXES_TEXT) or name.endswith(SUFFIX_PDF)


def find_document_files(folder: str) -> Dict[str, tuple]:
    """{path relative to the folder: (modification time, size)} of the
    documents in the folder and its subfolders (hidden ones skipped)."""
    d_file = {}
    for root, l_dir, l_name in os.walk(folder):
        l_dir[:] = sorted(name for name in l_dir if not name.startswith("."))
        for name in sorted(l_name):
            if name.startswith(".") or not is_document_file(name):
                continue
            path = os.path.join(root, name)
            try:
                stat = os.stat(path)
            except OSError:
                continue
            if stat.st_size > MAX_FILE_BYTES:
                continue
            d_file[os.path.relpath(path, folder).replace(os.sep, "/")] = (stat.st_mtime, stat.st_size)
    return d_file


def take_default_cache_folder(folder: str) -> str:
    """Cache folder of the index of a document folder in the cache of the
    user (~/.cache/cryspy, %LOCALAPPDATA%\\cryspy on Windows)."""
    base = os.environ.get("LOCALAPPDATA") if os.name == "nt" else None
    if not base:
        base = os.environ.get("XDG_CACHE_HOME") or os.path.join(os.path.expanduser("~"), ".cache")
    s_hash = hashlib.sha1(os.path.abspath(folder).encode("utf-8")).hexdigest()[:12]
    return os.path.join(base, "cryspy", "documents", s_hash)


class Bm25Index:
    """BM25 of the fragments with numpy: for every word the fragments it
    is in and how many times (sorted by word, as a sparse matrix)."""

    def __init__(self, l_text: List[str]):
        d_word = {}
        l_word_id, l_chunk_id, l_count = [], [], []
        length = numpy.zeros(len(l_text), dtype=float)
        for ind_chunk, text in enumerate(l_text):
            d_count = {}
            l_word = split_index_words(text)
            length[ind_chunk] = len(l_word)
            for word in l_word:
                d_count[word] = d_count.get(word, 0) + 1
            for word, count in d_count.items():
                l_word_id.append(d_word.setdefault(word, len(d_word)))
                l_chunk_id.append(ind_chunk)
                l_count.append(count)
        word_id = numpy.array(l_word_id, dtype=numpy.int64)
        order = numpy.argsort(word_id, kind="stable")
        self.d_word = d_word
        self.chunk_id = numpy.array(l_chunk_id, dtype=numpy.int64)[order]
        self.count = numpy.array(l_count, dtype=float)[order]
        # fragments of word i: self.chunk_id[start[i]:start[i + 1]]
        self.start = numpy.searchsorted(word_id[order], numpy.arange(len(d_word) + 1))
        self.length = length
        self.n_chunk = len(l_text)

    def score(self, query: str) -> numpy.ndarray:
        """BM25 score of every fragment for the query."""
        score = numpy.zeros(self.n_chunk, dtype=float)
        if self.n_chunk == 0:
            return score
        mean_length = max(float(self.length.mean()), 1.)
        for word in set(split_index_words(query)):
            ind_word = self.d_word.get(word)
            if ind_word is None:
                continue
            sl = slice(self.start[ind_word], self.start[ind_word + 1])
            chunk_id, count = self.chunk_id[sl], self.count[sl]
            n_with = len(chunk_id)
            idf = numpy.log(1. + (self.n_chunk - n_with + 0.5) / (n_with + 0.5))
            norm = BM25_K1 * (1. - BM25_B + BM25_B * self.length[chunk_id] / mean_length)
            numpy.add.at(score, chunk_id, idf * count * (BM25_K1 + 1.) / (count + norm))
        return score


class DocumentIndex:
    """Fragments of the documents of a folder, their BM25 index and the
    embedding vectors (when an embedding function is given)."""

    def __init__(self, folder: str, cache_folder: Optional[str] = None):
        self.folder = os.path.abspath(os.path.expanduser(folder))
        self.cache_folder = cache_folder or take_default_cache_folder(self.folder)
        self.d_file: Dict[str, dict] = {}  # {file: {"mtime", "size", "chunks", "error"}}
        self.l_chunk: List[DocumentChunk] = []
        self.bm25: Optional[Bm25Index] = None
        self.d_vector: Dict[str, numpy.ndarray] = {}  # {hash of the text: vector}
        self.embedding_name = ""
        self.flag_vectors_changed = False
        self.vector_matrix = (None, [], None)
        self.load()

    # files on disk
    def take_index_path(self) -> str:
        return os.path.join(self.cache_folder, INDEX_FILE)

    def take_embedding_path(self, name: str) -> str:
        s_hash = hashlib.sha1(name.encode("utf-8")).hexdigest()[:12]
        return os.path.join(self.cache_folder, f"{EMBEDDING_FILE_PREFIX:}{s_hash:}.npz")

    def load(self):
        try:
            with open(self.take_index_path(), "r", encoding="utf-8") as fid:
                d_index = json.load(fid)
        except (OSError, ValueError):
            return
        if d_index.get("version") != INDEX_VERSION or d_index.get("folder") != self.folder:
            return
        self.d_file = d_index.get("files", {})
        self.form_chunks()

    def save(self):
        os.makedirs(self.cache_folder, exist_ok=True)
        path = self.take_index_path()
        with open(path + ".tmp", "w", encoding="utf-8") as fid:
            json.dump({"version": INDEX_VERSION, "folder": self.folder, "files": self.d_file}, fid)
        os.replace(path + ".tmp", path)

    def form_chunks(self):
        self.l_chunk = [DocumentChunk(file_name, page, line, text)
                        for file_name in sorted(self.d_file)
                        for page, line, text in self.d_file[file_name].get("chunks", [])]
        self.bm25 = None

    def take_bm25(self) -> Bm25Index:
        if self.bm25 is None:
            self.bm25 = Bm25Index([chunk.text for chunk in self.l_chunk])
        return self.bm25

    def update(self, progress: Optional[Callable] = None) -> dict:
        """Read new and changed documents, forget removed ones.

        progress(info): optional, called before every file read, with info
        {"label": "Reading documents", "step", "n_step", "file"}; returning
        False stops the update (the files read are kept).

        Returns {"n_new", "n_changed", "n_removed", "errors": {file: text}}."""
        if not os.path.isdir(self.folder):
            raise FileNotFoundError(f"The document folder '{self.folder:}' is not found.")
        d_found = find_document_files(self.folder)
        l_removed = [name for name in self.d_file if name not in d_found]
        for name in l_removed:
            del self.d_file[name]
        l_read = [name for name, (mtime, size) in d_found.items()
                  if name not in self.d_file or self.d_file[name].get("mtime") != mtime
                  or self.d_file[name].get("size") != size]
        n_new = sum(name not in self.d_file for name in l_read)
        d_error = {}
        for ind, name in enumerate(l_read):
            if progress is not None and progress(
                    {"label": "Reading documents", "step": ind + 1, "n_step": len(l_read),
                     "file": name}) is False:
                break
            mtime, size = d_found[name]
            d_record = {"mtime": mtime, "size": size, "chunks": [], "error": ""}
            try:
                l_chunk = read_document_chunks(os.path.join(self.folder, name), name)
                d_record["chunks"] = [chunk.to_list() for chunk in l_chunk]
                if len(l_chunk) == 0:
                    d_record["error"] = "no text (a scanned PDF without a text layer?)"
            except ImportError as error:  # pypdf is not installed: read it later
                d_error[name] = str(error)
                continue
            except Exception as error:
                d_record["error"] = f"{type(error).__name__:}: {error:}"
            if d_record["error"]:
                d_error[name] = d_record["error"]
            self.d_file[name] = d_record
        if l_read or l_removed:
            self.form_chunks()
            self.save()
        return {"n_new": n_new, "n_changed": len(l_read) - n_new, "n_removed": len(l_removed),
                "errors": d_error}

    # embedding vectors
    def load_vectors(self, name: str):
        """Vectors of the embedding model 'name' from the cache."""
        if name == self.embedding_name:
            return
        self.save_vectors()
        self.embedding_name, self.d_vector, self.flag_vectors_changed = name, {}, False
        if name == "":
            return
        try:
            with numpy.load(self.take_embedding_path(name), allow_pickle=False) as data:
                self.d_vector = dict(zip((str(key) for key in data["hashes"]), data["vectors"]))
        except (OSError, ValueError, KeyError):
            pass

    def save_vectors(self):
        if not self.flag_vectors_changed or self.embedding_name == "":
            return
        set_hash = {chunk.hash for chunk in self.l_chunk}
        l_hash = [key for key in self.d_vector if key in set_hash]  # removed files are forgotten
        os.makedirs(self.cache_folder, exist_ok=True)
        path = self.take_embedding_path(self.embedding_name)
        vectors = (numpy.array([self.d_vector[key] for key in l_hash], dtype=numpy.float32)
                   if l_hash else numpy.zeros((0, 0), dtype=numpy.float32))
        with open(path + ".tmp", "wb") as fid:
            numpy.savez(fid, hashes=numpy.array(l_hash, dtype=str), vectors=vectors)
        os.replace(path + ".tmp", path)
        self.flag_vectors_changed = False

    def embed_chunks(self, l_ind: List[int], embed: Callable,
                     progress: Optional[Callable] = None) -> bool:
        """Vectors of the fragments l_ind that have none; False if stopped."""
        l_todo = [ind for ind in l_ind if self.l_chunk[ind].hash not in self.d_vector]
        for start in range(0, len(l_todo), EMBEDDING_BATCH):
            if progress is not None and progress(
                    {"label": "Embedding fragments", "step": start, "n_step": len(l_todo)}) is False:
                self.save_vectors()
                return False
            l_batch = [self.l_chunk[ind] for ind in l_todo[start:start + EMBEDDING_BATCH]]
            vectors = take_unit_vectors(embed([chunk.text for chunk in l_batch]))
            if vectors.shape[0] != len(l_batch):
                raise ValueError("The embedding function gave a wrong number of vectors.")
            for chunk, vector in zip(l_batch, vectors):
                self.d_vector[chunk.hash] = vector
            self.flag_vectors_changed = True
        self.save_vectors()
        return True

    def take_vector_matrix(self) -> tuple:
        """(indexes of the fragments with vectors, their vectors as rows),
        made again only when the fragments or the vectors change."""
        key = (id(self.l_chunk), len(self.d_vector), self.embedding_name)
        if self.vector_matrix[0] != key:
            l_ind = [ind for ind, chunk in enumerate(self.l_chunk) if chunk.hash in self.d_vector]
            matrix = (numpy.array([self.d_vector[self.l_chunk[ind].hash] for ind in l_ind])
                      if l_ind else numpy.zeros((0, 0), dtype=numpy.float32))
            self.vector_matrix = (key, l_ind, matrix)
        return self.vector_matrix[1], self.vector_matrix[2]

    def count_embedded(self) -> int:
        return sum(chunk.hash in self.d_vector for chunk in self.l_chunk)

    # search
    def search(self, query: str, n_results: int = 5, embed: Optional[Callable] = None,
               embed_query: Optional[Callable] = None) -> List[dict]:
        """Best fragments for the query, the best first (see search_documents)."""
        if len(self.l_chunk) == 0:
            return []
        score_bm25 = self.take_bm25().score(query)
        l_bm25 = [int(ind) for ind in numpy.argsort(-score_bm25, kind="stable")[:N_CANDIDATES]
                  if score_bm25[ind] > 0.]
        d_rrf = {ind: 1. / (RRF_K + rank + 1) for rank, ind in enumerate(l_bm25)}
        d_rank_embedding = {}
        if embed is not None:
            # the BM25 candidates get vectors now, other fragments after index_documents
            self.embed_chunks(l_bm25, embed)
            l_ind, matrix = self.take_vector_matrix()
            if l_ind:
                query_vector = take_unit_vectors((embed_query or embed)([query]))[0]
                if matrix.shape[1] == query_vector.shape[0]:
                    similarity = matrix @ query_vector
                    for rank, ind_sorted in enumerate(numpy.argsort(-similarity)[:N_CANDIDATES]):
                        ind = l_ind[int(ind_sorted)]
                        d_rank_embedding[ind] = rank + 1
                        d_rrf[ind] = d_rrf.get(ind, 0.) + 1. / (RRF_K + rank + 1)
        d_rank_bm25 = {ind: rank + 1 for rank, ind in enumerate(l_bm25)}
        l_best = sorted(d_rrf, key=lambda ind: -d_rrf[ind])[:n_results]
        return [{"file": self.l_chunk[ind].file, "page": self.l_chunk[ind].page,
                 "line": self.l_chunk[ind].line, "location": self.l_chunk[ind].location,
                 "text": self.l_chunk[ind].text, "score": float(d_rrf[ind]),
                 "rank_bm25": d_rank_bm25.get(ind), "rank_embedding": d_rank_embedding.get(ind)}
                for ind in l_best]


def take_unit_vectors(vectors) -> numpy.ndarray:
    """Vectors as rows of unit length (cosine similarity = dot product)."""
    matrix = numpy.atleast_2d(numpy.asarray(vectors, dtype=numpy.float32))
    norm = numpy.linalg.norm(matrix, axis=1, keepdims=True)
    return matrix / numpy.where(norm > 0, norm, 1.)


D_DOCUMENTS = {"folder": None, "cache_folder": None, "index": None,
               "embed": None, "embed_query": None, "embedding_name": ""}


def set_document_folder(folder: Optional[str], cache_folder: Optional[str] = None):
    """Folder with the documents of the user (PDF, txt, md); None: no
    documents. cache_folder keeps the index (default: the cache folder of
    the user, see take_default_cache_folder)."""
    if folder is not None:
        folder = os.path.abspath(os.path.expanduser(str(folder)))
    if folder != D_DOCUMENTS["folder"] or cache_folder != D_DOCUMENTS["cache_folder"]:
        close_document_index()
    D_DOCUMENTS["folder"] = folder
    D_DOCUMENTS["cache_folder"] = cache_folder


def take_document_folder() -> Optional[str]:
    return D_DOCUMENTS["folder"]


def set_document_embedding(embed: Optional[Callable], name: str = "",
                           embed_query: Optional[Callable] = None):
    """Function of the program that turns texts into vectors for the
    search by meaning; None: BM25 only.

    embed(list of texts) -> array (n_text, n_dimension) for the fragments,
    embed_query the same for the question (default: embed; models such as
    nomic-embed-text or e5 want different prefixes for them). name names
    the model: vectors of different models are kept apart."""
    index = D_DOCUMENTS["index"]
    if index is not None:
        index.save_vectors()
    D_DOCUMENTS["embed"] = embed
    D_DOCUMENTS["embed_query"] = embed_query
    D_DOCUMENTS["embedding_name"] = name if embed is not None else ""


def take_document_embedding_name() -> str:
    return D_DOCUMENTS["embedding_name"]


def close_document_index():
    """Save the vectors and forget the index (it is read again when needed)."""
    index = D_DOCUMENTS["index"]
    if index is not None:
        index.save_vectors()
    D_DOCUMENTS["index"] = None


def take_document_index(flag_update: bool = True, progress: Optional[Callable] = None) -> DocumentIndex:
    """The index of the document folder, updated with the changed files."""
    folder = D_DOCUMENTS["folder"]
    if folder is None:
        raise ValueError("No document folder is given (cryspy.set_document_folder).")
    index = D_DOCUMENTS["index"]
    if index is None:
        index = DocumentIndex(folder, D_DOCUMENTS["cache_folder"])
        D_DOCUMENTS["index"] = index
    if flag_update:
        index.update(progress=progress)
    index.load_vectors(D_DOCUMENTS["embedding_name"])
    return index


def update_document_index(progress: Optional[Callable] = None, flag_embed: bool = True) -> dict:
    """Read new and changed documents and, with an embedding function, make
    the vectors of all fragments (slow the first time: once per fragment).

    progress(info): see DocumentIndex.update; returning False stops.

    Returns the summary of DocumentIndex.update with "n_files", "n_chunks",
    "n_embedded" and "flag_stopped"."""
    flag_stopped = {"value": False}

    def take_progress(info):
        if progress is None:
            return True
        flag_go = progress(info) is not False
        flag_stopped["value"] = flag_stopped["value"] or not flag_go
        return flag_go
    index = take_document_index(flag_update=False)
    d_summary = index.update(progress=take_progress)
    index.load_vectors(D_DOCUMENTS["embedding_name"])
    embed = D_DOCUMENTS["embed"]
    if flag_embed and embed is not None and not flag_stopped["value"]:
        if not index.embed_chunks(list(range(len(index.l_chunk))), embed, progress=take_progress):
            flag_stopped["value"] = True
    d_summary.update(n_files=len(index.d_file), n_chunks=len(index.l_chunk),
                     n_embedded=index.count_embedded(), flag_stopped=flag_stopped["value"])
    return d_summary


def search_documents(query: str, n_results: int = 5) -> List[dict]:
    """Fragments of the documents of the user that fit the query best.

    The index is updated with the changed files first. The ranking is BM25
    merged with the similarity of embedding vectors when an embedding
    function is given (set_document_embedding).

    Returns [{"file", "page" (0 for text files), "line" (0 for PDF),
    "location" ("p. 3" or "line 40"), "text", "score", "rank_bm25",
    "rank_embedding" (None when the fragment is not in that ranking)}]."""
    index = take_document_index()
    return index.search(query, n_results=n_results, embed=D_DOCUMENTS["embed"],
                        embed_query=D_DOCUMENTS["embed_query"])


def describe_documents() -> dict:
    """{"folder", "n_files", "n_chunks", "n_embedded", "embedding_name",
    "files": {file: {"n_chunks", "n_pages", "error"}}} of the index."""
    index = take_document_index()
    d_out = {}
    for name, d_record in index.d_file.items():
        l_chunk = d_record.get("chunks", [])
        d_out[name] = {"n_chunks": len(l_chunk),
                       "n_pages": max((page for page, _, _ in l_chunk), default=0),
                       "error": d_record.get("error", "")}
    return {"folder": index.folder, "n_files": len(index.d_file), "n_chunks": len(index.l_chunk),
            "n_embedded": index.count_embedded(), "embedding_name": D_DOCUMENTS["embedding_name"],
            "files": d_out}
