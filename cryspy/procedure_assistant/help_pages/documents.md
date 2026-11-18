---
title: Search in your documents (papers, books, reports)
keywords: documents, document, papers, paper, article, literature, book, textbook, report, pdf, rag, search documents, my documents, source, reference, документы, статьи, статья, литература, книга, учебник, отчёт, отчет, источник
categories:
---
The assistant answers questions from your own documents: give a folder
with PDF, text (`.txt`) and Markdown (`.md`, `.rst`, `.tex`) files, the
subfolders are read too. The text is cut into fragments of about 180 words;
every fragment keeps its file and page (PDF) or first line (text files), so
that the answer cites the source as [file, page].

Tools:

- `search_documents` (query): the fragments that fit the question best;
- `list_documents`: the files read, their pages and fragments, files that
  could not be read (a scanned PDF without a text layer has no text);
- `index_documents`: reads the new and changed files and makes the
  embedding vectors of all fragments (slow the first time).

The search is hybrid. BM25 finds the fragments with the words of the
question (they are compared by their first 7 letters, so "magnetic" finds
"magnetization"); with an embedding model the fragments similar by meaning
are found too, also for a question in another language than the document.
Both rankings are merged (reciprocal rank fusion).

In a script:

```python
import cryspy
cryspy.set_document_folder("~/papers")
cryspy.update_document_index()
for d_found in cryspy.search_documents("extinction correction flipping ratio"):
    print(d_found["file"], d_found["location"], d_found["text"][:100])
```

CrysPy has no language model: the program gives the embedding function,
`cryspy.set_document_embedding(embed, name)` with `embed(list of texts) ->
array (n_text, n_dimension)`. CrysPy Editor: AI assistant → Documents...
(folder and embedding model: Ollama, an OpenAI-compatible server or a GGUF
file). The index is kept on disk and only changed files are read again.

PDF files need `pip install cryspy[pdf]` (pypdf).
