"""Extract Tao RMT PDF text page-by-page, write UTF-8 dumps for analysis."""
import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

import pypdf
from pathlib import Path

PDF_PATH = Path("C:/Users/Nate/OneDrive/Documents/profinite_transfer_operator/pdfs/Tao_Topics_Random_Matrix_Theory.pdf")
OUT_DIR = Path("C:/Collatz/_tao_rmt_pages")
OUT_DIR.mkdir(exist_ok=True, parents=True)

reader = pypdf.PdfReader(str(PDF_PATH))
n_pages = len(reader.pages)
print(f"Total pages: {n_pages}")

# Dump first 25 pages (TOC + preface) and full text in chunks
for i, page in enumerate(reader.pages):
    txt = page.extract_text() or ""
    out = OUT_DIR / f"page_{i:04d}.txt"
    out.write_text(txt, encoding="utf-8")

print(f"Wrote {n_pages} page files to {OUT_DIR}")
