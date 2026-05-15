"""Extract pages from the 4 primary papers via pypdf, UTF-8."""
import pypdf
import os

jobs = [
    ("C:/Users/Nate/OneDrive/Documents/closure hunt/AST_1995__232__243_0.pdf",
     "C:/Collatz/_voiculescu_pages"),
    ("C:/Users/Nate/OneDrive/Documents/closure hunt/memoirs.pdf",
     "C:/Collatz/_speicher_pages"),
    ("C:/Users/Nate/OneDrive/Documents/closure hunt/random-perturbations-of-matrix-cocycles.pdf",
     "C:/Collatz/_young_pages"),
    ("C:/Users/Nate/OneDrive/Documents/closure hunt/0806.0732v3.pdf",
     "C:/Collatz/_tsujii_pages"),
]

for src, dst in jobs:
    print(f"=== {src}")
    if not os.path.exists(src):
        print(f"   MISSING")
        continue
    reader = pypdf.PdfReader(src)
    n = len(reader.pages)
    print(f"   pages: {n}")
    for i in range(n):
        try:
            text = reader.pages[i].extract_text() or ""
        except Exception as e:
            text = f"[EXTRACT ERROR p{i+1}: {e}]"
        out = os.path.join(dst, f"page_{i+1:03d}.txt")
        with open(out, "w", encoding="utf-8") as f:
            f.write(text)
    print(f"   wrote {n} pages to {dst}")
