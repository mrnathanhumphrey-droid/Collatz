"""
w4_faure_extract.py — extract Faure 2009 PDF via pypdf for verbatim quote-hunt.

Targets the √3 / spectral-radius / partially expanding statement. Saves pages
as plaintext under C:/Collatz/experiments_output/faure_2009_pages/.
"""
import os
import sys
from pypdf import PdfReader

sys.stdout.reconfigure(encoding="utf-8")

PDF_PATH = r"C:/Users/Nate/OneDrive/Documents/faure_semiclassical/pdfs/Faure_2009_Semiclassical_Spectral_Gap_Partially_Expanding.pdf"
OUTDIR = r"C:/Collatz/experiments_output/faure_2009_pages"
os.makedirs(OUTDIR, exist_ok=True)


def main():
    r = PdfReader(PDF_PATH)
    n = len(r.pages)
    print(f"Faure 2009 PDF: {n} pages")
    print()
    # Save each page
    for i, pg in enumerate(r.pages):
        txt = pg.extract_text() or ""
        path = os.path.join(OUTDIR, f"page_{i+1:03d}.txt")
        with open(path, "w", encoding="utf-8") as f:
            f.write(txt)
    print(f"Saved {n} pages to {OUTDIR}")

    # Build a single concatenated text + search for √3 / sqrt(3) / 1.732 / spectral radius / partially expanding
    full = []
    for i, pg in enumerate(r.pages):
        full.append(f"=== PAGE {i+1} ===\n")
        full.append(pg.extract_text() or "")
        full.append("\n\n")
    full_str = "".join(full)
    full_path = os.path.join(OUTDIR, "_full_text.txt")
    with open(full_path, "w", encoding="utf-8") as f:
        f.write(full_str)
    print(f"Concatenated: {full_path}  ({len(full_str)} chars)")
    print()

    # Quick scan for key tokens
    tokens = ["√3", "sqrt(3)", "3^{1/2}", "3 ^ 1/2", "1.732", "1.7320", "spectral radius",
              "essential spectral radius", "partially expanding", "expanding factor",
              "anisotropic", "rho_0", "ρ_0", "ρ0", "rho 0", "Theorem", "theorem"]
    for tok in tokens:
        # case-insensitive
        idx = 0
        hits = 0
        while True:
            j = full_str.lower().find(tok.lower(), idx)
            if j < 0:
                break
            hits += 1
            if hits <= 4:
                snip = full_str[max(0, j-80):j+200].replace("\n", " ")
                print(f"  [{tok}] @ {j}: ...{snip}...")
            idx = j + 1
        if hits > 0:
            print(f"  >> {tok!r}: {hits} hits")
            print()


if __name__ == "__main__":
    main()
