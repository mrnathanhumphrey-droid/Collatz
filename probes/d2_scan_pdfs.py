"""D2 helper: read selected PDFs from closure-hunt corpus and dump text by page range."""
import sys, pypdf, os
sys.stdout.reconfigure(encoding='utf-8')

FOLDER = 'C:/Collatz/_d2_pdfs'

def dump(name, start, end, char_limit=3500):
    p = os.path.join(FOLDER, name)
    try:
        r = pypdf.PdfReader(p)
    except Exception as e:
        print(f'ERR {name}: {e}')
        return
    n = len(r.pages)
    print(f'\n###### FILE: {name} ({n} pages) ######\n')
    end = min(end, n)
    for i in range(start, end):
        print(f'=== page {i+1} of {name} ===')
        try:
            txt = r.pages[i].extract_text() or ''
            print(txt[:char_limit])
        except Exception as e:
            print(f'(extract error: {e})')
        print()

if __name__ == '__main__':
    out_path = sys.argv[1]
    sys.stdout = open(out_path, 'w', encoding='utf-8')
    fname = sys.argv[2]
    start = int(sys.argv[3])
    end = int(sys.argv[4])
    char_limit = int(sys.argv[5]) if len(sys.argv) >= 6 else 5000
    dump(fname, start, end, char_limit)
    sys.stdout.close()
