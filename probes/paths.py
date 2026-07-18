"""Repo-root path anchor — import this instead of hardcoding `C:/Collatz/...`.

The live L3-campaign probes already use CWD-relative paths (`outputs/...`, `logs/...`)
and run portably from the repo root: `python probes/<name>.py`. New code that needs
absolute paths should use these anchors so it runs on any machine / checkout.

    from paths import ROOT, OUTPUTS, LOGS, RESULTS, NOTES
    out = OUTPUTS / "my_dump.tsv"

NOTE: ~214 legacy probes/ scripts (the May–June R-numbered work) still contain
hardcoded `C:/Collatz` / `D:` paths and will not run off this machine without edits;
they are superseded and out of the active path (see ../RESEARCH_ARC.md).
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent   # probes/ -> repo root
OUTPUTS = ROOT / "outputs"
LOGS = ROOT / "logs"
RESULTS = ROOT / "results"
NOTES = ROOT / "notes"
PROBES = ROOT / "probes"

__all__ = ["ROOT", "OUTPUTS", "LOGS", "RESULTS", "NOTES", "PROBES"]
