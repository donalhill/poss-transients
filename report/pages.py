"""The page's text by section: each section N has its briefs in report/briefs/NN.md and its copy in report/copy/NN.md.
Joins them, in section order, into report/briefs.md and report/copy.md, the text the template, assembler, check, and
editor read.
usage: python report/pages.py"""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "report"


def build():
    for kind in ("briefs", "copy"):
        head = (REPORT / f"{kind}.head.md").read_text()
        body = "".join(f.read_text().rstrip("\n") + "\n\n" for f in sorted((REPORT / kind).glob("[0-9][0-9].md")))
        (REPORT / f"{kind}.md").write_text(head + body)


def slot_file(slot):
    """The copy file that holds a slot's text, for the editor."""
    for f in sorted((REPORT / "copy").glob("[0-9][0-9].md")):
        if re.search(rf"(?m)^### {re.escape(slot)}\s*$", f.read_text()):
            return f
    raise KeyError(slot)


if __name__ == "__main__":
    build()
    print("report/briefs.md, report/copy.md")
