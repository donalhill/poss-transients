"""Figure numbers on the assembled page: each caption opens with its label ("Figure 3.", "Table 1.", "Figure 2a."), and
each {fig:KEY} in the text becomes a link to that figure (labels from report/figure_labels.json, written by
template.py). An unknown key stops the build.
usage: python report/figrefs.py PAGE.html"""
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
page = Path(sys.argv[1]); html = page.read_text()
L = json.loads((ROOT / "report/figure_labels.json").read_text())


def cap(m):
    lab = f'<span class="fig-n">{m.group(2)}.</span> '
    return m.group(1) + (m.group(3) + lab if m.group(3) else lab)


html = re.sub(r'(<div class="(?:cap|lbcap)" data-label="([^"]+)">)(<p>)?', cap, html)
unknown = []


def ref(m):
    if m.group(1) not in L:
        unknown.append(m.group(1)); return m.group(0)
    x = L[m.group(1)]
    return f'<a class="fig-ref" href="#{x["anchor"]}">{x["label"]}</a>'


html = re.sub(r"\{fig:([\w.:-]+)\}", ref, html)
if unknown:
    raise SystemExit(f"figure references with no figure: {sorted(set(unknown))}")
page.write_text(html)
print(html.count('class="fig-n"'), "captions labelled,", html.count('class="fig-ref"'), "references")
