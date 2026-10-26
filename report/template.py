"""The page's template, generated from the briefs (report/briefs.md, from report/pages.py) so the page's structure is the briefs' structure: one section
per outline section, each slot placed by its kind (h2: the section's question; title; text; caption: under its
figures' chart areas; verdict), in the briefs' order. The look (head, stylesheet, header, and footer) comes from
report/base_template.html and docs/style.css. Writes report/template.html and report/figure_labels.json; the charts are
drawn by docs/app.js from docs/figures.json.
usage: python report/template.py"""
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "report"))
from parse import briefs_of                                                # noqa: E402

base = (ROOT / "report/base_template.html").read_text()
head = base[: base.index("</head>")]
B = briefs_of((ROOT / "report/briefs.md").read_text())
titles = dict(re.findall(r"(?m)^## (\d+)\. (.+)$", (ROOT / "report/briefs.md").read_text()))

# the controls, above the charts that use them
STEPS = '<div class="blink-btns steps" id="{id}"><button class="back" aria-label="Previous step">&#8249; Back</button><span class="pos" aria-live="polite"></span><button class="next">Next &#8250;</button></div>'
CONTROLS = {"fig.sky_shadow": '<div class="blink-btns"><button id="shadow-play">Play the shadow plates in date order</button></div>',
            "fig.shadow_plates": '<div class="blink-btns"><button id="plates-play">Play the exposure</button></div>',
            "fig.shadow_null": STEPS.format(id="ratio-steps"), "fig.align_method": STEPS.format(id="method-steps")}
TABLES = {"fig.dates_table", "fig.vasco_repeat"}   # drawn as HTML tables, not charts
DISC, WIDE = {"fig.sky", "fig.sky_shadow"}, {"fig.shadow_plates"}              # chart classes: the whole disc on screen; the plates
# the key above the asteroid crops
AST_LEGEND = ('<div class="legend"><span><i class="sw line" style="background:repeating-linear-gradient(90deg,#c0392b 0 5px,transparent 5px 8px)"></i>'
              'Asteroid’s predicted path on the archive’s date</span><span><i class="sw line" style="background:#c0392b"></i>Its predicted path one day later</span>'
              '<span><i class="sw" style="border:1.5px solid #c0392b;background:transparent"></i>The trail, enlarged</span></div>')
figs_of = lambda b: [f.strip() for f in b.get("figure", "").split(",") if f.strip()]
is_img = lambda b: b["kind"] == "caption" and figs_of(b) and all(f.startswith("img:") for f in figs_of(b))
# figures and tables numbered in page order, as in a paper: {fig:KEY} in the copy (KEY a figure's key, e.g. fig.sky or
# img:label_xe356.png, or its caption's slot; group:SLOT for the whole group a caption belongs to) becomes a link to it (report/figrefs.py); photographs in a group are one
# figure, each labelled a, b, ...
LABELS, COUNT = {}, {"Figure": 0, "Table": 0}


def number(kind, items):
    """The next figure (or table) number: its anchor id, and the label for each (caption slot, figure keys) in items."""
    COUNT[kind] += 1; n = COUNT[kind]; anchor = f"{kind.lower()}-{n}"
    labs = [f"{kind} {n}{'abcdefgh'[i] if len(items) > 1 else ''}" for i in range(len(items))]
    for (k, keys), lab in zip(items, labs):
        for key in [k, *keys]:
            LABELS[key] = {"label": lab, "anchor": anchor}
        LABELS[f"group:{k}"] = {"label": f"{kind} {n}", "anchor": anchor}   # the whole group, from any of its captions
    return anchor, labs


def chart(k, b):
    figs = figs_of(b)
    inner = "".join(f'<div class="chart-table" data-fig="{f}"></div>' if f in TABLES      # each figure in a strip that scrolls sideways on a phone
                    else CONTROLS.get(f, "") + f'<div class="fig-scroll"><svg class="chart-svg" data-fig="{f}" viewBox="0 0 760 420" role="img" aria-label="{f}"></svg></div>' for f in figs)
    cls = " disc" if set(figs) & DISC else " plates" if set(figs) & WIDE else " tail" if set(figs) & TABLES else ""   # tables: full width, caption below
    rows = len(figs) + sum(1 for f in figs if f in CONTROLS and f not in TABLES)          # the caption spans only the rows the figures use
    anchor, (lab,) = number("Table" if set(figs) & TABLES else "Figure", [(k, figs)])
    return f'<div class="chart{cls}" id="{anchor}" style="--rows:{rows}">{inner}<div class="cap" data-label="{lab}">{{{{{k}}}}}</div></div>'


def images(group):
    """Consecutive photographs, in pairs: asteroid crops in one light box under their key, other images
    each on its own light box."""
    if len(group) == 1 and figs_of(group[0][1])[0].startswith("img:plate_"):   # a whole plate: beside its caption, as a chart
        k, b = group[0]; anchor, (lab,) = number("Figure", [(k, figs_of(b))])
        return f'<div class="chart" id="{anchor}"><img class="plate-photo" src="img/{figs_of(b)[0][4:]}" alt="{b.get("alt", "")}"><div class="cap" data-label="{lab}">{{{{{k}}}}}</div></div>'
    anchor, labs = number("Figure", [(k, figs_of(b)) for k, b in group])
    if all(figs_of(b)[0].startswith("img:asteroid_") for _, b in group):
        figs = "".join(f'<figure><img src="img/{figs_of(b)[0][4:]}" alt="{b.get("alt", "")}"><div class="lbcap" data-label="{lab}">{{{{{k}|inline}}}}</div></figure>' for (k, b), lab in zip(group, labs))
        return f'<div class="lightbox" id="{anchor}">{AST_LEGEND}<div class="plates2">{figs}</div></div>'
    boxes = "".join(f'<div class="lightbox"><img src="img/{figs_of(b)[0][4:]}" alt="{b.get("alt", "")}"><div class="lbcap" data-label="{lab}">{{{{{k}}}}}</div></div>' for (k, b), lab in zip(group, labs))
    return f'<div class="labels2" id="{anchor}">{boxes}</div>'


# the page in the briefs' order: a section for each question (h2), its first text
# beside the heading as the introduction, then steps, each a text followed by its charts; photographs grouped in pairs;
# verdicts in the verdict box
# short text in one column of reading width: two columns split it mid-sentence
COPY = dict(re.findall(r"(?ms)^### (\S+)\n(.*?)(?=^### |\Z)", (ROOT / "report/copy.md").read_text()))
cols = lambda k: "" if len(COPY.get(k, "")) > 700 else " one"
body, order = [], [(k, b) for k, b in B.items()]
hero = [(k, b) for k, b in order if b["section"] == 1]
get = lambda kind: next((k for k, b in hero if b["kind"] == kind), None)
body.append('<header class="band"><div class="page"><div class="hero">'
            f'<div><h1>{{{{{get("title")}|inline}}}}</h1></div>'
            f'<div><div class="lede">{{{{hero-lede}}}}</div>'
            '<div class="who"><img class="face" src="img/author.jpg" alt="Donal Hill" width="64" height="64">'
            '<div><div class="byline"><span>Donal Hill</span><span>October 2026</span><a href="https://github.com/donalhill/poss-transients">Code and data</a></div>'
            f'<p class="author">{{{{hero-author|inline}}}}</p></div></div></div>'
            "</div></div></header>\n<main>")
rest = [(k, b) for k, b in order if b["section"] != 1 and b["kind"] != "footer"]
foot = [k for k, b in order if b["kind"] == "footer"]
# the sources: VASCO's papers and catalogue, other work on the question, and the data the analysis uses
SOURCES = [
    "Villarroel, B. et al. (2025). Aligned, multiple-transient events in the First Palomar Sky Survey. <em>Publications of the Astronomical Society of the Pacific</em> 137, 104504.",
    "Bruehl, S. &amp; Villarroel, B. (2025). Transients in the Palomar Observatory Sky Survey (POSS-I) may be associated with nuclear testing and reports of unidentified anomalous phenomena. <em>Scientific Reports</em>.",
    "Bruehl, S., Doherty, B., Streblyanska, A. &amp; Villarroel, B. (2026). Use of machine learning to enhance detection of transient astronomical phenomena in historical observatory images. <em>Scientific Reports</em> 16, 31890.",
    "VASCO’s catalogue, used here for comparison: B. Doherty and collaborators, <a href=\"https://github.com/dca-doherty/VASCO-ML\">github.com/dca-doherty/VASCO-ML</a> (data CC BY 4.0, code MIT). Solano, E., Villarroel, B. &amp; Rodrigo, C. (2022), <em>Monthly Notices of the Royal Astronomical Society</em> 515, 1380.",
    "Watters, W. A. et al. (2026). Critical evaluation of studies alleging evidence for technosignatures in the POSS1-E photographic plates. <em>Publications of the Astronomical Society of Australia</em>.",
    "Villarroel, B., Streblyanska, A., Bruehl, S. &amp; Geier, S. (2026). A response to “Critical evaluation of studies alleging evidence for technosignatures in the POSS1-E photographic plates” by Watters et al. (2026). arXiv:2602.15171.",
    "Doherty, B. (2026), arXiv:2604.00056.",
    "Edmunds, M. G. &amp; George, G. H. (1985). On the statistics of quasar alignments. <em>Monthly Notices of the Royal Astronomical Society</em> 213, 905.",
    "Lezama, J., Randall, G., Morel, J.-M. &amp; Grompone von Gioi, R. (2015). An unsupervised point alignment detection algorithm. <em>Image Processing On Line</em> 5, 296.",
    "Test times: US Department of Energy, <em>United States Nuclear Tests, July 1945 through September 1992</em> (DOE/NV-209, rev. 16). Asteroid positions: JPL Horizons and IMCCE SkyBoT. Astrometric reference: Gaia DR3. Plate times: Palomar’s log of the survey (VizieR VI/25). Dust: Schlegel, D. J., Finkbeiner, D. P. &amp; Davis, M. (1998), <em>The Astrophysical Journal</em> 500, 525.",
    "The Digitized Sky Surveys were produced at the Space Telescope Science Institute under U.S. Government grant NAG W-2166. The images of these surveys are based on photographic data obtained using the Oschin Schmidt Telescope on Palomar Mountain and the UK Schmidt Telescope."]
open_sec, open_step, sec_now, n, titled, step_n = False, False, None, 0, False, 0
def close_step():
    global open_step
    if open_step: body.append("  </div>"); open_step = False
def close_sec():
    global open_sec
    close_step()
    if open_sec: body.append("</div></section>"); open_sec = False
def new_sec(sid, sec):
    global open_sec
    close_sec(); body.append(f'<!-- {sec}. {titles.get(str(sec), "")} -->\n<section id="{sid}" class="band"><div class="page">'); open_sec = True
while n < len(rest):
    k, b = rest[n]; kind = b["kind"]
    if kind == "h2" or b["section"] != sec_now and not open_sec or b["section"] != sec_now and kind != "h2":
        if kind == "h2":
            new_sec(k.removesuffix("-h2"), b["section"]); sec_now = b["section"]; step_n = 0
            intro = rest[n + 1] if n + 1 < len(rest) and rest[n + 1][1]["kind"] == "text" else None
            body.append(f'  <div class="sec-head stack"><div><h2>{{{{{k}|inline}}}}</h2></div>' + (f'<div class="prose{cols(intro[0])}">{{{{{intro[0]}}}}}</div>' if intro else "") + "</div>")
            n += 2 if intro else 1; continue
        new_sec(f"s{b['section']}", b["section"]); sec_now = b["section"]
    if kind == "text" and titled:                                          # the text under a step's title, in the same step
        body.append(f'    <div class="step-text{cols(k)}">{{{{{k}}}}}</div></div>'); titled = False
    elif kind == "text":
        close_step(); body.append(f'  <div class="step"><div class="step-copy"><div class="step-text{cols(k)}">{{{{{k}}}}}</div></div>'); open_step = True
    elif kind == "verdict":
        close_step()
        body.append(f'  <div class="verdict one"><div class="vlab">Result</div><div class="vtext">{{{{{k}}}}}</div></div>' if k.endswith("-verdict")
                    else f'  <div class="step"><div class="prose wide">{{{{{k}}}}}</div></div>')
    elif kind == "title":                                                   # a numbered step heading; its text follows inside
        close_step(); step_n += 1
        body.append(f'  <div class="step"><div class="step-copy"><div class="step-h"><span class="step-n">{step_n}</span><h3>{{{{{k}|inline}}}}</h3></div>'); open_step = True; titled = True
    elif is_img(b):
        group = [(k, b)]
        while n + 1 < len(rest) and is_img(rest[n + 1][1]) and figs_of(rest[n + 1][1])[0].startswith("img:asteroid_") == figs_of(b)[0].startswith("img:asteroid_"):
            n += 1; group.append(rest[n])
        if not open_step: body.append('  <div class="step">'); open_step = True
        body.append("    " + images(group))
    elif kind == "caption":
        if not open_step: body.append('  <div class="step">'); open_step = True
        body.append("    " + chart(k, b))
    elif kind == "card":                                                    # the summary: one card per test, side by side
        close_step(); cards = [k]
        while n + 1 < len(rest) and rest[n + 1][1]["kind"] == "card":
            n += 1; cards.append(rest[n][0])
        body.append('  <div class="score three" aria-label="Summary">' + "".join(f'<div class="sc"><div class="o">{{{{{c}}}}}</div></div>' for c in cards) + "</div>")
    n += 1
close_sec()
body.append("</main>")
if foot:                                                                    # the footer: about this page, and the sources
    body.append('<footer id="about" class="band"><div class="page">\n  <div class="cols">\n    <div><h2>About this page</h2>'
                + "".join(f"{{{{{k}}}}}" for k in foot) + '</div>\n    <div><h2>Sources</h2>\n      <ul>\n'
                + "".join(f"        <li>{s}</li>\n" for s in SOURCES) + '      </ul>\n    </div>\n  </div>\n</div></footer>')
head = head + "<style>.verdict.one { grid-template-columns: 1fr; } .step-text.one, .sec-head.stack .prose.one { columns: auto; max-width: 46rem; } .step > .chart + .labels2, .step > .labels2 + .lightbox, .step > .lightbox + .chart, .step > .chart + .chart, .step > .labels2 + .chart { margin-top: 1.4rem; } .plate-photo { width: 100%; max-width: 520px; display: block; border-radius: 4px; } .chart > .cap { grid-row: 1 / span var(--rows, 1); } .fig-scroll, .chart-table { overflow-x: auto; -webkit-overflow-scrolling: touch; } @media (max-width: 640px) { .fig-scroll > svg.chart-svg { min-width: 620px; } .fig-scroll { margin: 0 -0.4rem; padding-bottom: 0.3rem; } }</style>\n"
html = (head + '<link rel="stylesheet" href="style.css">\n' + "</head>\n<body>\n" + "\n".join(body) +
        '\n<div class="tip" id="tip" hidden></div>\n'
        '<script src="https://cdnjs.cloudflare.com/ajax/libs/d3/7.9.0/d3.min.js"></script>\n<script src="app.js"></script>\n</body>\n</html>\n')
out = ROOT / "report/template.html"; out.write_text(html)
(ROOT / "report/figure_labels.json").write_text(json.dumps(LABELS, indent=1))
print(f"{out.relative_to(ROOT)}: {len(B)} slots in {len(body)} blocks")
