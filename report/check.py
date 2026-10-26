"""Holds the page's copy (report/copy.md, from report/pages.py) to the analysis, both ways, and to the house style.
Errors stop the build; notes are listed for the next drafting pass.
- Structure: every slot of the copy has a brief.
- Text to measurement: every {n:key} in the copy is a registry key (report/registry.py) with a value, and every key a
  brief names is quoted in that slot.
- Measurement to text: every registry key is named by a brief.
- Figures: every caption's figure is built, and every built figure has a caption.
- Style: no em dashes, no American spellings, no literal numbers outside {n:} (years, catalogue names, and the whitelist
  below excepted), and a note on lists of three without the Oxford comma.
usage: python report/check.py"""
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "report"))
from parse import briefs_of, slots_of                                         # noqa: E402
from registry import REG                                                     # noqa: E402

AMERICAN = r"\b(\w+iz(e|es|ed|ing|ation)|color\w*|center\w*|catalog|catalogs|analyz\w*|modeling|labeled|traveled|gray|meter|meters|defense)\b"
NOT_AMERICAN = {"size", "sizes", "sized", "prize", "seize", "citizen", "citizens", "capsize"}
STEERING = r"\b(surprisingly|interestingly|notably|importantly|crucially|clearly|striking|strikingly|remarkabl[ey]|dramatic\w*|robust|delve|leverage|landscape|tapestry|underscore\w*|showcase\w*|seamless\w*|holistic|paradigm|of course)\b"
ALLOWED_DIGITS = r"\b(68|90|95)%|\b(1[89]\d\d|20\d\d)s?\b|\b48-inch\b|\bDR[0-9]\b|\bPOSS-I+\b|\bJ2000\b|\bG < 17\b|\bXE\d{3}\b|\b[1-9]\d?(st|nd|rd|th)\b"


errors, notes = [], []
B = briefs_of((ROOT / "report/briefs.md").read_text())
C = slots_of((ROOT / "report/copy.md").read_text())
N = json.loads((ROOT / "report/numbers.json").read_text())
from registry import CONSTS                                                  # noqa: E402
reg = {r[0] for r in REG} | {c_[0] for c_ in CONSTS} | set(N["numbers"])    # with the values derived in registry.py

# structure
for k in C:
    if k not in B:
        errors.append(f"copy slot {k} has no brief")
undrafted = [k for k in B if k not in C or not C[k]]

# text to measurement
for k, b in B.items():
    for key in b["keys"]:
        if key not in reg:
            errors.append(f"brief {k} names {key}, which is not in the registry")
for k, text in C.items():
    used = set(re.findall(r"\{n:([\w.]+)\}", text))
    for key in used - reg:
        errors.append(f"slot {k} quotes {{n:{key}}}, which is not in the registry")
    for key in used & reg:
        if key not in N["numbers"]:
            errors.append(f"slot {k} quotes {key}, which has no value")
    for key in set(B.get(k, {}).get("keys", [])) - used:
        errors.append(f"slot {k} does not quote {key}, which its brief names")

# figures: every caption's figure is built, and every built figure has a caption
F = json.loads((ROOT / "report/figures.json").read_text())
from figures import FIGS                                                      # noqa: E402
named_figs = {f.strip() for b in B.values() for f in b.get("figure", "").split(",") if f.strip()}
for f in sorted(x for x in named_figs if x.startswith("img:")):              # images: photographs, beside the page
    if not (ROOT / "docs/img" / f[4:]).exists():
        errors.append(f"image {f[4:]} is named by a brief but missing from docs/img")
named_figs = {x for x in named_figs if not x.startswith("img:")}
for f in sorted(named_figs - set(FIGS)):
    errors.append(f"figure {f} is named by a brief but has no builder in report/figures.py")
for f in sorted(set(FIGS) - named_figs):
    errors.append(f"figure {f} is built but no caption's brief names it")
for f in sorted(named_figs & set(FIGS)):
    if f not in F["figures"]:
        errors.append(f"figure {f} has no data")

# measurement to text
named = {key for b in B.values() for key in b["keys"]}
for key in sorted(reg - named):
    errors.append(f"registry key {key} is named by no brief")

# style
for k, text in C.items():
    plain = re.sub(r"\{n:[\w.]+\}|\{fig:[\w.:-]+\}", "", text)
    if "—" in plain or " -- " in plain:
        errors.append(f"slot {k}: em dash")
    for m in re.finditer(AMERICAN, plain, re.I):
        if m.group(0).lower() not in NOT_AMERICAN:
            errors.append(f"slot {k}: American spelling '{m.group(0)}'")
    for m in re.finditer(STEERING, plain, re.I):
        errors.append(f"slot {k}: '{m.group(0)}' (writing-well: steering, inflated, or slop word)")
    rest = plain
    for lit in [x.strip() for x in B.get(k, {}).get("literals", "").split(";") if x.strip()]:   # facts the brief allows as typed
        rest = rest.replace(lit, "")
    rest = re.sub(ALLOWED_DIGITS, "", rest)
    for m in re.finditer(r"\d[\d,.]*", rest):
        notes.append(f"slot {k}: literal number '{m.group(0)}' (quote it from the registry, or whitelist it in check.py)")
    for m in re.finditer(r"\b(\w+), (\w+(?: \w+)?) (and|or) (\w+)", plain):
        notes.append(f"slot {k}: '{m.group(0)}' may lack the Oxford comma")

report = [f"check: {len(B)} briefs, {len(C)} slots drafted, {len(undrafted)} to draft; {len(N['numbers'])} numbers", ""]
report += [f"ERROR {e}" for e in errors] + [f"note  {n}" for n in notes]
if undrafted:
    report.append("to draft: " + ", ".join(undrafted))
(ROOT / "report/check.txt").write_text("\n".join(report) + "\n")
print("\n".join(report[:40]))
sys.exit(1 if errors else 0)
