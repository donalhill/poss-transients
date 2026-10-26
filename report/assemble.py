"""Fill the page's template (report/template.py) with the copy (report/copy.md, from report/pages.py) and the numbers (report/registry.py).

The copy's format: blocks starting '### slot-id', followed by text (paragraphs separated by blank lines, bullets as '- ').
Template placeholders: {{slot}} (block: paragraphs and lists), {{slot|inline}} (one run of text),
{{slot|headline}} (first line as the big figure, the rest as a gloss).
Numbers in the copy are written as {n:key} and filled from the numbers file's "numbers" table, which the pipeline
writes; a key with no value stops the build.
usage: python report/assemble.py COPY TEMPLATE NUMBERS OUT"""
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def parse(text):
    slots, cur, buf = {}, None, []
    for line in text.splitlines():
        m = re.match(r"^###\s+([\w-]+)\s*$", line)
        if m:
            if cur:
                slots[cur] = "\n".join(buf).strip()
            cur, buf = m.group(1), []
        elif cur:
            buf.append(line)
    if cur:
        slots[cur] = "\n".join(buf).strip()
    return slots


def esc(s):
    return re.sub(r"&(?!(?:[a-zA-Z]+|#\d+);)", "&amp;", s)


def block(s):
    out, para, items = [], [], []
    def flush():
        if para:
            out.append("<p>" + " ".join(para) + "</p>"); para.clear()
        if items:
            out.append("<ul>" + "".join(f"<li>{i}</li>" for i in items) + "</ul>"); items.clear()
    for line in s.splitlines():
        t = line.strip()
        if not t:
            flush()
        elif t.startswith("- "):
            if para:
                out.append("<p>" + " ".join(para) + "</p>"); para.clear()
            items.append(t[2:].strip())
        else:
            if items:
                out.append("<ul>" + "".join(f"<li>{i}</li>" for i in items) + "</ul>"); items.clear()
            para.append(t)
    flush()
    # keep a lead-in sentence ending in a colon together with the list that follows, so a column break cannot split them
    return re.sub(r"(<p>[^<]*:</p>)(<ul>.*?</ul>)", r'<div class="keep">\1\2</div>', "".join(out))


def inline(s):
    return " ".join(t.strip() for t in s.splitlines() if t.strip())


def headline(s):
    lines = [t.strip() for t in s.splitlines() if t.strip()]
    big, gloss = lines[0], " ".join(lines[1:])
    return f'<div class="big">{big}</div>' + (f'<div class="gloss">{gloss}</div>' if gloss else "")


EQS = {"track-score": '<i>μ</i> = <span class="frac"><span>π <i>N</i> <i>n</i> (2<i>n</i><i>w</i>)<sup><i>k</i>−2</sup> <i>L</i><sup><i>k</i></sup></span>'
                        '<span><i>k</i> (<i>k</i> − 2)!</span></span>'}


def numbers(path):
    """The numbers table ({"numbers": ...})."""
    return json.loads(Path(path).read_text())["numbers"]


def render(copy_path, template, numbers_path, edit=False):
    """Fill the template. With edit=True every slot is wrapped in an element tagged with its id, and every number in a
    locked span tagged with its key, for the editor (report/edit.py)."""
    slots = {k: esc(v) for k, v in parse(Path(copy_path).read_text()).items()}
    nums, unknown = numbers(numbers_path), []

    def num(m):
        if m.group(1) not in nums:
            unknown.append(m.group(1))
            return f"[[{m.group(1)}]]"
        v = esc(str(nums[m.group(1)]))
        return f'<span class="num" data-num="{m.group(1)}" title="From the pipeline: {m.group(1)}" contenteditable="false">{v}</span>' if edit else v
    tpl = Path(template).read_text()
    missing = []

    def sub(m):
        key, mode = m.group(1), m.group(2) or "block"
        if key not in slots:
            missing.append(key)
            return f"[[{key}]]"
        out = {"inline": inline, "headline": headline, "block": block}[mode](slots[key])
        # {2_circle_symbol} in the copy becomes an inline copy of the numbered step marker
        out = re.sub(r"\{(\d+)_circle_symbol\}", r'<span class="step-ref" data-step="\1" contenteditable="false">\1</span>', out)
        out = re.sub(r"\{n:([\w.-]+)\}", num, out)
        # {eq:name} on a line of its own becomes that typeset equation, centred
        out = re.sub(r"<p>\{eq:([\w-]+)\}</p>", lambda m: f'<p class="eq-line"><span class="eq" data-eq="{m.group(1)}" contenteditable="false">{EQS[m.group(1)]}</span></p>', out)
        if edit:
            tag = "span" if mode == "inline" else "div"
            out = f'<{tag} class="slot" data-slot="{key}" data-mode="{mode}">{out}</{tag}>'
        return out

    html = re.sub(r"\{\{([\w-]+)(?:\|(inline|headline))?\}\}", sub, tpl)
    html = re.sub(r"\{n:([\w.-]+)\}", lambda m: re.sub(r"<[^>]+>", "", num(m)), html)   # numbers in the template itself, always plain
    # exactly one standard deviation reads as such: "1.0 standard deviations" becomes "1 standard deviation"
    if not edit:                                                                  # the editor shows the copy as written
        html = re.sub(r"\b1\.0((?:</span>)?)( |&nbsp;)standard deviations", r"1\1\2standard deviation", html)
    unused = sorted(set(slots) - set(re.findall(r"\{\{([\w-]+)", tpl)))
    if unknown:
        raise SystemExit(f"numbers with no value in {numbers_path}: {sorted(set(unknown))}")
    return html, missing, unused


def main(copy_path, template, numbers_path, out):
    html, missing, unused = render(copy_path, template, numbers_path)
    out = Path(out); out.parent.mkdir(parents=True, exist_ok=True); out.write_text(html)
    print(f"wrote {out}; missing slots: {missing or 'none'}; unused slots: {unused or 'none'}")
    return missing


if __name__ == "__main__":
    main(*sys.argv[1:5])
