"""Fill site/template.html with the copy in site/copy.md and write docs/index.html.

copy.md format: blocks starting '### slot-id', followed by text (paragraphs separated by blank lines, bullets as '- ').
Template placeholders: {{slot}} (block: paragraphs and lists), {{slot|inline}} (one run of text),
{{slot|headline}} (first line as the big figure, the rest as a gloss)."""
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


def render(copy_path, edit=False):
    """Fill the template. With edit=True every slot is wrapped in an element tagged with its id, for site/edit.py."""
    slots = {k: esc(v) for k, v in parse(Path(copy_path).read_text()).items()}
    tpl = (ROOT / "site/template.html").read_text()
    missing = []

    def sub(m):
        key, mode = m.group(1), m.group(2) or "block"
        if key not in slots:
            missing.append(key)
            return f"[[{key}]]"
        out = {"inline": inline, "headline": headline, "block": block}[mode](slots[key])
        # {2_circle_symbol} in the copy becomes an inline copy of the numbered step marker
        out = re.sub(r"\{(\d+)_circle_symbol\}", r'<span class="step-ref" data-step="\1" contenteditable="false">\1</span>', out)
        if edit:
            tag = "span" if mode == "inline" else "div"
            out = f'<{tag} class="slot" data-slot="{key}" data-mode="{mode}">{out}</{tag}>'
        return out

    html = re.sub(r"\{\{([\w-]+)(?:\|(inline|headline))?\}\}", sub, tpl)
    unused = sorted(set(slots) - set(re.findall(r"\{\{([\w-]+)", tpl)))
    return html, missing, unused


def main(copy_path):
    html, missing, unused = render(copy_path)
    (ROOT / "docs/index.html").write_text(html)
    print(f"wrote docs/index.html; missing slots: {missing or 'none'}; unused slots: {unused or 'none'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ROOT / "site/copy.md")
