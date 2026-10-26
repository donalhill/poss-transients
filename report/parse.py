"""Parsers shared by the report's scripts: the copy's slots and the briefs."""
import re


def slots_of(text):
    """'### id' blocks: {id: text}."""
    out, cur, buf = {}, None, []
    for line in text.splitlines():
        m = re.match(r"^###\s+([\w-]+)\s*$", line)
        if m:
            if cur:
                out[cur] = "\n".join(buf).strip()
            cur, buf = m.group(1), []
        elif cur is not None:
            buf.append(line)
    if cur:
        out[cur] = "\n".join(buf).strip()
    return out


def briefs_of(text):
    """Briefs by slot, each with its section and fields."""
    out, sec = {}, None
    for block in re.split(r"(?m)^(?=## |### )", text):
        h = block.splitlines()[0] if block.strip() else ""
        if h.startswith("## "):
            m = re.match(r"## (\d+)\.", h); sec = int(m.group(1)) if m else sec
        elif h.startswith("### "):
            f = {}
            for name, val in re.findall(r"(?m)^- (\w+): (.*)$", block):          # a field given twice is merged
                f[name] = f"{f[name]}{';' if name == 'literals' else ','} {val}" if name in f else val
            f["keys"] = [k.strip() for k in f.get("keys", "").split(",") if k.strip()]
            f["section"] = sec; out[h[4:].strip()] = f
    return out
