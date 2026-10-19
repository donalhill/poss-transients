"""Local editor for the page: edit text in place, and pin comments to text, charts, images, or tables.

    python site/edit.py            then open http://localhost:8766

Text edits are written straight into site/copy.md (the previous version is kept in logs/copy_history/) and the public
page docs/index.html is rebuilt. Comments go to site/notes.json, mirrored as site/notes.md for reading."""
import json
import re
import sys
from datetime import datetime
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from assemble import ROOT, main as rebuild, parse, render  # noqa: E402

COPY = ROOT / "site/copy.md"
NOTES = ROOT / "site/notes.json"
HISTORY = ROOT / "logs/copy_history"
PORT = 8766


def load_notes():
    return json.loads(NOTES.read_text()) if NOTES.exists() else []


def save_notes(notes):
    NOTES.write_text(json.dumps(notes, indent=1, ensure_ascii=False))
    lines = ["# Notes on the page", "", "Open notes first. Each is pinned to the element named after the @.", ""]
    for n in sorted(notes, key=lambda n: (n["status"] != "open", n["id"])):
        where = n["anchor"] + (f' at {n["x"]:.0%} across, {n["y"]:.0%} down' if n.get("x") is not None else "")
        lines.append(f'## {n["id"]}. [{n["status"]}] @ {where}')
        if n.get("snippet"):
            lines.append(f'> {n["snippet"]}')
        lines += [n["text"], f'_{n["time"]}_', ""]
    (ROOT / "site/notes.md").write_text("\n".join(lines))


def set_slot(slot, text):
    src = COPY.read_text()
    slots = parse(src)
    if slot not in slots:
        raise KeyError(slot)
    HISTORY.mkdir(parents=True, exist_ok=True)
    (HISTORY / f"copy_{datetime.now():%Y%m%d_%H%M%S}.md").write_text(src)
    pat = re.compile(r"(^### " + re.escape(slot) + r"[ \t]*\n)(.*?)(?=^### |\Z)", re.S | re.M)
    new = pat.sub(lambda m: m.group(1) + text.strip() + "\n\n", src, count=1)
    COPY.write_text(new)
    (HISTORY / f"copy_{datetime.now():%Y%m%d_%H%M%S}_saved.md").write_text(new)   # the version just saved, too


class Handler(SimpleHTTPRequestHandler):
    def __init__(self, *a, **k):
        super().__init__(*a, directory=str(ROOT / "docs"), **k)

    def log_message(self, *a):
        pass

    def send_json(self, obj, code=200):
        body = json.dumps(obj).encode()
        self.send_response(code)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        if self.path.split("?")[0] in ("/", "/index.html"):
            html, missing, _ = render(COPY, edit=True)
            html = html.replace("</body>", EDITOR + "\n</body>")
            body = html.encode()
            self.send_response(200)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Cache-Control", "no-store")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)
        elif self.path == "/api/notes":
            self.send_json(load_notes())
        else:
            super().do_GET()

    def end_headers(self):
        self.send_header("Cache-Control", "no-store")
        super().end_headers()

    def do_POST(self):
        data = json.loads(self.rfile.read(int(self.headers.get("Content-Length", 0))) or b"{}")
        try:
            if self.path == "/api/save":
                for slot, text in data.items():
                    set_slot(slot, text)
                rebuild(COPY)
                self.send_json({"ok": True, "saved": list(data)})
            elif self.path == "/api/notes":
                notes = load_notes()
                note = {k: data.get(k) for k in ("anchor", "snippet", "x", "y", "text")}
                note.update(id=max([n["id"] for n in notes] + [0]) + 1, status="open", time=f"{datetime.now():%Y-%m-%d %H:%M}")
                notes.append(note)
                save_notes(notes)
                self.send_json(note)
            elif self.path == "/api/notes/status":
                notes = load_notes()
                for n in notes:
                    if n["id"] == data["id"]:
                        n["status"] = data["status"]
                save_notes(notes)
                self.send_json({"ok": True})
            else:
                self.send_json({"error": "unknown"}, 404)
        except Exception as e:  # report problems to the page instead of failing silently
            self.send_json({"error": str(e)}, 400)


EDITOR = r"""
<style>
  .slot { outline: 1px dashed transparent; outline-offset: 3px; border-radius: 2px; transition: outline-color .15s; }
  body.ed-text .slot:hover { outline-color: #2462a8; cursor: text; }
  body.ed-text .slot[contenteditable="true"]:focus { outline: 2px solid #2462a8; background: rgba(36,98,168,.04); }
  .slot.dirty { outline: 2px solid #e0a020 !important; }
  body.ed-text .num { background: rgba(16,128,107,.12); border-radius: 2px; cursor: help; }
  body.ed-note, body.ed-note * { cursor: crosshair !important; }
  body.ed-note .ed-hover { outline: 2px solid #b8501a !important; outline-offset: 2px; }
  #ed-bar { position: fixed; right: 16px; bottom: 16px; z-index: 1000; display: flex; gap: 6px; align-items: center;
    background: #13171b; color: #fff; padding: 8px; border-radius: 10px; font: 13px/1.2 "IBM Plex Sans", sans-serif;
    box-shadow: 0 6px 24px rgba(0,0,0,.25); }
  #ed-bar button { font: inherit; color: #fff; background: #2a3138; border: 0; border-radius: 6px; padding: 7px 11px; cursor: pointer; }
  #ed-bar button[aria-pressed="true"] { background: #2462a8; }
  #ed-bar #ed-save.ready { background: #e0a020; color: #13171b; }
  #ed-status { padding: 0 6px; color: #b9c2ca; min-width: 6rem; }
  .ed-pin { position: absolute; z-index: 999; width: 24px; height: 24px; margin: -12px 0 0 -12px; border-radius: 50%;
    background: #b8501a; color: #fff; font: 600 12px/24px "IBM Plex Sans", sans-serif; text-align: center;
    box-shadow: 0 0 0 2px #fff, 0 2px 8px rgba(0,0,0,.3); cursor: pointer; }
  .ed-pin.resolved { background: #8a949c; opacity: .6; }
  #ed-pop { position: absolute; z-index: 1001; width: 320px; background: #fff; color: #13171b; border-radius: 10px;
    box-shadow: 0 8px 30px rgba(0,0,0,.25); padding: 12px; font: 14px/1.4 "IBM Plex Sans", sans-serif; }
  #ed-pop .where { font-size: 12px; color: #5b656e; margin-bottom: 6px; }
  #ed-pop textarea { width: 100%; box-sizing: border-box; min-height: 90px; font: inherit; padding: 8px; border: 1px solid #cfd6dc; border-radius: 6px; }
  #ed-pop .row { display: flex; gap: 6px; justify-content: flex-end; margin-top: 8px; }
  #ed-pop button { font: inherit; border: 0; border-radius: 6px; padding: 6px 12px; cursor: pointer; background: #e7ebee; }
  #ed-pop button.primary { background: #2462a8; color: #fff; }
  #ed-list { position: fixed; right: 16px; bottom: 70px; z-index: 1000; width: 360px; max-height: 60vh; overflow: auto;
    background: #fff; border-radius: 10px; box-shadow: 0 8px 30px rgba(0,0,0,.25); padding: 8px 12px; font: 13px/1.4 "IBM Plex Sans", sans-serif; color: #13171b; }
  #ed-list .n { border-bottom: 1px solid #e7ebee; padding: 8px 0; }
  #ed-list .n.resolved { opacity: .5; }
  #ed-list .a { color: #5b656e; font-size: 12px; }
  #ed-list button { font: inherit; font-size: 12px; border: 0; background: #e7ebee; border-radius: 5px; padding: 3px 8px; cursor: pointer; margin-top: 4px; }
</style>
<div id="ed-bar" role="toolbar" aria-label="Editor">
  <button id="ed-text" aria-pressed="true" title="Click any text to edit it">Edit text</button>
  <button id="ed-note" aria-pressed="false" title="Click anything to pin a comment to it">Comment</button>
  <button id="ed-save" title="Save text edits (⌘S)">Save</button>
  <button id="ed-notes" title="Show comments">Comments</button>
  <span id="ed-status">Ready</span>
</div>
<script>
(() => {
  const $ = s => document.querySelector(s), $$ = s => [...document.querySelectorAll(s)];
  const status = t => { $("#ed-status").textContent = t; };
  const dirty = new Set();
  let mode = "text", notes = [];

  // ---------- text editing ----------
  const clean = html => html.replace(/<span[^>]*data-step="(\d+)"[^>]*>.*?<\/span>/g, "{$1_circle_symbol}").replace(/<span[^>]*data-num="([\w.-]+)"[^>]*>.*?<\/span>/g, "{n:$1}").replace(/<(?!\/?em>)[^>]+>/g, "").replace(/&nbsp;/g, " ").replace(/&amp;/g, "&")
    .replace(/&lt;/g, "<").replace(/&gt;/g, ">").replace(/\s+/g, " ").trim();
  function slotText(el) {
    const m = el.dataset.mode;
    if (m === "inline") return clean(el.innerHTML);
    if (m === "headline") return [clean((el.querySelector(".big") || el).innerHTML), clean((el.querySelector(".gloss") || {innerHTML: ""}).innerHTML)].filter(Boolean).join("\n");
    const out = [];
    const kids = [...el.children].flatMap(c => c.classList.contains("keep") ? [...c.children] : [c]);
    for (const c of kids) {
      if (c.tagName === "UL") out.push([...c.children].map(li => "- " + clean(li.innerHTML)).join("\n"));
      else { const t = clean(c.innerHTML); if (t) out.push(t); }
    }
    if (!el.children.length) out.push(clean(el.innerHTML));
    return out.join("\n\n");
  }
  function setMode(m) {
    mode = m;
    document.body.classList.toggle("ed-text", m === "text");
    document.body.classList.toggle("ed-note", m === "note");
    $("#ed-text").setAttribute("aria-pressed", m === "text");
    $("#ed-note").setAttribute("aria-pressed", m === "note");
    $$(".slot").forEach(s => s.contentEditable = m === "text" ? "true" : "false");
    status(m === "text" ? "Click text to edit" : "Click anything to comment");
  }
  document.addEventListener("input", e => {
    const s = e.target.closest && e.target.closest(".slot");
    if (s) { dirty.add(s.dataset.slot); s.classList.add("dirty"); $("#ed-save").classList.add("ready"); status(dirty.size + " unsaved"); }
  });
  async function save() {
    if (!dirty.size) { status("Nothing to save"); return; }
    const body = {};
    for (const k of dirty) { const el = $(`.slot[data-slot="${k}"]`); body[k] = slotText(el); }
    status("Saving…");
    const r = await fetch("/api/save", { method: "POST", body: JSON.stringify(body) }).then(r => r.json());
    if (r.error) { status("Save failed: " + r.error); return; }
    dirty.clear(); $$(".slot.dirty").forEach(s => s.classList.remove("dirty")); $("#ed-save").classList.remove("ready");
    status("Saved " + r.saved.length);
  }
  $("#ed-save").onclick = save;
  document.addEventListener("keydown", e => { if ((e.metaKey || e.ctrlKey) && e.key === "s") { e.preventDefault(); save(); } });
  window.addEventListener("beforeunload", e => { if (dirty.size) { e.preventDefault(); e.returnValue = ""; } });
  $("#ed-text").onclick = () => setMode("text");
  $("#ed-note").onclick = () => setMode(mode === "note" ? "text" : "note");

  // ---------- comments ----------
  function anchorOf(el) {
    const slot = el.closest(".slot");
    if (slot) return { node: slot, anchor: "text " + slot.dataset.slot, snippet: slot.textContent.trim().slice(0, 140) };
    const svgEl = el.closest("svg[id]");
    if (svgEl) return { node: svgEl, anchor: "chart " + svgEl.id, snippet: (svgEl.getAttribute("aria-label") || "").slice(0, 140) };
    const img = el.closest("img");
    if (img) return { node: img, anchor: "image " + img.getAttribute("src"), snippet: img.alt.slice(0, 140) };
    const tbl = el.closest("table");
    if (tbl) { const i = $$("table").indexOf(tbl); return { node: tbl, anchor: "table " + (tbl.id || "#" + i), snippet: (tbl.querySelector("th") || tbl).textContent.trim().slice(0, 80) }; }
    const blk = el.closest(".legend, .cap, .chart, .claim, .score, .verdict, .side, section");
    if (blk) { const sec = blk.closest("section"); return { node: blk, anchor: (blk.className || "block").split(" ")[0] + " in " + (sec ? "#" + sec.id : "page"), snippet: blk.textContent.trim().slice(0, 140) }; }
    return null;
  }
  function locate(n) {   // find the anchored element again, to place its pin
    const [kind, ...rest] = n.anchor.split(" "), key = rest.join(" ");
    if (kind === "text") return $(`.slot[data-slot="${key}"]`);
    if (kind === "chart") return document.getElementById(key);
    if (kind === "image") return $(`img[src="${key}"]`);
    if (kind === "table") return key.startsWith("#") && key.length > 1 && isNaN(key.slice(1)) ? $(key) : $$("table")[+key.slice(1)];
    const m = key.match(/#([\w-]+)/); const sec = m ? document.getElementById(m[1]) : document;
    return sec && sec.querySelector("." + kind);
  }
  function placePins() {
    $$(".ed-pin").forEach(p => p.remove());
    notes.forEach(n => {
      const el = locate(n); if (!el) return;
      const r = el.getBoundingClientRect(), p = document.createElement("div");
      p.className = "ed-pin" + (n.status === "open" ? "" : " resolved"); p.textContent = n.id; p.title = n.text;
      p.style.left = (scrollX + r.left + r.width * (n.x ?? 1)) + "px"; p.style.top = (scrollY + r.top + r.height * (n.y ?? 0)) + "px";
      p.onclick = e => { e.stopPropagation(); showList(n.id); };
      document.body.appendChild(p);
    });
  }
  let hoverEl = null;
  document.addEventListener("mousemove", e => {
    if (mode !== "note") return;
    const a = anchorOf(e.target); const node = a && a.node;
    if (hoverEl !== node) { hoverEl && hoverEl.classList.remove("ed-hover"); hoverEl = node; node && node.classList.add("ed-hover"); }
  });
  document.addEventListener("click", e => {
    if (mode !== "note" || e.target.closest("#ed-bar, #ed-pop, #ed-list, .ed-pin")) return;
    e.preventDefault(); e.stopPropagation();
    const a = anchorOf(e.target); if (!a) return;
    const r = a.node.getBoundingClientRect();
    openPop(a, (e.clientX - r.left) / r.width, (e.clientY - r.top) / r.height, e.pageX, e.pageY);
  }, true);
  function openPop(a, x, y, px, py) {
    $("#ed-pop") && $("#ed-pop").remove();
    const pop = document.createElement("div"); pop.id = "ed-pop";
    pop.innerHTML = `<div class="where">On ${a.anchor}</div><textarea placeholder="Your comment for Claude"></textarea>
      <div class="row"><button data-x>Cancel</button><button class="primary" data-ok>Add comment</button></div>`;
    pop.style.left = Math.min(px + 8, scrollX + innerWidth - 340) + "px"; pop.style.top = (py + 8) + "px";
    document.body.appendChild(pop); pop.querySelector("textarea").focus();
    pop.querySelector("[data-x]").onclick = () => pop.remove();
    pop.querySelector("[data-ok]").onclick = async () => {
      const text = pop.querySelector("textarea").value.trim(); if (!text) return;
      const n = await fetch("/api/notes", { method: "POST", body: JSON.stringify({ anchor: a.anchor, snippet: a.snippet, x, y, text }) }).then(r => r.json());
      notes.push(n); pop.remove(); placePins(); status("Comment " + n.id + " added");
    };
  }
  function showList(focus) {
    let box = $("#ed-list");
    if (box && focus === undefined) { box.remove(); return; }
    box && box.remove(); box = document.createElement("div"); box.id = "ed-list";
    box.innerHTML = notes.length ? "" : "<p>No comments yet. Click Comment, then click anything on the page.</p>";
    [...notes].sort((a, b) => (a.status !== "open") - (b.status !== "open") || a.id - b.id).forEach(n => {
      const d = document.createElement("div"); d.className = "n" + (n.status === "open" ? "" : " resolved");
      d.innerHTML = `<b>${n.id}.</b> ${n.text.replace(/</g, "&lt;")}<div class="a">${n.anchor} · ${n.status}</div>`;
      const go = document.createElement("button"); go.textContent = "Show"; go.onclick = () => { const el = locate(n); el && el.scrollIntoView({ block: "center" }); setTimeout(placePins, 400); };
      const st = document.createElement("button"); st.textContent = n.status === "open" ? "Resolve" : "Reopen";
      st.onclick = async () => { n.status = n.status === "open" ? "resolved" : "open"; await fetch("/api/notes/status", { method: "POST", body: JSON.stringify({ id: n.id, status: n.status }) }); placePins(); showList(n.id); };
      d.append(go, " ", st); box.appendChild(d);
      if (n.id === focus) setTimeout(() => d.scrollIntoView({ block: "nearest" }), 0);
    });
    document.body.appendChild(box);
  }
  $("#ed-notes").onclick = () => showList();
  addEventListener("resize", placePins);
  fetch("/api/notes").then(r => r.json()).then(n => { notes = n; setTimeout(placePins, 800); });
  setMode("text");
})();
</script>
"""

if __name__ == "__main__":
    srv = ThreadingHTTPServer(("127.0.0.1", PORT), Handler)
    print(f"Editor running at http://localhost:{PORT}  (edits go to site/copy.md, comments to site/notes.md)")
    srv.serve_forever()
