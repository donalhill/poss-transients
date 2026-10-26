# Slot briefs

One entry per text slot of the page, in reading order. report/check.py holds the copy to it: every slot has a brief,
every key a brief names is quoted in that slot, and every registry key appears in some brief.
Fields:
- kind: h2 (a section's question), title (a step or finding), text, caption, or verdict (a result sentence).
- keys: the registry keys it quotes (report/registry.py); the copy writes them as {n:key}.
- literals: numbers and names the copy may write as typed.
- figure: for a caption, the figure ids it describes (report/figures.py), drawn from report/figures.json.
- alt: for a photograph, its alternative text.

