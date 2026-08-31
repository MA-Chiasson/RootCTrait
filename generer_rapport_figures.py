"""generer_rapport_figures.py
Builds a single index report (rapport_figures.html) to review all the 3D figures
of the pipeline without loading them all at once.

Anti-crash principle: each Plotly figure is loaded (in an iframe) only when its
thumbnail is clicked. Hundreds of samples can therefore be reviewed without
saturating the browser; only the opened figures use memory.

Review features:
  - navigation by batch (tabs);
  - for each sample: Good / Doubtful / Bad buttons, stored in the browser
    (localStorage);
  - free note field per sample;
  - "Export CSV" button to retrieve all the judgments at the end;
  - progress counter per batch.

Placement: put this script at the ROOT (next to the results/ folder) and run it.
It scans results/<batch>/figures/*.html and writes rapport_figures.html there.
"""
import os, glob, json, re

RESULTS_ROOT = "results"
SORTIE = os.path.join(RESULTS_ROOT, "rapport_figures.html")

# Preferred display order of batches (others follow alphabetically)
PREFERRED_ORDER = []  # optional preferred display order of batches


def sample_number(name):
    m = re.search(r"(\d+)", name)
    return int(m.group(1)) if m else 0


def scan_figures():
    """Returns {batch: [(name, relative_path), ...]} sorted by number."""
    batches = {}
    if not os.path.isdir(RESULTS_ROOT):
        return batches
    for batch in os.listdir(RESULTS_ROOT):
        figdir = os.path.join(RESULTS_ROOT, batch, "figures")
        if not os.path.isdir(figdir):
            continue
        figs = []
        for fp in glob.glob(os.path.join(figdir, "*.html")):
            name = os.path.splitext(os.path.basename(fp))[0]
            rel = os.path.join(RESULTS_ROOT, batch, "figures", os.path.basename(fp)).replace("\\", "/")
            figs.append((name, rel))
        figs.sort(key=lambda t: sample_number(t[0]))
        if figs:
            batches[batch] = figs
    return batches


def order_batches(batches):
    ordered = [b for b in PREFERRED_ORDER if b in batches]
    ordered += sorted(b for b in batches if b not in PREFERRED_ORDER)
    return ordered


def build(batches):
    order = order_batches(batches)
    total = sum(len(v) for v in batches.values())

    # donnees injectees en JS
    data = {b: [{"name": n, "src": s} for n, s in batches[b]] for b in order}

    html = []
    html.append("""<!DOCTYPE html><html lang="fr"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Root figure review</title>
<style>
  :root { --good:#2ca02c; --doubtful:#ff7f0e; --bad:#d62728; --bg:#f7f9fc; }
  * { box-sizing: border-box; }
  body { font-family: system-ui, Arial, sans-serif; margin: 0; background: var(--bg); color: #222; }
  header { background: #1f3864; color: #fff; padding: 14px 20px; position: sticky; top: 0; z-index: 10; }
  header h1 { margin: 0 0 6px; font-size: 18px; }
  .barre { display: flex; gap: 8px; flex-wrap: wrap; align-items: center; }
  .onglet { background: #2e5496; color: #fff; border: none; padding: 7px 13px; border-radius: 5px;
            cursor: pointer; font-size: 13px; }
  .onglet.actif { background: #fff; color: #1f3864; font-weight: bold; }
  .onglet .cpt { opacity: .8; font-size: 11px; }
  #recherche { padding: 7px 10px; border: none; border-radius: 5px; font-size: 13px; min-width: 160px; }
  .expbtn { background: #2ca02c; color: #fff; border: none; padding: 7px 13px; border-radius: 5px;
            cursor: pointer; font-size: 13px; margin-left: auto; }
  .grille { display: grid; grid-template-columns: 1fr;
            gap: 14px; padding: 16px; max-width: 1100px; margin: 0 auto; }
  .carte { background: #fff; border-radius: 8px; box-shadow: 0 1px 4px rgba(0,0,0,.08);
           overflow: hidden; display: flex; flex-direction: column; }
  .carte h3 { margin: 0; padding: 8px 12px; font-size: 14px; background: #eef2f9; }
  .cadre { position: relative; height: 640px; background: #fafafa; border-top: 1px solid #eee;
           border-bottom: 1px solid #eee; }
  .cadre iframe { width: 100%; height: 100%; border: 0; }
  .placeholder { position: absolute; inset: 0; display: flex; align-items: center; justify-content: center;
                 cursor: pointer; color: #2e5496; font-size: 14px; flex-direction: column; gap: 6px; }
  .placeholder:hover { background: #eef4ff; }
  .judgments { display: flex; gap: 6px; padding: 8px 12px; }
  .jbtn { flex: 1; border: 1px solid #ccc; background: #fff; padding: 6px; border-radius: 5px;
          cursor: pointer; font-size: 12px; }
  .jbtn.good.on { background: var(--good); color: #fff; border-color: var(--good); }
  .jbtn.doubtful.on { background: var(--doubtful); color: #fff; border-color: var(--doubtful); }
  .jbtn.bad.on { background: var(--bad); color: #fff; border-color: var(--bad); }
  .note { width: 100%; border: 1px solid #eee; border-top: 0; padding: 6px 12px; font-size: 12px; }
  .carte[data-jug="good"] { outline: 3px solid var(--good); }
  .carte[data-jug="doubtful"] { outline: 3px solid var(--doubtful); }
  .carte[data-jug="bad"] { outline: 3px solid var(--bad); }
  footer { padding: 20px; text-align: center; color: #888; font-size: 12px; }
</style></head><body>""")

    html.append(f"""<header>
  <h1>Root figure review &mdash; {total} samples</h1>
  <div class="barre" id="onglets"></div>
  <div class="barre" style="margin-top:8px">
    <input id="recherche" placeholder="filter (number or name)...">
    <button class="expbtn" onclick="exportCsv()">Export CSV</button>
  </div>
</header>
<div id="contenu" class="grille"></div>
<footer>Judgments are stored in this browser. Export to CSV to keep them.</footer>
<script>
const DATA = {json.dumps(data, ensure_ascii=False)};
const ORDER = {json.dumps(order)};
let activeBatch = ORDER[0];

function makeKey(batch, name) {{ return "jug:" + batch + ":" + name; }}
function noteKey(batch, name) {{ return "note:" + batch + ":" + name; }}

function counter(batch) {{
  let n = 0;
  for (const it of DATA[batch]) if (localStorage.getItem(makeKey(batch, it.name))) n++;
  return n;
}}

function renderTabs() {{
  const box = document.getElementById("onglets");
  box.innerHTML = "";
  for (const b of ORDER) {{
    const btn = document.createElement("button");
    btn.className = "onglet" + (b === activeBatch ? " actif" : "");
    btn.innerHTML = b + ' <span class="cpt">(' + counter(b) + '/' + DATA[b].length + ')</span>';
    btn.onclick = () => {{ activeBatch = b; render(); }};
    box.appendChild(btn);
  }}
}}

function render() {{
  renderTabs();
  const c = document.getElementById("contenu");
  c.innerHTML = "";
  const filtre = document.getElementById("recherche").value.toLowerCase();
  for (const it of DATA[activeBatch]) {{
    if (filtre && !it.name.toLowerCase().includes(filtre)) continue;
    const jug = localStorage.getItem(makeKey(activeBatch, it.name)) || "";
    const note = localStorage.getItem(noteKey(activeBatch, it.name)) || "";
    const carte = document.createElement("div");
    carte.className = "carte";
    if (jug) carte.setAttribute("data-jug", jug);
    carte.innerHTML =
      '<h3>' + it.name + '</h3>' +
      '<div class="cadre"><div class="placeholder" onclick="loadFig(this,\\'' + it.src + '\\')">' +
        '<div>&#128065; click to display</div><div style="font-size:11px;opacity:.7">3D figure</div>' +
      '</div></div>' +
      '<div class="judgments">' +
        '<button class="jbtn good' + (jug==="good"?" on":"") + '" onclick="judge(this,\\'' + it.name + '\\',\\'good\\')">Good</button>' +
        '<button class="jbtn doubtful' + (jug==="doubtful"?" on":"") + '" onclick="judge(this,\\'' + it.name + '\\',\\'doubtful\\')">Doubtful</button>' +
        '<button class="jbtn bad' + (jug==="bad"?" on":"") + '" onclick="judge(this,\\'' + it.name + '\\',\\'bad\\')">Bad</button>' +
      '</div>' +
      '<input class="note" placeholder="note..." value="' + note.replace(/"/g,"&quot;") + '" ' +
        'oninput="setNote(\\'' + it.name + '\\', this.value)">';
    c.appendChild(carte);
  }}
}}

function loadFig(el, src) {{
  const f = document.createElement("iframe");
  f.src = src; f.loading = "lazy";
  el.parentNode.appendChild(f); el.remove();
}}

function judge(btn, name, val) {{
  const k = makeKey(activeBatch, name);
  const actuel = localStorage.getItem(k);
  const carte = btn.closest(".carte");
  carte.querySelectorAll(".jbtn").forEach(b => b.classList.remove("on"));
  if (actuel === val) {{
    localStorage.removeItem(k); carte.removeAttribute("data-jug");
  }} else {{
    localStorage.setItem(k, val); btn.classList.add("on"); carte.setAttribute("data-jug", val);
  }}
  renderTabs();
}}

function setNote(name, val) {{
  const k = noteKey(activeBatch, name);
  if (val) localStorage.setItem(k, val); else localStorage.removeItem(k);
}}

function exportCsv() {{
  let rows = [["batch","sample","judgment","note"]];
  for (const b of ORDER) for (const it of DATA[b]) {{
    const j = localStorage.getItem(makeKey(b, it.name)) || "";
    const n = localStorage.getItem(noteKey(b, it.name)) || "";
    if (j || n) rows.push([b, it.name, j, '"' + n.replace(/"/g,'""') + '"']);
  }}
  const csv = rows.map(l => l.join(",")).join("\\n");
  const blob = new Blob([csv], {{type:"text/csv"}});
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob); a.download = "figure_judgments.csv"; a.click();
}}

document.getElementById("recherche").addEventListener("input", render);
render();
</script></body></html>""")
    return "".join(html)


def main():
    batches = scan_figures()
    if not batches:
        print(f"No figure found in {RESULTS_ROOT}/<batch>/figures/.")
        print("Check that the pipeline ran with SAVE_FIGURES=1 and that this script")
        print("is placed at the root (next to the results/ folder).")
        return
    html = build(batches)
    with open(SORTIE, "w", encoding="utf-8") as f:
        f.write(html)
    total = sum(len(v) for v in batches.values())
    print(f"Report written: {SORTIE}  ({total} samples, {len(batches)} batches)")
    for b in order_batches(batches):
        print(f"  {b}: {len(batches[b])} figures")
    print("\nOpen rapport_figures.html in a browser.")
    print("Figures load on click (nothing is preloaded, no crash).")


if __name__ == "__main__":
    main()
