"""plant_report.py
Builds a figure report grouped by PLANT (results/plant_report.html): the figures
of the same plant at the different times are shown together, e.g. B2T2S1 and
B2T3S1 side by side under "B2 S1".

How plants are matched: a batch name such as B2T3 is read as block B2 and time
T3 (the T<number> part is the time). The sample name, without its batch prefix
(B2T3S1 -> S1), identifies the plant within its block. Batches without a
T<number> are shown as their own block with a single time.

Layout: one figure takes the full width, two are side by side, three or more
are shown two per row (three in a row is too small). A button stacks all the
figures one below the other.

Each figure keeps its own Good / Doubtful / Bad buttons and note. Judgments use
the same browser storage keys as figure_report.html, so a figure judged in one
report appears judged in the other (same browser, same results folder).

Run from the project root: python -m tools.plant_report
It scans results/<batch>/figures/*.html and writes results/plant_report.html.
"""
import os, re, json, glob
import os as _os, sys as _sys

PROJECT_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
if PROJECT_ROOT not in _sys.path:
    _sys.path.insert(0, PROJECT_ROOT)

from rootctrait.sorting import natural_key
from tools import figure_report as fr

RESULTS_ROOT = fr.RESULTS_ROOT
OUTPUT_HTML = os.path.join(RESULTS_ROOT, "plant_report.html")

TIME_RE = re.compile(r"T\d+", re.IGNORECASE)


def split_batch(batch):
    """B2T3 -> ('B2', 'T3'); a batch without T<number> -> (batch, '')."""
    m = list(TIME_RE.finditer(batch))
    if not m:
        return batch, ""
    t = m[-1]
    return (batch[:t.start()] + batch[t.end():]) or batch, t.group(0).upper()


def plant_id(batch, name):
    """Sample name without its batch prefix: B2T3S1 in batch B2T3 -> S1."""
    return name[len(batch):] if name.startswith(batch) and len(name) > len(batch) else name


def group_by_plant(batches, qc):
    """{block: {plant: [figure, ...]}} with the figures in time order."""
    blocks = {}
    for batch, figs in batches.items():
        block, time = split_batch(batch)
        rk = qc.get(batch, {})
        for name, src in figs:
            fig = {"batch": batch, "name": name, "time": time or batch, "src": src,
                   "flags": rk.get(name, (0, ""))[1]}
            blocks.setdefault(block, {}).setdefault(plant_id(batch, name), []).append(fig)
    out = {}
    for block in sorted(blocks, key=natural_key):
        plants = blocks[block]
        out[block] = [{"plant": p, "figs": sorted(plants[p], key=lambda f: natural_key(f["time"]))}
                      for p in sorted(plants, key=natural_key)]
    return out


PAGE = r"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Root figures by plant</title>
<style>
  :root { --good:#009E73; --doubtful:#E69F00; --bad:#D55E00; --bg:#f7f9fc; }
  * { box-sizing: border-box; }
  body { font-family: system-ui, Arial, sans-serif; margin: 0; background: var(--bg); color: #222; }
  header { background: #1f3864; color: #fff; padding: 14px 20px; position: sticky; top: 0; z-index: 10; }
  header h1 { margin: 0 0 6px; font-size: 18px; }
  .bar { display: flex; gap: 8px; flex-wrap: wrap; align-items: center; }
  .tab, .tool { background: #2e5496; color: #fff; border: none; padding: 7px 13px; border-radius: 5px;
                cursor: pointer; font-size: 13px; }
  .tab.active { background: #fff; color: #1f3864; font-weight: bold; }
  .tab .cpt { opacity: .8; font-size: 11px; }
  #search { padding: 7px 10px; border: none; border-radius: 5px; font-size: 13px; min-width: 160px; }
  .export { background: #009E73; color: #fff; border: none; padding: 7px 13px; border-radius: 5px;
            cursor: pointer; font-size: 13px; margin-left: auto; }
  #content { padding: 16px; max-width: 1900px; margin: 0 auto; display: flex; flex-direction: column; gap: 18px; }
  .plant { background: #fff; border-radius: 8px; box-shadow: 0 1px 4px rgba(0,0,0,.08); overflow: hidden; }
  .plant > h2 { margin: 0; padding: 9px 14px; font-size: 15px; background: #e3e9f4; display: flex;
                align-items: center; gap: 12px; }
  .plant > h2 .count { font-weight: normal; font-size: 12px; color: #555; }
  .plant > h2 button { margin-left: auto; font-size: 12px; padding: 4px 10px; border: 1px solid #2e5496;
                       background: #fff; color: #2e5496; border-radius: 4px; cursor: pointer; }
  .figs { display: grid; gap: 10px; padding: 10px; grid-template-columns: 1fr; }
  .figs.n2, .figs.n3plus { grid-template-columns: 1fr 1fr; }
  body.stacked .figs { grid-template-columns: 1fr !important; }
  .card { border: 1px solid #e5e9f0; border-radius: 6px; overflow: hidden; display: flex; flex-direction: column; }
  .card h3 { margin: 0; padding: 7px 12px; font-size: 13px; background: #f2f5fa; }
  .flags { font-weight: normal; font-size: 12px; color: #b03a2e; margin-left: 8px; }
  .frame { position: relative; height: 640px; background: #fafafa; border-bottom: 1px solid #eee; }
  .frame iframe { width: 100%; height: 100%; border: 0; }
  .placeholder { position: absolute; inset: 0; display: flex; align-items: center; justify-content: center;
                 cursor: pointer; color: #2e5496; font-size: 14px; flex-direction: column; gap: 6px; }
  .placeholder:hover { background: #eef4ff; }
  .judgments { display: flex; gap: 6px; padding: 8px 12px; }
  .jbtn { flex: 1; border: 1px solid #ccc; background: #fff; padding: 6px; border-radius: 5px;
          cursor: pointer; font-size: 12px; }
  .jbtn.good.on { background: var(--good); color: #fff; border-color: var(--good); }
  .jbtn.doubtful.on { background: var(--doubtful); color: #fff; border-color: var(--doubtful); }
  .jbtn.bad.on { background: var(--bad); color: #fff; border-color: var(--bad); }
  .note { width: 100%; border: 0; border-top: 1px solid #eee; padding: 6px 12px; font-size: 12px; }
  .card[data-jug="good"] { outline: 3px solid var(--good); }
  .card[data-jug="doubtful"] { outline: 3px solid var(--doubtful); }
  .card[data-jug="bad"] { outline: 3px solid var(--bad); }
  footer { padding: 20px; text-align: center; color: #888; font-size: 12px; }
</style></head><body>
<header>
  <h1>Root figures by plant &mdash; __NPLANTS__ plants, __NFIGS__ figures</h1>
  <div class="bar" id="tabs"></div>
  <div class="bar" style="margin-top:8px">
    <input id="search" placeholder="filter plant (e.g. S11)...">
    <button class="tool" id="layout" onclick="toggleLayout()">Layout: side by side</button>
    <button class="export" onclick="exportCsv()">Export CSV</button>
  </div>
</header>
<div id="content"></div>
<footer>Judgments are stored in this browser and shared with figure_report.html. Export to CSV to keep them.</footer>
<script>
const DATA = __DATA__;
const BLOCKS = Object.keys(DATA);
let active = BLOCKS[0];

function jkey(f) { return "jug:" + f.batch + ":" + f.name; }
function nkey(f) { return "note:" + f.batch + ":" + f.name; }
function esc(s) { return String(s).replace(/&/g,"&amp;").replace(/"/g,"&quot;").replace(/</g,"&lt;"); }

function counts(block) {
  let done = 0, total = 0;
  for (const p of DATA[block]) for (const f of p.figs) { total++; if (localStorage.getItem(jkey(f))) done++; }
  return [done, total];
}

function renderTabs() {
  const box = document.getElementById("tabs");
  box.innerHTML = "";
  for (const b of BLOCKS) {
    const [d, t] = counts(b);
    const btn = document.createElement("button");
    btn.className = "tab" + (b === active ? " active" : "");
    btn.innerHTML = b + ' <span class="cpt">(' + d + '/' + t + ' figures, ' + DATA[b].length + ' plants)</span>';
    btn.onclick = () => { active = b; render(); window.scrollTo(0, 0); };
    box.appendChild(btn);
  }
}

function figCard(f, pi, fi) {
  const jug = localStorage.getItem(jkey(f)) || "";
  const note = localStorage.getItem(nkey(f)) || "";
  const id = pi + "_" + fi;
  return '<div class="card" id="c' + id + '"' + (jug ? ' data-jug="' + jug + '"' : '') + '>' +
    '<h3>' + esc(f.time) + ' &middot; ' + esc(f.name) + (f.flags ? '<span class="flags">' + esc(f.flags) + '</span>' : '') + '</h3>' +
    '<div class="frame"><div class="placeholder" data-src="' + esc(f.src) + '" onclick="loadFig(this)">' +
      '<div>&#128065; click to display</div><div style="font-size:11px;opacity:.7">3D figure</div></div></div>' +
    '<div class="judgments">' +
      ['good','doubtful','bad'].map(v => '<button class="jbtn ' + v + (jug===v?' on':'') +
        '" onclick="judge(' + pi + ',' + fi + ',\'' + v + '\')">' + v[0].toUpperCase() + v.slice(1) + '</button>').join('') +
    '</div>' +
    '<input class="note" placeholder="note..." value="' + esc(note) + '" oninput="setNote(' + pi + ',' + fi + ',this.value)">' +
  '</div>';
}

function render() {
  renderTabs();
  const c = document.getElementById("content");
  const filter = document.getElementById("search").value.toLowerCase();
  let html = "";
  DATA[active].forEach((p, pi) => {
    if (filter && !p.plant.toLowerCase().includes(filter)) return;
    const n = p.figs.length;
    const cls = n === 1 ? "n1" : (n === 2 ? "n2" : "n3plus");
    html += '<section class="plant"><h2>' + esc(active) + ' ' + esc(p.plant) +
      '<span class="count">' + n + (n > 1 ? ' times' : ' time') + '</span>' +
      (n > 1 ? '<button onclick="loadAll(this)">display all</button>' : '') + '</h2>' +
      '<div class="figs ' + cls + '">' + p.figs.map((f, fi) => figCard(f, pi, fi)).join('') + '</div></section>';
  });
  c.innerHTML = html || '<p>No plant matches the filter.</p>';
}

function loadFig(el) {
  const f = document.createElement("iframe");
  f.src = el.dataset.src; f.loading = "lazy";
  el.parentNode.appendChild(f); el.remove();
}

function loadAll(btn) {
  btn.closest(".plant").querySelectorAll(".placeholder").forEach(loadFig);
  btn.remove();
}

function judge(pi, fi, val) {
  const f = DATA[active][pi].figs[fi];
  const card = document.getElementById("c" + pi + "_" + fi);
  const k = jkey(f);
  card.querySelectorAll(".jbtn").forEach(b => b.classList.remove("on"));
  if (localStorage.getItem(k) === val) {
    localStorage.removeItem(k); card.removeAttribute("data-jug");
  } else {
    localStorage.setItem(k, val); card.setAttribute("data-jug", val);
    card.querySelector(".jbtn." + val).classList.add("on");
  }
  renderTabs();
}

function setNote(pi, fi, val) {
  const k = nkey(DATA[active][pi].figs[fi]);
  if (val) localStorage.setItem(k, val); else localStorage.removeItem(k);
}

function toggleLayout() {
  document.body.classList.toggle("stacked");
  document.getElementById("layout").textContent =
    document.body.classList.contains("stacked") ? "Layout: stacked" : "Layout: side by side";
}

function exportCsv() {
  const rows = [["block","plant","time","batch","sample","judgment","note"]];
  for (const b of BLOCKS) for (const p of DATA[b]) for (const f of p.figs) {
    const j = localStorage.getItem(jkey(f)) || "";
    const n = localStorage.getItem(nkey(f)) || "";
    if (j || n) rows.push([b, p.plant, f.time, f.batch, f.name, j, '"' + n.replace(/"/g,'""') + '"']);
  }
  const blob = new Blob([rows.map(r => r.join(",")).join("\n")], {type: "text/csv"});
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob); a.download = "plant_judgments.csv"; a.click();
}

document.getElementById("search").addEventListener("input", render);
render();
</script></body></html>"""


def build(batches):
    from tools.qc_rank import load_ranking
    data = group_by_plant(batches, load_ranking(RESULTS_ROOT))
    nplants = sum(len(v) for v in data.values())
    nfigs = sum(len(p["figs"]) for v in data.values() for p in v)
    return (PAGE.replace("__DATA__", json.dumps(data, ensure_ascii=False))
                .replace("__NPLANTS__", str(nplants)).replace("__NFIGS__", str(nfigs))), data


def generate(results_root=None):
    """Write the plant report; returns its path, or None if no figure was found."""
    global RESULTS_ROOT, OUTPUT_HTML
    if results_root is not None:
        RESULTS_ROOT = results_root
        OUTPUT_HTML = os.path.join(results_root, "plant_report.html")
    fr.RESULTS_ROOT = RESULTS_ROOT
    batches = fr.scan_figures()
    if not batches:
        print(f"No figure found in {RESULTS_ROOT}/<batch>/figures/.")
        print("Check that the analysis ran with figures enabled (SAVE_FIGURES=1).")
        return None
    html, data = build(batches)
    with open(OUTPUT_HTML, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"Plant report written: {OUTPUT_HTML}")
    for block, plants in data.items():
        multi = sum(1 for p in plants if len(p["figs"]) > 1)
        print(f"  {block}: {len(plants)} plants ({multi} with several times)")
    return OUTPUT_HTML


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--results', default=None, help='results folder (default: results/ of the project)')
    if generate(ap.parse_args().results) is not None:
        print("\nOpen plant_report.html in a browser.")


if __name__ == "__main__":
    main()
