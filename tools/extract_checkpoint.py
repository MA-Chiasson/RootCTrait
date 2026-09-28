"""extract_checkpoint.py

Turn a run's checkpoint into a trait table immediately, without finishing the run.

The pipeline appends one JSON line per processed sample to a checkpoint file
(checkpoint_traits.jsonl) as it goes. If you stop a run, this tool reads whatever
has been written so far and produces a trait table of the samples already done.

Each checkpoint line holds: name, n_brut, n_ret, and T (the ~40 traits). The output
mirrors the normal table: ID, n_brut, n_retire, %retire, then the traits.

Run:  python extract_checkpoint.py                     (scan results/<batch>/)
      python extract_checkpoint.py path/to/checkpoint_traits.jsonl
      python extract_checkpoint.py --format xlsx        (or: csv | both)
Output: a *_partial table next to each checkpoint.
"""
import os
import sys
import glob
import json
import re
import pandas as pd

import os as _os, sys as _sys
# tools/ lives under the project root; resolve the root and make it importable
PROJECT_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
if PROJECT_ROOT not in _sys.path:
    _sys.path.insert(0, PROJECT_ROOT)


def load_factors():
    """Read the (trait -> scaling factor) map from run_pipeline.py so the partial
    table uses the SAME units as the normal output. Falls back to no scaling."""
    factors, order = {}, []
    src_path = os.path.join(PROJECT_ROOT, "run_pipeline.py")
    try:
        src = open(src_path, encoding="utf-8", errors="ignore").read()
        block = re.search(r"COLS\s*=\s*\[(.*?)\]", src, re.S).group(1)
        for name, fac in re.findall(r"\(\s*'([^']+)'\s*,\s*'[^']*'\s*,\s*([0-9.]+)\)", block):
            factors[name] = float(fac); order.append(name)
    except Exception as e:
        print(f"    (could not read scaling factors, values left raw: {e})")
    return factors, order


FACTORS, COL_ORDER = load_factors()


def results_root():
    pfile = os.environ.get("PARAMS", os.path.join(PROJECT_ROOT, "params.txt"))
    root = os.path.join(PROJECT_ROOT, "results")
    if os.path.exists(pfile):
        for line in open(pfile, encoding="utf-8", errors="ignore"):
            line = line.strip()
            if line.upper().startswith("RESULTS_ROOT"):
                root = line.split("=", 1)[1].strip()
    if not os.path.isabs(root):
        root = os.path.join(PROJECT_ROOT, root)
    return root


def read_checkpoint(path):
    """Read a JSONL checkpoint into a DataFrame. Skips any malformed line
    (e.g. a final line half-written when the run was interrupted)."""
    rows, bad = [], 0
    with open(path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                d = json.loads(line)
            except json.JSONDecodeError:
                bad += 1
                continue
            nb = d.get("n_brut"); nr = d.get("n_ret")
            row = {"ID": d.get("name"),
                   "n_brut": nb,
                   "n_retire": nr,
                   "%retire": (round(100 * nr / nb, 1) if nb else 0)}
            T = d.get("T", {})
            for k, v in T.items():
                row[k] = (v * FACTORS[k]) if (k in FACTORS and isinstance(v, (int, float))) else v
            rows.append(row)
    if bad:
        print(f"    ({bad} incomplete line(s) skipped)")
    if not rows:
        return None
    df = pd.DataFrame(rows)
    lead = [c for c in ["ID", "n_brut", "n_retire", "%retire"] if c in df.columns]
    ordered = [c for c in COL_ORDER if c in df.columns]              # same order as normal output
    rest = [c for c in df.columns if c not in lead + ordered]
    return df[lead + ordered + rest]


def write(df, base, fmt):
    if fmt in ("xlsx", "both"):
        df.to_excel(base + ".xlsx", index=False); print(f"  -> {base}.xlsx")
    if fmt in ("csv", "both"):
        df.to_csv(base + ".csv", index=False); print(f"  -> {base}.csv")


def main():
    args = [a for a in sys.argv[1:]]
    fmt = "both"
    if "--format" in args:
        i = args.index("--format"); fmt = args[i + 1].lower(); del args[i:i + 2]

    # explicit file, or scan the results folder for checkpoints
    if args:
        paths = [args[0]]
    else:
        root = results_root()
        paths = sorted(glob.glob(os.path.join(root, "*", "*checkpoint*.jsonl")))
        paths += sorted(glob.glob(os.path.join(root, "*checkpoint*.jsonl")))
    if not paths:
        print("No checkpoint file found.")
        return

    for p in paths:
        if not os.path.exists(p):
            print(f"Not found: {p}"); continue
        df = read_checkpoint(p)
        if df is None:
            print(f"  {p}: empty (no sample done yet)"); continue
        batch = os.path.basename(os.path.dirname(p)) or "checkpoint"
        print(f"{batch}: {len(df)} samples done so far")
        base = os.path.join(os.path.dirname(p), f"traits_{batch}_partial")
        write(df, base, fmt)


if __name__ == "__main__":
    main()
