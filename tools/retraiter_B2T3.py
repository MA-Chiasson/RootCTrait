"""merge_batches.py

Merge all per-batch trait tables into a single master table.

Generic: it scans the results folder, takes every trait table it finds, and stacks
them into one table. The batch name is taken verbatim from each table's subfolder
name — no naming convention is assumed (temporal or not). A 'batch' column records
where each row came from; the per-sample ID already present in each table is kept.

Columns are aligned on the UNION of all batches (prudent mode): a column missing
from a batch is filled with blanks rather than dropped, and any column not present
in every batch is reported.

Run:  python merge_batches.py               (writes xlsx + csv)
      python merge_batches.py --format xlsx (or: csv | both)
Output: <results>/merged_traits.xlsx and/or .csv
"""
import os
import sys
import glob
import pandas as pd

import os as _os, sys as _sys
# tools/ lives under the project root; resolve the root and make it importable
PROJECT_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
if PROJECT_ROOT not in _sys.path:
    _sys.path.insert(0, PROJECT_ROOT)


def results_root():
    """Read RESULTS_ROOT from params.txt if present, else default to 'results'."""
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


def read_batch_table(path):
    """Read one trait table, skipping the title row and the units row.
    Returns a DataFrame with the batch's samples, or None if unreadable."""
    try:
        df = pd.read_excel(path, header=1)          # row 2 = headers (row 1 = title)
    except Exception as e:
        print(f"  ! could not read {path}: {e}")
        return None
    if "ID" not in df.columns:
        print(f"  ! {path}: no 'ID' column, skipped")
        return None
    df = df[df["ID"].notna()].copy()                # drop the units row (ID is blank)
    df = df.dropna(how="all")
    return df


def main():
    fmt = "both"
    if "--format" in sys.argv:
        fmt = sys.argv[sys.argv.index("--format") + 1].lower()
    root = results_root()
    if not os.path.isdir(root):
        print(f"Results folder not found: {root}")
        return

    # every xlsx under results/<batch>/ is a trait table
    paths = sorted(glob.glob(os.path.join(root, "*", "*.xlsx")))
    paths = [p for p in paths if "merged_traits" not in os.path.basename(p)]
    if not paths:
        print(f"No trait table found under {root}/<batch>/")
        return

    frames = []
    col_presence = {}   # column -> number of batches containing it
    batch_names = []
    for p in paths:
        batch = os.path.basename(os.path.dirname(p))   # subfolder name, verbatim
        df = read_batch_table(p)
        if df is None or len(df) == 0:
            continue
        df.insert(0, "batch", batch)                   # batch name as-is, first column
        frames.append(df)
        batch_names.append(batch)
        for c in df.columns:
            col_presence[c] = col_presence.get(c, 0) + 1
        print(f"  + {batch}: {len(df)} samples, {df.shape[1]-1} columns")

    if not frames:
        print("Nothing to merge.")
        return

    # prudent union: concat aligns on the union of columns, missing -> NaN
    merged = pd.concat(frames, ignore_index=True, sort=False)
    # keep a stable, readable column order: batch, ID, then the rest as first seen
    lead = [c for c in ["batch", "ID", "n_brut", "n_retire", "%retire"] if c in merged.columns]
    rest = [c for c in merged.columns if c not in lead]
    merged = merged[lead + rest]

    n_batches = len(frames)
    partial = [c for c, n in col_presence.items() if n < n_batches and c != "batch"]
    print(f"\nMerged {n_batches} batches -> {len(merged)} samples, {merged.shape[1]} columns")
    if partial:
        print("Columns NOT present in every batch (filled with blanks):")
        for c in sorted(partial):
            print(f"    {c}  (in {col_presence[c]}/{n_batches} batches)")
    else:
        print("All batches share the same columns.")

    base = os.path.join(root, "merged_traits")
    if fmt in ("xlsx", "both"):
        merged.to_excel(base + ".xlsx", index=False)
        print(f"  -> {base}.xlsx")
    if fmt in ("csv", "both"):
        merged.to_csv(base + ".csv", index=False)
        print(f"  -> {base}.csv")


if __name__ == "__main__":
    main()
