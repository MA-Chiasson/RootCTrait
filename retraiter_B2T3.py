"""retraiter_B2T3.py

Dedicated re-processing for the B2T3 batch, whose .mat files each hold several
versions of the root mask (a noisy main mask, a denoised SD version, and a cleaned
ROOT... version). The generic loader picks the largest 3D array, which cannot tell
these same-shape versions apart and may load the polluted one. This script selects
the CLEANED version explicitly, then runs the CURRENT pipeline (run_pipeline.process,
RootCTrait 2.4.0 or later) on it, with the same settings, outputs and safeguards as a
normal batch: shared parameters from params.txt, per-sample time limit, checkpoint
(an interrupted run resumes where it stopped), failure log, quality control columns,
RSML export and review figures.

Selection order (choisir_version): cleaned ROOT... variant (longest name = most
cleaning) -> SD -> main pattern -> densest. A traceability CSV records the choice.

Layout: B2T3 files in <DATA_ROOT>/B2T3/ (e.g. B2T3S12.mat). Settings come from
params.txt (or the file given by the PARAMS environment variable).

Run:  python retraiter_B2T3.py
Output, in <RESULTS_ROOT>/B2T3/:
        traits_B2T3.xlsx           trait table (same columns as a normal batch)
        selection_versions.csv     mask variable used for each sample, and why
        checkpoint_traits.jsonl    resume state (delete it to recompute everything)
        failures.jsonl             samples that timed out or failed
        params_used.json           settings used
        rsml/                      one RSML file per sample (if EXPORT_RSML=1)
        figures/                   review figures (if SAVE_FIGURES=1)
Then, as for any batch: python -m tools.qc_rank  and  python -m tools.figure_report
"""
import os
import re
import glob
import csv
import json
import time
import numpy as np
import h5py

import run_pipeline as rp   # current algorithm and batch machinery
from rootctrait.sorting import natural_key

BATCH = "B2T3"
BASE_PAT = re.compile(r"^ROOT(_SEL)?_B2T3G\d+S\d+$")


def _is_binary(ds, step=8):
    """Quick test on a strided subsample (1 voxel in step**3): True when it holds only
    0 and 1 with some foreground. Avoids reading whole volumes, which is slow for
    several 3D versions per file (and on synced folders such as OneDrive)."""
    if ds.dtype.kind not in 'biuf':
        return False
    sub = np.asarray(ds[::step, ::step, ::step])
    u = np.unique(sub)
    return bool(set(u.tolist()) <= {0, 1} and (sub > 0).any())


def _count(ds, slab=64):
    """Number of foreground voxels, read slab by slab (only when needed)."""
    return int(sum(np.count_nonzero(ds[i:i + slab]) for i in range(0, ds.shape[0], slab)))


def choisir_version(fpath):
    """Return (variable_name, reason) for the best mask in the file."""
    try:
        with h5py.File(fpath, 'r') as f:
            keys = [k for k in f.keys() if not k.startswith('#')]
            bins = [k for k in keys
                    if isinstance(f[k], h5py.Dataset) and len(f[k].shape) == 3 and _is_binary(f[k])]
            if not bins:
                raise ValueError("no usable binary mask")
            epures = [k for k in bins if k.upper().startswith("ROOT") and not BASE_PAT.match(k)]
            if epures:
                longest = max(len(k) for k in epures)
                top = [k for k in epures if len(k) == longest]
                k = top[0] if len(top) == 1 else max(top, key=lambda k_: _count(f[k_]))
                return k, "cleaned"
            if 'SD' in bins:
                return 'SD', "SD (no cleaned version)"
            princ = [k for k in bins if BASE_PAT.match(k)]
            if princ:
                return princ[0], "main (no cleaned nor SD)"
            return max(bins, key=lambda k_: _count(f[k_])), "fallback (densest)"
    except OSError:
        # older .mat format (not HDF5): let the generic loader pick the array
        return None, "generic loader (non-HDF5 file)"


def main():
    data_dir = os.path.join(rp.DATA_ROOT, BATCH)
    res_dir = os.path.join(rp.RESULTS_ROOT, BATCH)
    os.makedirs(res_dir, exist_ok=True)
    checkpoint = os.path.join(res_dir, "checkpoint_traits.jsonl")
    failures_path = os.path.join(res_dir, "failures.jsonl")
    files = sorted(glob.glob(os.path.join(data_dir, "*.mat")), key=lambda p: natural_key(os.path.basename(p)))
    if not files:
        print(f"No .mat file in {data_dir}")
        return

    # same settings as every other batch (params.txt), recorded for reproducibility
    eff = rp._effective_params()
    base_ctx = {"data_dir": data_dir, "pattern": "{name}.mat", "axis_order": None,
                "fig_dir": os.path.join(res_dir, "figures"),
                "rsml_dir": os.path.join(res_dir, "rsml"), **eff}
    used = dict(eff, voxel_size=list(rp.VOXEL_SIZE), batch=BATCH,
                data_dir=os.path.abspath(data_dir), pattern="{name}.mat",
                mask_selection="selection_versions.csv", timeout=rp.TIMEOUT,
                generated_at=time.strftime('%Y-%m-%d %H:%M:%S'))
    with open(os.path.join(res_dir, "params_used.json"), "w", encoding="utf-8") as fh:
        json.dump(used, fh, indent=2, ensure_ascii=False)

    store = rp.load_store(checkpoint)
    names = [os.path.splitext(os.path.basename(fp))[0] for fp in files]
    todo = [n for n in names if n not in store]
    print(f"=== {BATCH} re-processing: {len(files)} files | {len(names) - len(todo)} already done | "
          f"{len(todo)} to do ===")

    # mask selection, file by file; choices already recorded in the CSV are reused
    sel_csv = os.path.join(res_dir, "selection_versions.csv")
    selection = {}
    if os.path.exists(sel_csv):
        with open(sel_csv, newline="") as fh:
            for row in csv.DictReader(fh):
                v = row["variable_used"]
                selection[row["sample"]] = (None if v == "(largest)" else v, row["reason"])

    def save_selection():
        with open(sel_csv, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["sample", "variable_used", "reason"])
            for n in names:
                if n in selection:
                    v, r = selection[n]
                    w.writerow([n, v if v else "(largest)", r])

    for i, name in enumerate(todo, 1):
        if name not in selection:
            print(f"  [{i}/{len(todo)}] {name:12s} choosing the mask version...", flush=True)
            selection[name] = choisir_version(os.path.join(data_dir, name + ".mat"))
            save_selection()
        var, reason = selection[name]
        print(f"  [{i}/{len(todo)}] {name:12s} mask: {var if var else '(largest)'} [{reason}]", flush=True)
        ctx = dict(base_ctx, var_name=var)
        t0 = time.time()
        status, payload = rp.run_with_timeout(name, ctx, rp.TIMEOUT)
        rp._handle_result(name, status, payload, t0, store, checkpoint, failures_path)

    results = [(n, store[n]['n_raw'], store[n]['n_removed'], store[n]['T'], store[n].get('seconds'))
               for n in names if n in store]
    out = os.path.join(res_dir, f"traits_{BATCH}.xlsx")
    rp.write_excel(results, out, f"Root traits - {BATCH} - cleaned mask version selected")
    print(f"  -> {out}  ({len(results)}/{len(files)} samples)")
    print(f"  -> {sel_csv}")
    if len(results) < len(files):
        print(f"  -> {len(files) - len(results)} not completed (see failures.jsonl; run again to resume them)")


if __name__ == "__main__":
    main()
