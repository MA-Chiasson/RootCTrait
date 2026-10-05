"""sensitivity.py

One-at-a-time sensitivity analysis of the RootCTrait thresholds.

Each threshold is moved below and above its default while all the others keep their
default value, the full pipeline is rerun on a random subset of samples, and every
trait is compared with its default value. For each threshold, variant and trait the
summary reports:

  rho        Spearman rank correlation across samples between the variant and the
             default (1 = the ranking of the samples is unchanged, which is what
             matters for association analysis);
  med_rel    median absolute relative change of the trait value, in percent.

Thresholds read from params.txt (PRUNE_VOX, MIN_SEG_LEN_MM, BC_MIN, LIN_MAX,
LEN_MAX, DENS_MAX, DENS_LEN_MAX, RESCUE_MIN_MM) and thresholds fixed in the code (collar layer, hypocotyl angle, collar
climb fraction, high branch margin) are all covered.

Batches whose files hold several versions of the mask use the variable recorded in
RESULTS_ROOT/<batch>/selection_versions.csv (written by retraiter_B2T3.py), so that
the same mask is analysed as in the trait table.

Run from the project root (not from tools/), for example:
    python -m tools.sensitivity --batch block1_t1 --n 30 --workers 4
    python -m tools.sensitivity --batch block1_t1 --batch block2_t2 --n 20
Runs are saved as they end (sensitivity_partial.jsonl): rerun the same command to
resume an interrupted analysis; delete that file to start over.
Output (in RESULTS_ROOT/sensitivity/):
    sensitivity_traits.csv   one row per sample and variant, all traits
    sensitivity_summary.csv  one row per threshold, variant and trait
"""
import argparse
import functools
import multiprocessing as mp
import os
import random
import sys
import time
import json
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

import run_pipeline as rp                                    # noqa: E402
from rootctrait import detection_hypocotyle as dh            # noqa: E402
from rootctrait import root_traits_full as rtf               # noqa: E402
from rootctrait.io_volume import load_volume                 # noqa: E402

# threshold -> (default, low, high)
FILE_PARAMS = {
    'PRUNE_VOX': (rp.PRUNE_VOX, 3, 7),
    'MIN_SEG_LEN_MM': (rp.MIN_SEG_LEN_MM, 1.5, 2.5),
    'BC_MIN': (rp.BC_MIN, 2, 4),
    'LIN_MAX': (rp.LIN_MAX, 0.6, 0.8),
    'LEN_MAX': (rp.LEN_MAX, 10.0, 20.0),
    'DENS_MAX': (rp.DENS_MAX, 30.0, 40.0),
    'DENS_LEN_MAX': (rp.DENS_LEN_MAX, 4.0, 8.0),
    'RESCUE_MIN_MM': (rp.RESCUE_MIN_MM, 3.0, 8.0),
}
CODE_PARAMS = {
    'collar_layer': (0.20, 0.15, 0.25),
    'hypocotyl_angle_deg': (45.0, 35.0, 55.0),
    'climb_fraction': (0.7, 0.6, 0.8),
    'high_branch_mm': (3.0, 2.0, 4.0),
}
CTX_KEY = {'PRUNE_VOX': 'prune_vox', 'MIN_SEG_LEN_MM': 'min_seg_len_mm',
           'BC_MIN': 'bc_min', 'LIN_MAX': 'lin_max', 'LEN_MAX': 'len_max',
           'DENS_MAX': 'dens_max', 'DENS_LEN_MAX': 'dens_len_max', 'RESCUE_MIN_MM': 'rescue_min_mm'}


def variants():
    """('default', None, None) followed by every one-at-a-time low and high variant."""
    out = [('default', None, None)]
    for name, (_, lo, hi) in {**FILE_PARAMS, **CODE_PARAMS}.items():
        out += [(name, 'low', lo), (name, 'high', hi)]
    return out


def _base_ctx(batch):
    return {'data_dir': os.path.join(rp.DATA_ROOT, batch['name']),
            'pattern': batch['pattern'], 'axis_order': batch['axis_order'],
            'var_name': batch.get('var_name'),
            'voxel_size': tuple(rp.VOXEL_SIZE), 'prune_vox': rp.PRUNE_VOX,
            'min_seg_len_mm': rp.MIN_SEG_LEN_MM, 'bc_min': rp.BC_MIN,
            'lin_max': rp.LIN_MAX, 'len_max': rp.LEN_MAX,
            'drop_orphans': rp.DROP_ORPHANS, 'dens_max': rp.DENS_MAX,
            'dens_len_max': rp.DENS_LEN_MAX, 'rescue_min_mm': rp.RESCUE_MIN_MM,
            'save_figures': False, 'fig_dir': None}


def _sample_ctx(batch, sample):
    ctx = _base_ctx(batch)
    # mask variable chosen per sample by a dedicated script (selection_versions.csv)
    sel = batch.get('selection', {})
    stem = os.path.splitext(batch['pattern'].format(name=sample))[0]   # file name without .mat
    var = sel.get(sample, sel.get(stem))
    if var:
        ctx['var_name'] = var
    return ctx


def _variant(batch, sample, variant, BW, pre, ctx0):
    """One variant on a mask already loaded and prepared. The code-level thresholds
    are injected by wrapping the pipeline functions inside this worker only."""
    param, side, value = variant
    ctx = dict(ctx0)
    code = {k: v[0] for k, v in CODE_PARAMS.items()}
    if param in CTX_KEY:
        ctx[CTX_KEY[param]] = value
    elif param in code:
        code[param] = value
    rp.detect_base = functools.partial(_ORIG['detect_base'], layer=code['collar_layer'])
    rp.collar_and_hypocotyl = functools.partial(
        _ORIG['collar_and_hypocotyl'], angle_max=code['hypocotyl_angle_deg'],
        frac=code['climb_fraction'], haut_min_mm=code['high_branch_mm'])
    base = {'batch': batch['name'], 'sample': sample, 'parameter': param, 'side': side, 'value': value}
    try:
        R = rp.analyse_mask(BW, ctx, pre=pre)
    except Exception as e:                                   # keep going, log the failure
        return dict(base, error=repr(e))
    row = dict(base, n_raw=R['n_raw'], n_removed=R['n_removed'])
    row.update({k: R['T'].get(k) for k, _, _ in rp.COLS})
    return row


def _run_sample(stask):
    """All the variants of one sample. The mask is loaded once, and the steps that do
    not depend on the thresholds (distance map, skeleton, pruning per length) are
    computed once and shared by the variants: same results as separate runs, several
    times faster."""
    batch, sample, vlist = stask
    ctx0 = _sample_ctx(batch, sample)
    path = os.path.join(ctx0['data_dir'], ctx0['pattern'].format(name=sample))
    V = load_volume(path, axis_order=ctx0['axis_order'], var_name=ctx0.get('var_name'))
    BW = V > 0.5 * float(np.max(V))
    del V
    pre = rp.prepare_mask(BW, ctx0['voxel_size'])
    return [_variant(batch, sample, v, BW, pre, ctx0) for v in vlist]


_ORIG = {'detect_base': rp.detect_base,
         'collar_and_hypocotyl': dh.collar_and_hypocotyl,
         'compute_all_traits': rtf.compute_all_traits}


def _child(stask, q):
    q.put(_run_sample(stask))


def _run_isolated(stask, timeout):
    """Run all the variants of one sample in their own process, with a time limit."""
    batch, sample, vlist = stask
    q = mp.Queue()
    p = mp.Process(target=_child, args=(stask, q), daemon=True)
    p.start()
    try:
        return q.get(timeout=timeout)
    except Exception:
        reason = 'timeout' if p.is_alive() else f'process ended without a result (exit code {p.exitcode}, out of memory?)'
        return [{'batch': batch['name'], 'sample': sample, 'parameter': v[0], 'side': v[1],
                 'value': v[2], 'error': reason} for v in vlist]
    finally:
        if p.is_alive():
            p.terminate()
        p.join(5)


def summarize(df):
    from scipy.stats import spearmanr
    traits = [k for k, _, _ in rp.COLS]
    key = ['batch', 'sample']
    ref = df[df['parameter'] == 'default'].set_index(key)[traits]
    rows = []
    for (param, side), g in df[df['parameter'] != 'default'].groupby(['parameter', 'side']):
        g = g.set_index(key)[traits]
        common = ref.index.intersection(g.index)
        for t in traits:
            a = pd.to_numeric(ref.loc[common, t], errors='coerce')
            b = pd.to_numeric(g.loc[common, t], errors='coerce')
            ok = a.notna() & b.notna()
            rho = spearmanr(a[ok], b[ok])[0] if ok.sum() >= 3 and a[ok].nunique() > 1 else np.nan
            rel = (100 * (b[ok] - a[ok]).abs() / a[ok].abs().replace(0, np.nan)).median()
            rows.append({'parameter': param, 'side': side, 'trait': t, 'n': int(ok.sum()),
                         'rho': rho, 'med_rel_pct': rel})
    return pd.DataFrame(rows)


def load_selection(batch_name):
    """{sample: variable} from RESULTS_ROOT/<batch>/selection_versions.csv, written by a
    dedicated re-processing script (e.g. retraiter_B2T3.py) for batches whose files hold
    several versions of the mask; empty when the file does not exist."""
    import csv
    path = os.path.join(rp.RESULTS_ROOT, batch_name, 'selection_versions.csv')
    if not os.path.exists(path):
        return {}
    with open(path, newline='') as f:
        return {r['sample']: r['variable_used'] for r in csv.DictReader(f)
                if r['variable_used'] and r['variable_used'] != '(largest)'}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[2])
    ap.add_argument('--batch', action='append', required=True,
                    help='batch name as defined in params.txt (repeatable)')
    ap.add_argument('--n', type=int, default=30, help='samples drawn per batch (default 30)')
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--workers', type=int, default=1)
    args = ap.parse_args()

    defs = {b['name']: b for b in rp.BATCH_DEFS}
    rng = random.Random(args.seed)
    tasks = []
    for name in args.batch:
        if name not in defs:
            sys.exit(f"Batch '{name}' is not defined in params.txt")
        b = dict(defs[name], selection=load_selection(name))
        if b['selection']:
            print(f"{name}: mask variables read from selection_versions.csv")
        samples = rp.list_samples(os.path.join(rp.DATA_ROOT, name), b['pattern'])
        pick = sorted(rng.sample(samples, min(args.n, len(samples))))
        print(f"{name}: {len(pick)} of {len(samples)} samples")
        tasks += [(b, s, v) for s in pick for v in variants()]
    print(f"{len(tasks)} runs ({len(variants())} variants per sample)")

    # Every run is written to the partial file as soon as it ends, so an interrupted
    # analysis resumes where it stopped (same --seed = same samples). Each run is a
    # separate process with the per sample time limit of params.txt (TIMEOUT): a run
    # that hangs or crashes (for example out of memory) is recorded as failed and the
    # analysis goes on, and its memory is freed at once.
    out = os.path.join(rp.RESULTS_ROOT, 'sensitivity')
    os.makedirs(out, exist_ok=True)
    part = os.path.join(out, 'sensitivity_partial.jsonl')
    key = lambda b, smp, v: f"{b}|{smp}|{v[0]}|{v[1]}"
    done = {}
    if os.path.exists(part):
        with open(part, encoding='utf-8') as f:
            for line in f:
                try:
                    r = json.loads(line)
                except ValueError:
                    continue
                if not r.get('error'):
                    done[key(r['batch'], r['sample'], (r['parameter'], r['side']))] = r
    todo = [t for t in tasks if key(t[0]['name'], t[1], t[2][:2]) not in done]
    rows = list(done.values())
    if done:
        print(f"{len(done)} runs already done (resumed from {part}), {len(todo)} to do")
    # group the remaining variants by sample: one process per sample
    groups = {}
    for b, smp, v in todo:
        groups.setdefault((b['name'], smp), (b, smp, []))[2].append(v)
    stasks = list(groups.values())
    nvar = len(variants())
    limit = max(rp.TIMEOUT, 120 * nvar)
    t0 = time.time()
    print(f"{len(stasks)} samples to process, all their variants at once "
          f"(time limit {limit} s per sample); progress is shown after every sample", flush=True)
    with open(part, 'a', encoding='utf-8') as fp, ThreadPoolExecutor(max(1, args.workers)) as ex:
        futs = [ex.submit(_run_isolated, st, limit) for st in stasks]
        for k, fu in enumerate(as_completed(futs), 1):
            res = fu.result()
            for r in res:
                rows.append(r)
                fp.write(json.dumps(r, default=float) + '\n')
            fp.flush()
            el = time.time() - t0
            nerr = sum(1 for r in res if r.get('error'))
            err = f"  {nerr} FAILED: {next(r['error'] for r in res if r.get('error'))[:60]}" if nerr else ''
            print(f"  {k}/{len(stasks)} samples  {res[0]['batch']} {res[0]['sample']} ({len(res)} variants){err}"
                  f"  | elapsed {el / 60:.1f} min, about {el / k * (len(stasks) - k) / 60:.0f} min left", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(out, 'sensitivity_traits.csv'), index=False)
    if 'error' in df.columns:
        bad = df[df['error'].notna()]
        print(f"{len(bad)} failed runs (see the error column)")
        df = df[df['error'].isna()]
    summ = summarize(df)
    summ.to_csv(os.path.join(out, 'sensitivity_summary.csv'), index=False)
    worst = summ.groupby(['parameter', 'side'])['rho'].min().sort_values()
    print('\nLowest rank correlation per variant (over all traits):')
    print(worst.round(3).to_string())
    print(f"\n-> {out}")


if __name__ == '__main__':
    main()
