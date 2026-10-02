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
LEN_MAX) and thresholds fixed in the code (collar layer, hypocotyl angle, collar
climb fraction, high branch margin) are all covered.

Run from the project root, for example:
    python -m tools.sensitivity --batch block1_t1 --n 30 --workers 4
    python -m tools.sensitivity --batch block1_t1 --batch block2_t2 --n 20
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

import numpy as np
import pandas as pd

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

import run_pipeline as rp                                    # noqa: E402
from rootctrait import detection_hypocotyle as dh            # noqa: E402
from rootctrait import root_traits_full as rtf               # noqa: E402

# threshold -> (default, low, high)
FILE_PARAMS = {
    'PRUNE_VOX': (rp.PRUNE_VOX, 3, 7),
    'MIN_SEG_LEN_MM': (rp.MIN_SEG_LEN_MM, 1.5, 2.5),
    'BC_MIN': (rp.BC_MIN, 2, 4),
    'LIN_MAX': (rp.LIN_MAX, 0.6, 0.8),
    'LEN_MAX': (rp.LEN_MAX, 10.0, 20.0),
}
CODE_PARAMS = {
    'collar_layer': (0.20, 0.15, 0.25),
    'hypocotyl_angle_deg': (45.0, 35.0, 55.0),
    'climb_fraction': (0.7, 0.6, 0.8),
    'high_branch_mm': (3.0, 2.0, 4.0),
}
CTX_KEY = {'PRUNE_VOX': 'prune_vox', 'MIN_SEG_LEN_MM': 'min_seg_len_mm',
           'BC_MIN': 'bc_min', 'LIN_MAX': 'lin_max', 'LEN_MAX': 'len_max'}


def variants():
    """('default', None, None) followed by every one-at-a-time low and high variant."""
    out = [('default', None, None)]
    for name, (_, lo, hi) in {**FILE_PARAMS, **CODE_PARAMS}.items():
        out += [(name, 'low', lo), (name, 'high', hi)]
    return out


def _base_ctx(batch):
    return {'data_dir': os.path.join(rp.DATA_ROOT, batch['name']),
            'pattern': batch['pattern'], 'axis_order': batch['axis_order'],
            'voxel_size': tuple(rp.VOXEL_SIZE), 'prune_vox': rp.PRUNE_VOX,
            'min_seg_len_mm': rp.MIN_SEG_LEN_MM, 'bc_min': rp.BC_MIN,
            'lin_max': rp.LIN_MAX, 'len_max': rp.LEN_MAX,
            'drop_orphans': rp.DROP_ORPHANS, 'save_figures': False, 'fig_dir': None}


def _run_one(task):
    """Process one sample with one variant. The code-level thresholds are injected by
    wrapping the pipeline functions inside this worker only."""
    batch, sample, (param, side, value) = task
    ctx = _base_ctx(batch)
    code = {k: v[0] for k, v in CODE_PARAMS.items()}
    if param in CTX_KEY:
        ctx[CTX_KEY[param]] = value
    elif param in code:
        code[param] = value
    rp.detect_base = functools.partial(_ORIG['detect_base'], layer=code['collar_layer'])
    rp.collar_and_hypocotyl = functools.partial(
        _ORIG['collar_and_hypocotyl'], angle_max=code['hypocotyl_angle_deg'],
        frac=code['climb_fraction'], haut_min_mm=code['high_branch_mm'])
    try:
        n_raw, n_rem, T = rp.process(sample, ctx)
    except Exception as e:                                   # keep going, log the failure
        return {'batch': batch['name'], 'sample': sample, 'parameter': param,
                'side': side, 'value': value, 'error': repr(e)}
    row = {'batch': batch['name'], 'sample': sample, 'parameter': param, 'side': side,
           'value': value, 'n_raw': n_raw, 'n_removed': n_rem}
    row.update({k: T.get(k) for k, _, _ in rp.COLS})
    return row


_ORIG = {'detect_base': rp.detect_base,
         'collar_and_hypocotyl': dh.collar_and_hypocotyl,
         'compute_all_traits': rtf.compute_all_traits}


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
        b = defs[name]
        samples = rp.list_samples(os.path.join(rp.DATA_ROOT, name), b['pattern'])
        pick = sorted(rng.sample(samples, min(args.n, len(samples))))
        print(f"{name}: {len(pick)} of {len(samples)} samples")
        tasks += [(b, s, v) for s in pick for v in variants()]
    print(f"{len(tasks)} runs ({len(variants())} variants per sample)")

    rows = []
    if args.workers > 1:
        with mp.Pool(args.workers) as pool:
            for k, r in enumerate(pool.imap_unordered(_run_one, tasks), 1):
                rows.append(r)
                if k % 20 == 0:
                    print(f"  {k}/{len(tasks)}", flush=True)
    else:
        for k, t in enumerate(tasks, 1):
            rows.append(_run_one(t))
            if k % 20 == 0:
                print(f"  {k}/{len(tasks)}", flush=True)

    out = os.path.join(rp.RESULTS_ROOT, 'sensitivity')
    os.makedirs(out, exist_ok=True)
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
