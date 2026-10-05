"""qc_rank.py
Ranks the processed samples for visual review, so that the review starts with the
samples most likely to hold a processing error.

Each sample receives flags from explicit rules applied to the quality control values
written by the pipeline (checkpoint_traits.jsonl of each batch):

  pivot return      upward return of the primary root tip above 3 mm (Section 2.11 of
                    the article; same threshold as the review figure)
  removal high      percentage of segments removed by decontamination
  rescued high      number of detached roots rescued
  collar raised     how far the collar was raised
  hypocotyl long    length of the excluded hypocotyl
  path / depth      ratio of primary root length to rooting depth (a primary path
                    that wanders along laterals or through a contact is long for its depth)
  trait outlier     LRP, PM, TRL, NRL, ANGO2 or DRP far from the rest of the batch

Except for the pivot return, a value is flagged when it lies above the batch median
by more than K robust standard deviations (1.4826 x median absolute deviation),
K = 3 by default; trait outliers are flagged on both sides. Batches are compared
separately, because scan quality differs between scanning sessions.

Samples are sorted by number of flags, then by their largest robust z score.
Nothing is changed in the traits: the ranking only orders the review.

Run from the project root:
    python -m tools.qc_rank                 # all batches in results/
    python -m tools.qc_rank --k 3.5         # stricter threshold
Writes results/qc_ranking.csv. If it exists, tools/figure_report.py orders the
figures of each batch by this ranking and shows the flags on each card.
"""
import os, sys, json, glob, argparse, csv
import numpy as np

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)
from rootctrait.legacy import upgrade_record

PIVOT_FLAG_MM = 3.0
HIGH_ONLY = [('removal high', '%removed'), ('rescued high', 'N_RESCUED'),
             ('collar raised', 'COLLAR_RAISE'), ('hypocotyl long', 'HYPOCOTYL_LEN'),
             ('path / depth', 'LRP_PM')]
TWO_SIDED = ['LRP', 'PM', 'TRL', 'NRL', 'ANGO2', 'DRP']


def load_batch(res_dir):
    recs = {}
    path = os.path.join(res_dir, 'checkpoint_traits.jsonl')
    if not os.path.exists(path):
        return []
    with open(path, encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                r = upgrade_record(json.loads(line))
            except Exception:
                continue
            recs[r['name']] = r          # the last record of a sample wins
    rows = []
    for name, r in recs.items():
        T = dict(r.get('T') or {})
        nb, nr = r.get('n_raw'), r.get('n_removed')
        T['%removed'] = 100.0 * nr / nb if nb else None
        lrp, pm = T.get('LRP'), T.get('PM')
        T['LRP_PM'] = lrp / pm if (isinstance(lrp, (int, float)) and isinstance(pm, (int, float)) and pm > 0) else None
        rows.append(dict(name=name, T=T))
    return rows


def robust_z(values):
    v = np.array([np.nan if x is None else float(x) for x in values], float)
    ok = np.isfinite(v)
    z = np.full(len(v), np.nan)
    if ok.sum() < 5:
        return z
    med = np.median(v[ok]); mad = 1.4826 * np.median(np.abs(v[ok] - med))
    if mad <= 0:          # more than half the values equal: use the mean absolute deviation
        mad = 1.2533 * np.mean(np.abs(v[ok] - med))
    if mad <= 0:
        return z
    z[ok] = (v[ok] - med) / mad
    return z


def rank_batch(rows, k=3.0):
    out = []
    zs = {}
    for _, key in HIGH_ONLY:
        zs[key] = robust_z([r['T'].get(key) for r in rows])
    for key in TWO_SIDED:
        zs[key] = robust_z([r['T'].get(key) for r in rows])
    for i, r in enumerate(rows):
        T = r['T']; flags = []; zmax = 0.0
        pr = T.get('PIVOT_RETURN')
        if isinstance(pr, (int, float)) and pr > PIVOT_FLAG_MM:
            flags.append(f'pivot return {pr:.1f} mm')
        for label, key in HIGH_ONLY:
            z = zs[key][i]
            if np.isfinite(z) and z > k:
                flags.append(f'{label} (z={z:.1f})'); zmax = max(zmax, z)
        for key in TWO_SIDED:
            z = zs[key][i]
            if np.isfinite(z) and abs(z) > k:
                flags.append(f'{key} {"high" if z > 0 else "low"} (z={z:.1f})'); zmax = max(zmax, abs(z))
        out.append(dict(ID=r['name'], n_flags=len(flags), max_z=round(zmax, 2), flags='; '.join(flags),
                        pivot_return=T.get('PIVOT_RETURN'), pct_removed=T.get('%removed'),
                        n_rescued=T.get('N_RESCUED'), collar_raise=T.get('COLLAR_RAISE'),
                        hypocotyl=T.get('HYPOCOTYL_LEN'), LRP=T.get('LRP'), PM=T.get('PM')))
    out.sort(key=lambda d: (-d['n_flags'], -d['max_z'], d['ID']))
    return out


def run(results_root, k=3.0, out_csv=None):
    out_csv = out_csv or os.path.join(results_root, 'qc_ranking.csv')
    allrows = []
    for res_dir in sorted(glob.glob(os.path.join(results_root, '*'))):
        if not os.path.isdir(res_dir):
            continue
        rows = load_batch(res_dir)
        if not rows:
            continue
        ranked = rank_batch(rows, k)
        for rank, d in enumerate(ranked, 1):
            allrows.append(dict(batch=os.path.basename(res_dir), rank=rank, **d))
        nflag = sum(1 for d in ranked if d['n_flags'])
        print(f'  {os.path.basename(res_dir)}: {len(ranked)} samples, {nflag} flagged')
    if not allrows:
        print(f'No checkpoint_traits.jsonl found in {results_root}/<batch>/.')
        return None
    fields = list(allrows[0].keys())
    with open(out_csv, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for d in allrows:
            w.writerow({k_: (round(v, 3) if isinstance(v, float) else v) for k_, v in d.items()})
    print(f'Ranking written: {out_csv}')
    return out_csv


def load_ranking(results_root):
    """{batch: {ID: (rank, flags)}} from results/qc_ranking.csv, or {} if absent."""
    path = os.path.join(results_root, 'qc_ranking.csv')
    if not os.path.exists(path):
        return {}
    out = {}
    with open(path, encoding='utf-8') as f:
        for d in csv.DictReader(f):
            out.setdefault(d['batch'], {})[d['ID']] = (int(d['rank']), d['flags'])
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('--results', default=os.path.join(PROJECT_ROOT, 'results'))
    ap.add_argument('--k', type=float, default=3.0, help='robust z threshold (default 3)')
    a = ap.parse_args()
    run(a.results, a.k)


if __name__ == '__main__':
    main()
