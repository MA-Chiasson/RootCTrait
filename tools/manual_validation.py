"""manual_validation.py
Validation of RootCTrait against manual tracing on real CT masks.

Two steps, run from the project root:

1. export: writes, for each selected sample, the root mask as a point cloud (.ply, in mm)
   to be traced by hand in CloudCompare, without looking at the RootCTrait outputs.
       python -m tools.manual_validation export --selection manual_validation_selection.csv
   Output: results/manual_validation/ply/<ID>.ply

2. compare: reads the manual measurements and compares them with RootCTrait, computed
   with the current code on the same masks.
       python -m tools.manual_validation compare --selection manual_validation_selection.csv
   Input, in results/manual_validation/manual/:
       lengths.csv          columns ID, primary_length_mm (length of the traced primary
                            root polyline, from the Properties panel of CloudCompare)
       <ID>_points.txt      points picked with the Point List Picking tool and exported as
                            ASCII, in this order: (1) collar, (2) deepest point of the root
                            system, then, for every first order lateral longer than 1 cm,
                            its insertion on the primary root followed by its tip
   Output, in results/manual_validation/:
       comparison_samples.csv, comparison_laterals.csv, comparison_summary.csv,
       manual_validation.png / .pdf

Coordinates are in mm in the frame of the input volume: x and y are the horizontal axes
X and Z, z is the depth (downward), as in the RSML export, so that manual and automatic
points can be matched directly.

The selection file has the columns batch and ID (sample identifier, as in the trait
table); batch must be defined in params.txt. Masks of a batch with several mask versions
use the variable recorded in RESULTS_ROOT/<batch>/selection_versions.csv.
"""
import argparse
import csv
import os
import re
import sys

import numpy as np
import scipy.ndimage as ndimage

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)
import run_pipeline as rp                                    # noqa: E402
from rootctrait.io_volume import load_volume                 # noqa: E402

OUT = os.path.join(rp.RESULTS_ROOT, 'manual_validation')
MIN_LAT_MM = 10.0        # first order laterals longer than this (chord, mm) are compared
MATCH_BASE_MM = 4.0      # insertion points closer than this can be matched
MATCH_TIP_MM = 8.0       # and their tips closer than this


# ------------------------------------------------------------------ helpers
def read_selection(path):
    with open(path, newline='', encoding='utf-8') as f:
        return [(r['batch'], r['ID']) for r in csv.DictReader(f)]


def _batch(name):
    defs = {b['name']: b for b in rp.BATCH_DEFS}
    if name not in defs:
        sys.exit(f"Batch '{name}' is not defined in params.txt")
    return defs[name]


def _sample_name(batch, ID):
    """Sample name as used by the batch pattern (e.g. ID B3T2S55 with pattern
    B3T2{name}.mat gives S55)."""
    pre, _, post = batch['pattern'].partition('{name}')
    stem = os.path.splitext(post)[0]
    name = ID
    if pre and name.startswith(pre):
        name = name[len(pre):]
    if stem and name.endswith(stem):
        name = name[:-len(stem)]
    return name


def _var_name(batch_name, sample, ID):
    path = os.path.join(rp.RESULTS_ROOT, batch_name, 'selection_versions.csv')
    if not os.path.exists(path):
        return None
    with open(path, newline='') as f:
        sel = {r['sample']: r['variable_used'] for r in csv.DictReader(f)}
    v = sel.get(sample, sel.get(ID))
    return None if v in (None, '', '(largest)') else v


def load_mask(batch_name, ID):
    b = _batch(batch_name)
    sample = _sample_name(b, ID)
    path = os.path.join(rp.DATA_ROOT, batch_name, b['pattern'].format(name=sample))
    V = load_volume(path, axis_order=b['axis_order'],
                    var_name=b.get('var_name') or _var_name(batch_name, sample, ID))
    return V > 0.5 * float(np.max(V))


def write_ply(path, P):
    """Binary little endian PLY of float32 points."""
    P = np.asarray(P, np.float32)
    with open(path, 'wb') as f:
        f.write((f"ply\nformat binary_little_endian 1.0\nelement vertex {len(P)}\n"
                 "property float x\nproperty float y\nproperty float z\nend_header\n").encode())
        f.write(P.tobytes())


def to_mm(vox, vs):
    """Voxel indices (depth, X, Z) to mm points (x = X, y = Z, z = depth)."""
    v = np.asarray(vox, float) * vs
    return np.c_[v[:, 1], v[:, 2], v[:, 0]]


def read_points(path):
    pts = []
    with open(path, encoding='utf-8', errors='ignore') as f:
        for line in f:
            nums = re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', line)
            if len(nums) >= 3:
                pts.append([float(x) for x in nums[-3:]])
    return np.array(pts)


def angle_vertical(a, b):
    d = np.asarray(b, float) - np.asarray(a, float)
    return float(np.degrees(np.arccos(min(1.0, abs(d[2]) / (np.linalg.norm(d) + 1e-12)))))


# ------------------------------------------------------------------ export
def cmd_export(args):
    vs = np.asarray(rp.VOXEL_SIZE, float)
    os.makedirs(os.path.join(OUT, 'ply'), exist_ok=True)
    for batch_name, ID in read_selection(args.selection):
        BW = load_mask(batch_name, ID)
        if args.full:
            vox = np.argwhere(BW)
        else:                                   # surface voxels: lighter, same shape
            vox = np.argwhere(BW & ~ndimage.binary_erosion(BW))
        out = os.path.join(OUT, 'ply', f'{ID}.ply')
        write_ply(out, to_mm(vox, vs))
        print(f'  {ID}: {len(vox)} points -> {out}', flush=True)
    os.makedirs(os.path.join(OUT, 'manual'), exist_ok=True)
    tmpl = os.path.join(OUT, 'manual', 'lengths.csv')
    if not os.path.exists(tmpl):
        with open(tmpl, 'w', newline='') as f:
            w = csv.writer(f); w.writerow(['ID', 'primary_length_mm'])
            for _, ID in read_selection(args.selection):
                w.writerow([ID, ''])
        print(f'  template -> {tmpl}')


# ------------------------------------------------------------------ compare
def rootctrait_side(batch_name, ID):
    """Traits and first order laterals computed by RootCTrait on the mask."""
    BW = load_mask(batch_name, ID)
    R = rp.analyse_mask(BW, rp._effective_params())
    vs = np.asarray(rp.VOXEL_SIZE, float)
    off = np.array([s.start for s in R['crop']])
    roots = R['roots']
    order1 = {s['seg_id'] for s in roots if s['order'] == 1}
    lats = []
    for s in roots:
        if s['order'] == 2 and s.get('parent_seg') in order1:
            P = to_mm(s['coords'] + off, vs)
            if np.linalg.norm(P[-1] - P[0]) >= MIN_LAT_MM:
                lats.append((P[0], P[-1]))
    return R['T'], lats


def match(manual, auto):
    """Greedy matching of laterals on insertion and tip distance."""
    pairs, used = [], set()
    cand = []
    for i, (mb, mt) in enumerate(manual):
        for j, (ab, at) in enumerate(auto):
            db, dt = np.linalg.norm(mb - ab), np.linalg.norm(mt - at)
            if db <= MATCH_BASE_MM and dt <= MATCH_TIP_MM:
                cand.append((db + 0.5 * dt, i, j))
    mi = set()
    for _, i, j in sorted(cand):
        if i in mi or j in used:
            continue
        mi.add(i); used.add(j); pairs.append((i, j))
    return pairs


def _stats(m, a):
    m, a = np.asarray(m, float), np.asarray(a, float)
    ok = np.isfinite(m) & np.isfinite(a)
    m, a = m[ok], a[ok]
    if len(m) < 2:
        return dict(n=len(m))
    r = np.corrcoef(m, a)[0, 1]
    return dict(n=len(m), r2=r * r, mae=float(np.mean(np.abs(a - m))),
                mape=float(np.mean(np.abs(a - m) / np.abs(m)) * 100), bias=float(np.mean(a - m)))


def cmd_compare(args):
    global OUT
    man_dir = args.manual_dir or os.path.join(OUT, 'manual')
    if args.manual_dir:                      # second observer: results written next to its files
        OUT = args.manual_dir
    lengths = {}
    with open(os.path.join(man_dir, 'lengths.csv'), newline='', encoding='utf-8') as f:
        for r in csv.DictReader(f):
            if r['primary_length_mm'].strip():
                lengths[r['ID']] = float(r['primary_length_mm'].replace(',', '.'))
    srows, lrows = [], []
    TP = FP = FN = 0
    for batch_name, ID in read_selection(args.selection):
        pf = os.path.join(man_dir, f'{ID}_points.txt')
        if not os.path.exists(pf) or ID not in lengths:
            print(f'  {ID}: manual data missing, skipped'); continue
        P = read_points(pf)
        if len(P) < 2 or (len(P) - 2) % 2:
            print(f'  {ID}: {len(P)} points; expected collar, deepest point, then pairs (insertion, tip). Skipped')
            continue
        collar, deep = P[0], P[1]
        man_lat = [(P[k], P[k + 1]) for k in range(2, len(P), 2)
                   if np.linalg.norm(P[k + 1] - P[k]) >= MIN_LAT_MM]
        T, auto_lat = rootctrait_side(batch_name, ID)
        pairs = match(man_lat, auto_lat)
        tp = len(pairs); fp = len(auto_lat) - tp; fn = len(man_lat) - tp
        TP += tp; FP += fp; FN += fn
        for i, j in pairs:
            lrows.append(dict(ID=ID, angle_manual=angle_vertical(*man_lat[i]),
                              angle_rootctrait=angle_vertical(*auto_lat[j]),
                              insertion_distance_mm=float(np.linalg.norm(man_lat[i][0] - auto_lat[j][0]))))
        ang_m = np.mean([angle_vertical(*l) for l in man_lat]) if man_lat else np.nan
        ang_a = np.mean([angle_vertical(*l) for l in auto_lat]) if auto_lat else np.nan
        srows.append(dict(batch=batch_name, ID=ID,
                          PRL_manual=lengths[ID], PRL_rootctrait=T['PRL'],
                          MD_manual=float(deep[2] - collar[2]), MD_rootctrait=T['MD'],
                          laterals_manual=len(man_lat), laterals_rootctrait=len(auto_lat), matched=tp,
                          angle_manual=ang_m, angle_rootctrait=ang_a))
        print(f'  {ID}: PRL {lengths[ID]:.1f} / {T["PRL"]:.1f} mm, laterals {len(man_lat)} / {len(auto_lat)} '
              f'({tp} matched)', flush=True)
    if not srows:
        sys.exit('No sample with complete manual data.')
    os.makedirs(OUT, exist_ok=True)
    for name, rows in (('comparison_samples.csv', srows), ('comparison_laterals.csv', lrows)):
        if rows:
            with open(os.path.join(OUT, name), 'w', newline='') as f:
                w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    col = lambda rows, k: [r[k] for r in rows]
    summ = []
    for lab, a, b in (('PRL (mm)', 'PRL_manual', 'PRL_rootctrait'), ('MD (mm)', 'MD_manual', 'MD_rootctrait'),
                      ('First order laterals > 1 cm (n)', 'laterals_manual', 'laterals_rootctrait'),
                      ('Mean lateral angle (deg)', 'angle_manual', 'angle_rootctrait')):
        summ.append(dict(measure=lab, **_stats(col(srows, a), col(srows, b))))
    if lrows:
        summ.append(dict(measure='Angle of matched laterals (deg)',
                         **_stats(col(lrows, 'angle_manual'), col(lrows, 'angle_rootctrait'))))
    prec = TP / (TP + FP) if TP + FP else np.nan
    rec = TP / (TP + FN) if TP + FN else np.nan
    f1 = 2 * prec * rec / (prec + rec) if prec + rec else np.nan
    summ.append(dict(measure='Lateral detection', n=TP + FN, precision=prec, recall=rec, f1=f1))
    keys = []
    for r in summ:
        keys += [k for k in r if k not in keys]
    with open(os.path.join(OUT, 'comparison_summary.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(summ)
    for r in summ:
        print('  ' + ', '.join(f'{k}={v:.3g}' if isinstance(v, float) else f'{k}={v}' for k, v in r.items()))
    _figure(srows, lrows, summ)
    print(f'-> {OUT}')


def _figure(srows, lrows, summ):
    try:
        import matplotlib; matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        return
    plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Liberation Sans', 'DejaVu Sans'],
                         'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False})
    panels = [('PRL_manual', 'PRL_rootctrait', 'Primary root length (mm)'),
              ('MD_manual', 'MD_rootctrait', 'Rooting depth (mm)'),
              ('laterals_manual', 'laterals_rootctrait', 'First order laterals (n)')]
    fig, axs = plt.subplots(1, 4, figsize=(7.2, 2.1))
    for ax, (a, b, t), s, L in zip(axs, panels, summ, 'abc'):
        x = np.array([r[a] for r in srows], float); y = np.array([r[b] for r in srows], float)
        lo, hi = np.nanmin([x, y]), np.nanmax([x, y]); pad = 0.08 * (hi - lo + 1e-9)
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color='#9AA0A6', lw=0.7, ls=(0, (3, 2)))
        ax.scatter(x, y, s=14, color='#48627A', edgecolor='white', lw=0.5)
        ax.set_xlabel('Manual'); ax.set_ylabel('RootCTrait'); ax.set_title(t, fontsize=8, loc='left')
        ax.text(0.97, 0.04, f"R² = {s.get('r2', np.nan):.2f}\nMAE = {s.get('mae', np.nan):.1f}",
                transform=ax.transAxes, ha='right', va='bottom', fontsize=6.5)
        ax.text(-0.02, 1.02, L, transform=ax.transAxes, fontsize=10, fontweight='bold', ha='right', va='bottom')
    ax = axs[3]
    if lrows:
        x = np.array([r['angle_manual'] for r in lrows]); y = np.array([r['angle_rootctrait'] for r in lrows])
        ax.plot([0, 90], [0, 90], color='#9AA0A6', lw=0.7, ls=(0, (3, 2)))
        ax.scatter(x, y, s=8, color='#D48A48', edgecolor='white', lw=0.3)
        s = summ[4] if len(summ) > 5 else {}
        ax.text(0.97, 0.04, f"MAE = {s.get('mae', np.nan):.1f}°", transform=ax.transAxes, ha='right', va='bottom', fontsize=6.5)
    ax.set_xlabel('Manual'); ax.set_ylabel('RootCTrait'); ax.set_title('Lateral angle (deg)', fontsize=8, loc='left')
    ax.text(-0.02, 1.02, 'd', transform=ax.transAxes, fontsize=10, fontweight='bold', ha='right', va='bottom')
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(OUT, f'manual_validation.{ext}'), dpi=300)


def main():
    ap = argparse.ArgumentParser(description='Validation of RootCTrait against manual tracing.')
    sub = ap.add_subparsers(dest='cmd', required=True)
    e = sub.add_parser('export'); e.add_argument('--selection', required=True)
    e.add_argument('--full', action='store_true', help='all mask voxels instead of the surface only')
    c = sub.add_parser('compare'); c.add_argument('--selection', required=True)
    c.add_argument('--manual-dir', default=None,
                   help='folder of the manual files (default results/manual_validation/manual); '
                        'use another folder for a second observer')
    a = ap.parse_args()
    cmd_export(a) if a.cmd == 'export' else cmd_compare(a)


if __name__ == '__main__':
    main()
