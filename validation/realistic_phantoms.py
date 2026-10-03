"""realistic_phantoms.py

Synthetic root systems of known geometry, built at four levels of increasing
difficulty, to test the steps of RootCTrait that the simple phantoms of
validate_phantoms.py do not exercise.

  L1 shape      curved, tapering laterals with order 3 branches, thick primary base
  L2 contacts   L1 + laterals crossing other laterals, primary root turning back up
  L3 artifacts  L2 + real surface sheets transplanted from a scan, isolated fragments,
                rough mask surface
  L4 stem       L3 + hypocotyl column above the collar with a high branch

Levels are cumulative. For each phantom the true value of every tested trait is
computed from the centre lines used to build it, and the masks of the real roots,
of the artifacts and of the hypocotyl are kept separately, so that decontamination
and hypocotyl exclusion can be scored voxel by voxel on the skeleton.

Surface sheets are transplanted from the real scan B2T2S11 (validation/data).

Run: python validation/realistic_phantoms.py [--n 5] [--out validation/realistic]
"""
import os, sys, csv, argparse, json
import numpy as np
import scipy.ndimage as ndimage
from scipy.spatial import cKDTree
from skimage.morphology import skeletonize

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from rootctrait.graph_extraction import prune_skeleton
from rootctrait.root_decomposition import decompose_root_system
from rootctrait.decontamination import decontaminate, keep_base_component
from rootctrait.root_traits_full import compute_all_traits, root_mask
from rootctrait.detection_hypocotyle import collar_and_hypocotyl
import run_pipeline as rp

VS = np.array([0.39, 0.39, 0.2])          # mm, (depth, X, Z)
BOX = np.array([100.0, 70.0, 70.0])       # mm
SHAPE = tuple(int(b / v) for b, v in zip(BOX, VS))
CENTRE = np.array([0.0, BOX[1] / 2, BOX[2] / 2])
DS = 0.25                                  # centre line step, mm


# ---------------------------------------------------------------- drawing
def draw_tube(V, P, R):
    """Union of capsules between successive centre line points P (mm), radii R (mm)."""
    for a, b, ra, rb in zip(P[:-1], P[1:], R[:-1], R[1:]):
        r = max(ra, rb)
        lo = np.maximum(np.floor((np.minimum(a, b) - r) / VS).astype(int) - 1, 0)
        hi = np.minimum(np.ceil((np.maximum(a, b) + r) / VS).astype(int) + 1, np.array(V.shape) - 1)
        if np.any(hi < lo):
            continue
        g = np.mgrid[lo[0]:hi[0] + 1, lo[1]:hi[1] + 1, lo[2]:hi[2] + 1].reshape(3, -1).T
        X = (g + 0.5) * VS
        d = b - a; L2 = d @ d
        t = np.clip(((X - a) @ d) / (L2 + 1e-12), 0, 1)
        dist = np.linalg.norm(X - (a + t[:, None] * d), axis=1)
        rr = ra + t * (rb - ra)
        m = dist <= rr
        g = g[m]
        V[g[:, 0], g[:, 1], g[:, 2]] = True


def unit(v):
    return v / (np.linalg.norm(v) + 1e-12)


def rotate_towards(d, target, ang):
    """Rotate unit vector d towards target by at most ang radians."""
    c = np.clip(d @ target, -1, 1); th = np.arccos(c)
    if th < 1e-9:
        return d
    k = min(ang, th) / th
    return unit(d * (1 - k) + target * k)


def path_len(P):
    return float(np.linalg.norm(np.diff(P, axis=0), axis=1).sum())


def angle_vertical(p0, p1):
    v = p1 - p0
    return float(np.degrees(np.arccos(np.clip(abs(v[0]) / (np.linalg.norm(v) + 1e-12), 0, 1))))


DOWN = np.array([1.0, 0.0, 0.0])


# ---------------------------------------------------------------- phantom
def build(level, seed):
    rng = np.random.default_rng(seed)
    V_root = np.zeros(SHAPE, bool); V_art = np.zeros(SHAPE, bool); V_hyp = np.zeros(SHAPE, bool)
    collar = CENTRE + np.array([14.0, rng.uniform(-3, 3), rng.uniform(-3, 3)])
    roots = []                                   # dicts: P, R, order

    # primary root: thick base plateau, gentle wander, taper
    Lp = rng.uniform(62, 78)
    n = int(Lp / DS); P = [collar.copy()]; d = DOWN.copy()
    ph = rng.uniform(0, 2 * np.pi, 2)
    for i in range(1, n):
        s = i * DS
        wob = 0.18 * np.array([0, np.sin(s / 9 + ph[0]), np.cos(s / 11 + ph[1])])
        d = unit(DOWN + wob)
        P.append(P[-1] + d * DS)
    P = np.array(P)
    hook = 0.0
    if level >= 2 and seed % 2 == 0:             # primary turns back up at the pot bottom
        d = unit(P[-1] - P[-2]); turn = []
        side = unit(np.array([0, rng.uniform(-1, 1), rng.uniform(-1, 1)]))
        q = P[-1].copy()
        for k in range(int(4.5 / DS)):           # half turn, radius about 1.5 mm
            d = rotate_towards(d, side, np.radians(180 * DS / 4.5))
            q = q + d * DS; turn.append(q.copy())
        for k in range(int(6.0 / DS)):           # then straight up, slightly inclined
            d = rotate_towards(d, -DOWN, np.radians(10))
            q = q + d * DS; turn.append(q.copy())
        P = np.vstack([P, turn])
        hook = float(P[:, 0].max() - P[-1, 0])
    s = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))])
    R = np.where(s < 6.0, 1.6, 1.6 - (s - 6.0) / (s[-1] - 6.0) * 1.0)
    roots.append(dict(P=P, R=R, order=1, nchild=0))

    # laterals
    nl = int(rng.integers(8, 13))
    pos = np.sort(rng.uniform(3.0, 0.82 * Lp, nl)); pos[0] = rng.uniform(2.0, 3.5)  # one basal root
    phis = rng.uniform(0, 2 * np.pi, nl)
    lat_idx = []
    for j, (sp, phi) in enumerate(zip(pos, phis)):
        i0 = int(np.searchsorted(s, sp)); p0 = P[i0]
        a0 = np.radians(rng.uniform(80, 88) if j == 0 else rng.uniform(55, 85))
        d = unit(np.array([np.cos(a0), np.sin(a0) * np.cos(phi), np.sin(a0) * np.sin(phi)]))
        L = rng.uniform(12, 30) * (1 - 0.4 * sp / Lp)
        kappa = np.radians(rng.uniform(0.5, 1.5)) * DS
        Q = [p0.copy()]
        for k in range(int(L / DS)):
            d = rotate_towards(d, DOWN, kappa)
            d = unit(d + rng.normal(0, 0.04, 3) * np.array([0.3, 1, 1]))
            Q.append(Q[-1] + d * DS)
        Q = np.array(Q)
        sq = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(Q, axis=0), axis=1))])
        roots.append(dict(P=Q, R=0.55 - 0.15 * sq / sq[-1], order=2, nchild=0)); lat_idx.append(len(roots) - 1); roots[0]['nchild'] += 1
    # order 3 branches on the longer laterals
    for li in list(lat_idx):
        Q = roots[li]['P']; Lq = path_len(Q)
        if Lq < 16:
            continue
        for f in rng.uniform(0.3, 0.65, int(rng.integers(1, 3))):
            i0 = int(f * (len(Q) - 1)); p0 = Q[i0]
            t = unit(Q[min(i0 + 4, len(Q) - 1)] - Q[max(i0 - 4, 0)])
            side = unit(np.cross(t, rng.normal(0, 1, 3)))
            side[0] = -abs(side[0]) * 0.3                         # never deeper than the parent tip
            d = unit(0.5 * t + side)
            L = rng.uniform(3.5, 7.0)
            B = np.array([p0 + d * DS * k for k in range(int(L / DS) + 1)])
            roots.append(dict(P=B, R=np.full(len(B), 0.42), order=3, nchild=0)); roots[li]['nchild'] += 1
    # crossings: laterals aimed through the middle of another lateral
    if level >= 2:
        for li in rng.choice(lat_idx[1:], 2, replace=False):
            A = roots[li]['P']; mid = A[len(A) // 2]
            i0 = int(np.searchsorted(s, min(s[-1] * 0.85, roots[li]['P'][0, 0] - collar[0] + 6)))
            p0 = P[i0]; d = unit(mid - p0)
            L = np.linalg.norm(mid - p0) + rng.uniform(5, 10)
            C = np.array([p0 + d * DS * k for k in range(int(L / DS) + 1)])
            roots.append(dict(P=C, R=np.full(len(C), 0.5), order=2, nchild=0)); roots[0]['nchild'] += 1
    for r in roots:
        draw_tube(V_root, r['P'], r['R'])

    # artifacts
    if level >= 3:
        # surface sheets transplanted from a real scan (B2T2S11): the mask voxels that
        # the sheet rule removed there, pasted at the same voxel size, mirrored at random,
        # and placed so that they touch the root system as in the real data
        sheets = np.load(os.path.join(ROOT, 'validation', 'data', 'real_sheets_B2T2S11.npz'))
        shallow = sorted(lat_idx, key=lambda i: roots[i]['P'][0, 0])
        anchors = [roots[shallow[0]]['P'][int(0.5 * (len(roots[shallow[0]]['P']) - 1))],   # large sheet on a basal root
                   collar + np.array([1.0, rng.uniform(-3, 3), rng.uniform(-3, 3)])]       # small sheet at the collar
        for key, anc in zip(('arr_0', 'arr_1'), anchors):
            G = sheets[key].copy()
            for ax in (1, 2):
                if rng.random() < 0.5:
                    G[:, ax] = G[:, ax].max() - G[:, ax]
            cen = G.mean(0)
            g = np.round(G - cen + anc / VS).astype(int)
            g = g[(g >= 0).all(1) & (g < np.array(SHAPE)).all(1)]
            V_art[g[:, 0], g[:, 1], g[:, 2]] = True
        for k in range(6):                       # isolated fragments
            for _ in range(50):
                c = CENTRE + np.array([rng.uniform(15, 80), rng.uniform(-28, 28), rng.uniform(-28, 28)])
                if not V_root[tuple((c / VS).astype(int))]:
                    break
            r = rng.uniform(0.6, 1.0)
            draw_tube(V_art, np.array([c, c + unit(rng.normal(0, 1, 3)) * rng.uniform(1, 3)]), np.array([r, r]))
    if level >= 4:                               # hypocotyl column with a high branch
        top = collar + np.array([-10.0, rng.uniform(-1, 1), rng.uniform(-1, 1)])
        H = np.array([collar + (top - collar) * t for t in np.linspace(0, 1, int(10 / DS))])
        r_h = rng.uniform(1.15, 1.45)
        draw_tube(V_hyp, H, np.full(len(H), r_h))
        b0 = collar + (top - collar) * 0.55; d = unit(np.array([-0.3, rng.uniform(-1, 1), rng.uniform(-1, 1)]))
        Bh = np.array([b0 + d * DS * k for k in range(int(6 / DS) + 1)])
        draw_tube(V_hyp, Bh, np.full(len(Bh), 0.6))
    labels = np.zeros(SHAPE, np.uint8)           # 1 root, 2 artifact, 3 hypocotyl
    labels[V_art] = 2; labels[V_hyp] = 3; labels[V_root] = 1
    if level >= 3:                               # rough surface: random one voxel bumps
        V = labels > 0
        surf = np.argwhere(V & ~ndimage.binary_erosion(V))
        pick = surf[rng.random(len(surf)) < 0.04]
        q = np.clip(pick + rng.integers(-1, 2, pick.shape), 0, np.array(SHAPE) - 1)
        free = labels[q[:, 0], q[:, 1], q[:, 2]] == 0
        labels[q[free, 0], q[free, 1], q[free, 2]] = labels[pick[free, 0], pick[free, 1], pick[free, 2]]
    V = labels > 0

    # truth, computed on the centre lines of the real roots
    lat = [r for r in roots if r['order'] >= 2]
    o2 = [r for r in roots if r['order'] == 2]
    allP = np.vstack([r['P'] for r in roots])
    T = dict(LRP=path_len(roots[0]['P']), TRL=sum(path_len(r['P']) for r in roots),
             LTRL=sum(path_len(r['P']) for r in lat), NRL=len(lat), NT=sum(1 for r in roots if r['nchild'] == 0),
             MaxO=max(r['order'] for r in roots),
             MLRL=float(np.mean([path_len(r['P']) for r in lat])),
             ANGO2=float(np.mean([angle_vertical(r['P'][0], r['P'][-1]) for r in o2])),
             PM=float(allP[:, 0].max() - collar[0]),
             DRP=float(np.mean(2 * roots[0]['R'])),
             DRS=float(np.mean([2 * np.median(r['R']) for r in lat])),
             TOR=float(np.mean([path_len(r['P']) / np.linalg.norm(r['P'][-1] - r['P'][0]) for r in lat])),
             VRT=float((labels == 1).sum() * VS.prod()), PIVOT_RETURN=hook)
    return V, labels, T, collar


# ---------------------------------------------------------------- pipeline (run_pipeline.analyse_mask)
CTX = dict(voxel_size=tuple(VS), prune_vox=5, min_seg_len_mm=2.0, bc_min=3, lin_max=0.7, len_max=15,
           drop_orphans=True, dens_max=35, dens_len_max=6.0, rescue_min_mm=5.0)


def run(V, **over):
    R = rp.analyse_mask(V, dict(CTX, **over))
    return dict(T=R['T'], off=np.array([s.start for s in R['crop']]), base2=np.asarray(R['base2']),
                segs=R['segs'], removed=R['removed'], hypo=R['hypocotyl'], roots=R['roots'],
                n_rescued=R['n_rescued'])


def seg_points(segs, off):
    return np.vstack([s['coords'] for s in segs]) + off if segs else np.zeros((0, 3), int)


TRAITS = ['LRP', 'TRL', 'LTRL', 'NRL', 'NT', 'MaxO', 'MLRL', 'ANGO2', 'PM', 'DRP', 'DRS', 'TOR', 'VRT']
TO_MM = dict(VRT=1.0)                            # compute_all_traits works in mm and mm3


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--n', type=int, default=5)
    ap.add_argument('--out', default=os.path.join(ROOT, 'validation', 'realistic'))
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    rows, qc = [], []
    for level in (1, 2, 3, 4):
        for k in range(a.n):
            seed = 100 * level + k
            V, lab, truth, collar = build(level, seed)
            out = run(V)
            T = out['T']
            for tr in TRAITS:
                t = truth[tr]; mval = T.get(tr)
                rows.append(dict(level=level, phantom=seed, trait=tr, true=t, measured=mval,
                                 err_pct=(100 * (mval - t) / t) if t else np.nan))
            # skeleton based scores
            def lab_of(P):
                return lab[P[:, 0], P[:, 1], P[:, 2]] if len(P) else np.zeros(0, int)
            # scores on the voxels of the initial skeleton tree, each labelled by the
            # object it lies in; a voxel counts as kept when it is in the final roots
            allp = np.unique(seg_points(out['segs'], out['off']), axis=0); la = lab_of(allp)
            kept = np.unique(seg_points(out['roots'], out['off']), axis=0)
            hyp = np.unique(seg_points(out['hypo'], out['off']), axis=0)
            key = lambda P: set(map(tuple, P))
            K = key(kept); H = key(hyp)
            in_k = np.array([tuple(p) in K for p in allp]); in_h = np.array([tuple(p) in H for p in allp])
            lk = lab_of(kept)
            n_art = int((la == 2).sum()); n_hyp = int((la == 3).sum()); n_root = int((la == 1).sum())
            q = dict(level=level, phantom=seed, rescued=out['n_rescued'],
                     artifact_removed_pct=100 * ((la == 2) & ~in_k).sum() / n_art if n_art else np.nan,
                     root_removed_pct=100 * ((la == 1) & ~in_k).sum() / n_root,
                     artifact_kept_pct=100 * (lk == 2).sum() / max(1, len(lk)),
                     hypocotyl_excluded_pct=100 * ((la == 3) & in_h).sum() / n_hyp if n_hyp else np.nan,
                     hypocotyl_kept_pct=100 * ((la == 3) & in_k).sum() / n_hyp if n_hyp else np.nan,
                     collar_error_mm=float((out['base2'][0] + out['off'][0]) * VS[0] - collar[0]),
                     pivot_return_true=truth['PIVOT_RETURN'], pivot_return=T.get('PIVOT_RETURN'))
            qc.append(q)
            print(level, seed, {tr: round(rows[-len(TRAITS) + i]['err_pct'], 1) for i, tr in enumerate(TRAITS)},
                  {k_: (round(v, 2) if isinstance(v, float) else v) for k_, v in q.items() if k_ not in ('level', 'phantom')}, flush=True)
    for name, data in (('traits.csv', rows), ('cleaning.csv', qc)):
        with open(os.path.join(a.out, name), 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(data[0].keys())); w.writeheader(); w.writerows(data)
    print('\nmean absolute error (%) per level')
    for tr in TRAITS:
        print(f'{tr:>6}', ' '.join(f'{np.nanmean([abs(r["err_pct"]) for r in rows if r["trait"] == tr and r["level"] == L]):6.1f}' for L in (1, 2, 3, 4)))


if __name__ == '__main__':
    main()
