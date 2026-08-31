"""Decontamination of false roots.

Two geometric steps, uniform across all samples (GWAS compatible):

1. Sheet rule: a segment is removed if it is simultaneously
     - parallel     : bc >= BC_MIN  (at least BC_MIN parallel, offset neighbors)
     - low linearity: lin < LIN_MAX (sheet-like neighborhood rather than a line)
     - short        : length < LEN_MAX mm

2. Floating fragments: after step 1, only the segments still connected to the
   collar are kept (connected component containing the base). Segments that were
   attached to the system only through removed pollution become orphans and are
   discarded in turn.

bc  : number of neighbors (centroid within rpar) that are parallel (|cos| > 0.9)
      and laterally offset (|cos of the offset direction| < 0.5).
lin : linearity (l1 - l2) / l1 of the PCA of the skeleton neighborhood around the segment.

Axis convention: col0 = Y (depth), col1 = X, col2 = Z.
"""
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

BC_MIN = 3
LIN_MAX = 0.7
LEN_MAX = 15.0
RPAR = 6.0
RNB = 6.0


def _seg_dir(coords_phys):
    if len(coords_phys) < 2:
        return np.array([1.0, 0.0, 0.0])
    X = coords_phys - coords_phys.mean(0)
    d = np.linalg.svd(X, full_matrices=False)[2][0]
    return d / (np.linalg.norm(d) + 1e-9)


def segment_features(segments, skv_full, voxel_size, rpar=RPAR, rnb=RNB):
    """Return bc, lin and length per segment."""
    vs = np.asarray(voxel_size, float)
    n = len(segments)
    length = np.array([s['length_mm'] for s in segments], float)
    if n == 0:
        return dict(bc=np.zeros(0, int), lin=np.zeros(0), length=length)
    cents = np.array([(s['coords'] * vs).mean(0) for s in segments])
    dirs = np.array([_seg_dir(s['coords'].astype(float) * vs) for s in segments])

    bc = np.zeros(n, int)
    ct = cKDTree(cents)
    for i in range(n):
        cnt = 0
        for j in ct.query_ball_point(cents[i], rpar):
            if j == i:
                continue
            if abs(dirs[j] @ dirs[i]) < 0.9:
                continue
            o = cents[j] - cents[i]
            no = np.linalg.norm(o)
            if no > 1e-6 and abs((o / no) @ dirs[i]) < 0.5:
                cnt += 1
        bc[i] = cnt

    lin = np.zeros(n)
    if len(skv_full):
        P = skv_full.astype(float) * vs
        vt = cKDTree(P)
        for i in range(n):
            nb = np.asarray(vt.query_ball_point(cents[i], rnb))
            if len(nb) >= 3:
                Q = P[nb] - P[nb].mean(0)
                ev = np.linalg.svd(Q, full_matrices=False)[1] ** 2
                ev = ev / ev.sum()
                lin[i] = (ev[0] - ev[1]) / (ev[0] + 1e-9)
    return dict(bc=bc, lin=lin, length=length)


def keep_base_component(segments, base):
    """Keep only the segments connected to the collar (connected component of the base).
    Two segments are connected if their voxels touch (distance <= sqrt(3)).
    Returns (connected_segments, orphan_segments)."""
    n = len(segments)
    if n <= 1:
        return list(segments), []
    allvox = np.vstack([s['coords'] for s in segments])
    owner = np.concatenate([[i] * len(s['coords']) for i, s in enumerate(segments)])
    pairs = cKDTree(allvox).query_pairs(r=1.7321, output_type='ndarray')
    if len(pairs):
        oi = owner[pairs[:, 0]]; oj = owner[pairs[:, 1]]; m = oi != oj
        A = csr_matrix((np.ones(int(m.sum()) * 2),
                        (np.r_[oi[m], oj[m]], np.r_[oj[m], oi[m]])), shape=(n, n))
    else:
        A = csr_matrix((n, n))
    _, lab = connected_components(A, directed=False)
    # main component = the one containing the pivot (order 1).
    # Fallback to the component closest to the collar if no order 1 exists.
    order1 = [i for i, s in enumerate(segments) if s.get('order') == 1]
    if order1:
        main = lab[order1[0]]
    else:
        cents = np.array([s['coords'].mean(0) for s in segments])
        main = lab[np.argmin(np.linalg.norm(cents - np.asarray(base), axis=1))]
    kept = [s for s, l in zip(segments, lab) if l == main]
    orphan = [s for s, l in zip(segments, lab) if l != main]
    return kept, orphan


def decontaminate(segments, skv_full, voxel_size, base=None,
                  bc_min=BC_MIN, lin_max=LIN_MAX, len_max=LEN_MAX, drop_orphans=True):
    """Remove the sheets, then the floating fragments.
    Returns (kept_segments, removed_segments, features)."""
    f = segment_features(segments, skv_full, voxel_size)
    rule_removed = (f['bc'] >= bc_min) & (f['lin'] < lin_max) & (f['length'] < len_max)
    kept = [s for s, r in zip(segments, rule_removed) if not r]
    removed = [s for s, r in zip(segments, rule_removed) if r]
    if drop_orphans and base is not None and len(kept) > 1:
        kept, orphan = keep_base_component(kept, base)
        removed += orphan
    return kept, removed, f
