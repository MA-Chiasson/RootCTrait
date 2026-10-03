"""Decontamination of false roots.

Two geometric steps, uniform across all samples (GWAS compatible):

1. Sheet rule: a segment is removed if it is simultaneously
     - parallel     : bc >= BC_MIN  (at least BC_MIN parallel, offset neighbors)
     - low linearity: lin < LIN_MAX (sheet-like neighborhood rather than a line)
     - short        : length < LEN_MAX mm

   Dense rule: a segment is also removed if the skeleton around it is dense
   (dens >= DENS_MAX voxels within RDENS mm, a mesh rather than a line) and it is
   short (length < DENS_LEN_MAX mm). This catches the parts of surface layers whose
   skeleton is a mesh rather than a ladder of parallel lines.

2. Floating fragments: after step 1, only the segments still connected to the
   collar are kept (connected component containing the base). Segments that were
   attached to the system only through removed pollution become orphans. An orphan
   group that looks like a root (total length >= RESCUE_MIN_MM, median dens below
   DENS_MAX) is rescued: a genuine root whose only link to the system ran through a
   removed layer. The others are discarded. Rescued groups are returned in
   features['rescued'] and reattached to the system by reattach().

bc  : number of neighbors (centroid within rpar) that are parallel (|cos| > 0.9)
      and laterally offset (|cos of the offset direction| < 0.5).
lin : linearity (l1 - l2) / l1 of the PCA of the skeleton neighborhood around the segment.
dens: median, over the voxels of the segment, of the number of skeleton voxels within
      RDENS mm.

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
DENS_MAX = 35
DENS_LEN_MAX = 6.0
RESCUE_MIN_MM = 5.0
RDENS = 2.0


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
    dens = np.zeros(n)
    if n:
        U = np.unique(np.vstack([sg['coords'] for sg in segments]), axis=0).astype(float) * vs
        ut = cKDTree(U)
        for i, sgm in enumerate(segments):
            dens[i] = np.median(ut.query_ball_point(sgm['coords'] * vs, RDENS, return_length=True))
    return dict(bc=bc, lin=lin, length=length, dens=dens)


# Half of the 26-neighbourhood: each touching voxel pair is visited once.
_HALF_26 = np.array([d for d in np.ndindex(3, 3, 3)
                     if (np.array(d) - 1).tolist() > [0, 0, 0]]) - 1


def segment_adjacency(segments):
    """Sparse segment graph for connectivity. Two segments touch when one of their
    voxels is identical to, or 26-connected with, a voxel of the other (distance
    <= sqrt(3) in voxel units). The graph has the same connected components as this
    touching relation (a voxel shared by several segments links them to each other
    rather than pairwise to every neighbour). Memory grows linearly with the number
    of voxels (sorted voxel keys and binary search), so it is safe for very large
    root systems."""
    n = len(segments)
    V = np.vstack([s['coords'] for s in segments]).astype(np.int64)
    owner = np.concatenate([np.full(len(s['coords']), i, dtype=np.int64)
                            for i, s in enumerate(segments)])
    V = V - V.min(0) + 1                      # pad by one voxel on every side
    dims = tuple(int(d) for d in V.max(0) + 2)
    key = np.ravel_multi_index(V.T, dims)
    order = np.argsort(key, kind='stable')
    skey, sown = key[order], owner[order]
    rows, cols = [], []
    same = np.where(skey[1:] == skey[:-1])[0]     # one voxel shared by two segments
    rows.append(sown[same]); cols.append(sown[same + 1])
    for off in _HALF_26:
        nk = np.ravel_multi_index((V + off).T, dims)
        pos = np.minimum(np.searchsorted(skey, nk), len(skey) - 1)
        hit = skey[pos] == nk
        rows.append(owner[hit]); cols.append(sown[pos[hit]])
    r = np.concatenate(rows); c = np.concatenate(cols); m = r != c
    r, c = r[m], c[m]
    return csr_matrix((np.ones(2 * len(r)), (np.r_[r, c], np.r_[c, r])), shape=(n, n))


def keep_base_component(segments, base):
    """Keep only the segments connected to the collar (connected component of the base).
    Two segments are connected if their voxels touch (distance <= sqrt(3)).
    Returns (connected_segments, orphan_segments)."""
    n = len(segments)
    if n <= 1:
        return list(segments), []
    _, lab = connected_components(segment_adjacency(segments), directed=False)
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
                  bc_min=BC_MIN, lin_max=LIN_MAX, len_max=LEN_MAX, drop_orphans=True,
                  dens_max=DENS_MAX, dens_len_max=DENS_LEN_MAX, rescue_min_mm=RESCUE_MIN_MM):
    """Remove the sheets (sheet rule and dense rule), then the floating fragments,
    except the root-like ones, which are rescued.
    Returns (kept_segments, removed_segments, features); features['rescued'] lists
    the rescued segments, which are NOT in kept_segments (see reattach)."""
    f = segment_features(segments, skv_full, voxel_size)
    rule_removed = (f['bc'] >= bc_min) & (f['lin'] < lin_max) & (f['length'] < len_max)
    if dens_max is not None and dens_max > 0:
        rule_removed |= (f['dens'] >= dens_max) & (f['length'] < dens_len_max)
    kept = [s for s, r in zip(segments, rule_removed) if not r]
    removed = [s for s, r in zip(segments, rule_removed) if r]
    rescued = []
    if drop_orphans and base is not None and len(kept) > 1:
        kept, orphan = keep_base_component(kept, base)
        if orphan and rescue_min_mm and rescue_min_mm > 0:
            dens_of = {id(s): d for s, d in zip(segments, f['dens'])}
            _, lab = connected_components(segment_adjacency(orphan), directed=False)
            for c in np.unique(lab):
                comp = [o for o, l in zip(orphan, lab) if l == c]
                if (sum(o['length_mm'] for o in comp) >= rescue_min_mm and
                        np.median([dens_of[id(o)] for o in comp]) < (dens_max or np.inf)):
                    rescued += comp
            rid = set(id(o) for o in rescued)
            orphan = [o for o in orphan if id(o) not in rid]
        removed += orphan
    f['rule'] = rule_removed
    f['rescued'] = rescued
    return kept, removed, f


def reattach(segments, kept, rescued, voxel_size):
    """Join each rescued fragment to the kept system by the shortest path on the
    skeleton of the decomposition (which runs through the removed layer that hid the
    connection). Returns (bridge_voxels, joined, not_joined): the voxels of the
    bridging paths, and the rescued segments that could or could not be joined."""
    from scipy.sparse.csgraph import dijkstra
    from .graph_extraction import build_skel_graph
    if not rescued or not kept:
        return np.zeros((0, 3), int), [], list(rescued)
    vs = np.asarray(voxel_size, float)
    allv = np.unique(np.vstack([s['coords'] for s in segments]), axis=0)
    lo = allv.min(0); shp = tuple(allv.max(0) - lo + 1)
    sk = np.zeros(shp, bool); sk[tuple((allv - lo).T)] = True
    graph, vox, v2n = build_skel_graph(sk, voxel_size=vs)
    node = lambda C: v2n[tuple((np.asarray(C) - lo).T)]
    src = np.unique(np.concatenate([node(s['coords']) for s in kept]))
    dist, pred, _ = dijkstra(graph, directed=False, indices=src, min_only=True,
                             return_predecessors=True)
    _, lab = connected_components(segment_adjacency(rescued), directed=False)
    bridges, joined, not_joined = [], [], []
    for c in np.unique(lab):
        comp = [o for o, l in zip(rescued, lab) if l == c]
        nodes = np.unique(np.concatenate([node(o['coords']) for o in comp]))
        k = nodes[np.argmin(dist[nodes])]
        if not np.isfinite(dist[k]):
            not_joined += comp
            continue
        path = []
        while k >= 0:
            path.append(k)
            k = pred[k]
        bridges.append(vox[path] + lo)
        joined += comp
    B = np.vstack(bridges) if bridges else np.zeros((0, 3), int)
    return B, joined, not_joined
