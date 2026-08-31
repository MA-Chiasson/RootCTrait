"""detection_hypocotyle.py

Bounded collar and hypocotyl detection, refined and validated version.

Principle (validated visually on both early and late time points):

1. VERTICAL COLUMN (vertical_column): the hypocotyl column is identified by a
   vertical seed above the collar, close to the axis. There is a single hypocotyl
   per plant, so only the connected column closest to the collar is kept.
   Propagation along the column requires verticality, which avoids absorbing the
   horizontal roots at collar level.

2. BOUNDED COLLAR (bounded_climb): the collar starts at the thickest point of the
   top (detect_base, in the pipeline), then CLIMBS along the continuous thick
   column (material near the axis, radius >= 0.7 x collar radius), STOPPING as soon
   as it reaches the hypocotyl. This places the start of the primary root at the
   true top of the fleshy base, without climbing into the stem. Returns base2 (the
   raised collar).

3. HIGH BRANCHES (add_high_branches): the branches hanging HIGH on the stem (their
   lowest point stays more than haut_min above the RAISED collar base2) are added
   to the hypocotyl, whatever their angle. Roots that come back down to the raised
   collar level stay roots. Height is measured relative to base2 (the raised
   collar), which correctly separates stem branches from basal roots.

High-level function: collar_and_hypocotyl(sk, edt, kept, base, voxel_size)
returns (base2, hypo_ids, gain_mm, column_pts).

The orphan cleanup that follows (removing fragments detached once the hypocotyl is
removed) is done in the pipeline by a single call to the existing keep_base_component
function (decontamination module).

Axis convention: col0 = depth (Y), col1 = X, col2 = Z. voxel (0.39, 0.39, 0.2).
"""
import numpy as np


def _mean_vertical_angle(seg, vs):
    """Mean angle (deg) of the segment with the vertical, length-weighted.
    0 = vertical, 90 = horizontal."""
    P = seg['coords'].astype(float) * vs
    d = np.diff(P, axis=0)
    sl = np.linalg.norm(d, axis=1)
    ok = sl > 1e-9
    if ok.sum() == 0:
        return 90.0
    cos = np.abs(d[ok, 0]) / sl[ok]
    return float(np.average(np.degrees(np.arccos(np.clip(cos, 0, 1))), weights=sl[ok]))


def vertical_column(segments, base, voxel_size,
                    dist_max_mm=10.0, marge_mm=2.0, angle_max=45.0):
    """Hypocotyl column: vertical seed + VERTICAL propagation above the collar.
    Returns (ids, seg_by_id, children_of)."""
    vs = np.asarray(voxel_size, float)
    base_y = base[0] * vs[0]
    base_xz = np.array([base[1] * vs[1], base[2] * vs[2]])

    def has_material_above(s):
        d = s['coords'].astype(float)[:, 0] * vs[0] - base_y
        return d.min() < -marge_mm and d.max() < marge_mm

    def is_seed(s):
        P = s['coords'].astype(float) * vs
        d = P[:, 0] - base_y
        if d.max() >= 0 or d.min() >= -marge_mm:
            return False
        if _mean_vertical_angle(s, vs) >= angle_max:
            return False
        dist_xz = np.min(np.sqrt((P[:, 1] - base_xz[0]) ** 2 + (P[:, 2] - base_xz[1]) ** 2))
        return dist_xz < dist_max_mm

    seg_by_id = {s['seg_id']: s for s in segments}
    children_of = {}
    for s in segments:
        children_of.setdefault(s['parent_seg'], []).append(s['seg_id'])

    seeds = [s['seg_id'] for s in segments if is_seed(s)]
    if not seeds:
        return set(), seg_by_id, children_of

    def chain_root(sid):
        ch = [sid]
        p = seg_by_id[sid]['parent_seg']
        while p in seg_by_id and p not in ch:
            ch.append(p)
            p = seg_by_id[p]['parent_seg']
        return ch[-1]

    def dist_to_collar(sid):
        P = seg_by_id[sid]['coords'].astype(float) * vs
        return np.min(np.sqrt((P[:, 1] - base_xz[0]) ** 2 + (P[:, 2] - base_xz[1]) ** 2))

    groups = {}
    for g in seeds:
        groups.setdefault(chain_root(g), []).append(g)
    best = min(groups.keys(), key=lambda r: min(dist_to_collar(g) for g in groups[r]))
    ids = set(groups[best])

    changed = True
    while changed:
        changed = False
        for s in segments:
            if s['seg_id'] in ids or not has_material_above(s):
                continue
            if _mean_vertical_angle(s, vs) >= angle_max:   # VERTICAL propagation
                continue
            if s['parent_seg'] in ids or any(e in ids for e in children_of.get(s['seg_id'], [])):
                ids.add(s['seg_id'])
                changed = True
    return ids, seg_by_id, children_of


def bounded_climb(sk, edt, base, hypo_vox, voxel_size=(0.39, 0.39, 0.2),
                  frac=0.7, dist_max=4.0, step_mm=1.0):
    """Climb the collar along the continuous thick column above the thickest point,
    stopping as soon as a depth slice is entirely hypocotyl. Returns
    (base2, gain_mm, column_pts).
    hypo_vox: set of (y, x, z) tuples of the skeleton voxels of the vertical column."""
    vs = np.asarray(voxel_size, float)
    sc = np.argwhere(sk)
    rad_all = edt[sc[:, 0], sc[:, 1], sc[:, 2]]
    depths = sc[:, 0] * vs[0]
    xz = sc[:, [1, 2]] * vs[[1, 2]]
    base_xz = base[[1, 2]] * vs[[1, 2]]
    dist = np.sqrt(((xz - base_xz) ** 2).sum(1))
    base_y = base[0] * vs[0]
    rad_collar = edt[base[0], base[1], base[2]]
    thr = frac * rad_collar
    is_hypo = np.array([tuple(v) in hypo_vox for v in sc])
    near = dist < dist_max
    lo = base_y
    pts = []
    while lo > depths.min() - step_mm:
        hi = lo
        lo = lo - step_mm
        m = near & (depths >= lo) & (depths < hi) & (rad_all >= thr)
        if not m.any():
            break
        if is_hypo[m].all():
            break
        m_ok = m & (~is_hypo)
        if not m_ok.any():
            break
        pts.append(sc[m_ok])
    if pts:
        allp = np.vstack(pts)
        dp = allp[:, 0] * vs[0]
        k = np.argmin(dp)
        return allp[k], float(base_y - dp[k]), allp
    return np.asarray(base), 0.0, np.empty((0, 3), int)


def add_high_branches(segments, column_ids, base, base2, voxel_size,
                      marge_mm=2.0, haut_min_mm=3.0):
    """Starting from the vertical column, ADD the attached branches hanging high
    above the RAISED collar base2 (their lowest point stays more than haut_min above
    base2), whatever their angle. Roots that come back down to base2 level are not
    added. Returns the final set of hypocotyl seg_ids."""
    vs = np.asarray(voxel_size, float)
    base_y = base[0] * vs[0]
    base2_y = base2[0] * vs[0]
    children_of = {}
    for s in segments:
        children_of.setdefault(s['parent_seg'], []).append(s['seg_id'])

    def above_base(s):
        d = s['coords'].astype(float)[:, 0] * vs[0] - base_y
        return d.min() < -marge_mm and d.max() < marge_mm

    def hangs_high(s):
        d = s['coords'].astype(float)[:, 0] * vs[0] - base2_y
        return d.max() < -haut_min_mm

    ids = set(column_ids)
    changed = True
    while changed:
        changed = False
        for s in segments:
            if s['seg_id'] in ids or not above_base(s) or not hangs_high(s):
                continue
            if s['parent_seg'] in ids or any(e in ids for e in children_of.get(s['seg_id'], [])):
                ids.add(s['seg_id'])
                changed = True
    return ids


def collar_and_hypocotyl(sk, edt, kept, base, voxel_size):
    """High-level function. Returns (base2, hypo_ids, gain_mm, column_pts).

    base2       : raised (bounded) collar, to use as the collar for the traits.
    hypo_ids    : seg_ids of the hypocotyl segments to exclude (column + high branches).
    gain_mm     : how far the collar was raised.
    column_pts  : voxels of the thick column (for pivot extension / figure).
    """
    col_ids, _, _ = vertical_column(kept, base, voxel_size)
    # The collar always climbs along the thick column (it goes up to the top of the
    # fleshy base). If a hypocotyl is present, it stops at its entrance; otherwise it
    # climbs to the end of the thick column.
    hypo_vox = set(tuple(v) for s in kept if s['seg_id'] in col_ids for v in s['coords'])
    base2, gain, column_pts = bounded_climb(sk, edt, base, hypo_vox, voxel_size)
    if col_ids:
        hypo_ids = add_high_branches(kept, col_ids, base, base2, voxel_size)
    else:
        hypo_ids = set()
    return base2, hypo_ids, gain, column_pts
