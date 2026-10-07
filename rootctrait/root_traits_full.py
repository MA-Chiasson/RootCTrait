"""Full computation of root traits on the ordered-segment decomposition.
Input: segments (dict order/coords/length_mm/parent_seg/seg_id), primary path,
base node, binary mask BW, physical distance map edt, voxel_size, and the full
list of skeleton voxels (for topology).
Output: dict of traits in native units (mm, mm2, mm3, degrees). Conversion to cm
is done when writing the table.
Axis convention: col0 = Y (depth), col1 = X, col2 = Z.
"""
import numpy as np
from scipy.spatial import ConvexHull, cKDTree
from skimage import measure


def _len_phys(P):
    return float(np.sum(np.linalg.norm(np.diff(P, axis=0), axis=1))) if len(P) >= 2 else 0.0


def _seg_diam_mm(seg, edt):
    c = seg['coords']
    return 2.0 * float(np.median(edt[c[:, 0], c[:, 1], c[:, 2]]))


def _angle_vs_vertical(C):
    """Angle (deg) of the segment global vector relative to the Y axis (vertical)."""
    if len(C) < 2:
        return np.nan
    d = C[-1] - C[0]
    n = np.linalg.norm(d)
    if n < 1e-9:
        return np.nan
    return float(np.degrees(np.arccos(np.clip(abs(d[0]) / n, -1.0, 1.0))))


def _angle_initial_vs_vertical(C, vs, mm=4.0):
    """Angle (deg) of the INITIAL PORTION of a segment (first `mm` mm after its
    insertion point) relative to the vertical (col0 = depth). Captures the starting
    direction of the root, more stable than the global vector because it is
    insensitive to curvature and to distal-tip noise."""
    P = C.astype(float) * vs
    if len(P) < 2:
        return np.nan
    cum = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))])
    idx = int(np.searchsorted(cum, mm))
    idx = max(1, min(idx, len(P) - 1))
    d = P[idx] - P[0]
    n = np.linalg.norm(d)
    if n < 1e-9:
        return np.nan
    return float(np.degrees(np.arccos(np.clip(abs(d[0]) / n, -1.0, 1.0))))


def _insertion_angle(seg, parent, vs, mm=4.0, win=6):
    """Angle (deg) between the lateral's initial direction (first `mm` mm after its
    insertion) and the local direction of its parent root at the attachment point.
    0 = the lateral leaves in the parent's direction, 90 = perpendicular. This is
    an emergence angle relative to the parent, independent of system orientation."""
    C = seg['coords'].astype(float) * vs
    if len(C) < 2:
        return np.nan
    cum = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(C, axis=0), axis=1))])
    idx = int(np.searchsorted(cum, mm)); idx = max(1, min(idx, len(C) - 1))
    d_lat = C[idx] - C[0]
    nl = np.linalg.norm(d_lat)
    Cp = parent['coords'].astype(float) * vs
    if nl < 1e-9 or len(Cp) < 2:
        return np.nan
    d_lat = d_lat / nl
    j = int(np.argmin(np.linalg.norm(Cp - C[0], axis=1)))     # attachment on parent
    lo = max(0, j - win); hi = min(len(Cp), j + win + 1)
    chunk = Cp[lo:hi]
    if len(chunk) < 2:
        return np.nan
    X = chunk - chunk.mean(0)
    d_par = np.linalg.svd(X, full_matrices=False)[2][0]
    d_par = d_par / (np.linalg.norm(d_par) + 1e-9)
    cos = abs(float(d_lat @ d_par))
    return float(np.degrees(np.arccos(np.clip(cos, 0.0, 1.0))))


HOOK_MM = 3.0   # upward return at the pivot tip above which it is flagged in the review figure (mm)


def pivot_return_mm(primary_path, voxel_size):
    """Upward return of the pivot after its deepest point, in mm (0 when the deepest
    point is the tip). A primary root that reaches the pot bottom or wall can turn
    back up, but the same shape also arises when the pivot follows a neighbouring
    upward root through a contact. The return is therefore not corrected: the whole
    path is kept for the traits, and the return is reported as a quality control
    value so that such samples can be checked in the review figure."""
    vs = np.asarray(voxel_size, float)
    P = np.asarray(primary_path, int)
    if len(P) < 2:
        return 0.0
    deepest = int(np.argmax(P[:, 0]))
    return float((P[deepest, 0] - P[-1, 0]) * vs[0])


def split_pivot_return(primary_path, voxel_size, flag_mm=HOOK_MM):
    """Split the pivot at its deepest point when its upward return exceeds `flag_mm`,
    for display only. Returns (descending_part, returning_part); returning_part is
    empty when the return is below the threshold. The traits always use the whole path."""
    P = np.asarray(primary_path, int)
    if len(P) > 2 and pivot_return_mm(P, voxel_size) > flag_mm:
        deepest = int(np.argmax(P[:, 0]))
        return P[:deepest + 1], P[deepest:]
    return P, P[:0]


def root_mask(BW, all_segments, root_segments, primary_path, voxel_size, chunk=500000):
    """Binary mask restricted to the cleaned root system. Every foreground voxel of
    BW is assigned to its nearest skeleton voxel (physical distance) among all the
    segments of the decomposition; it is kept when that skeleton voxel belongs to the
    cleaned roots (or to the primary path). Material attached to the hypocotyl, to
    removed sheets or to orphan fragments is therefore excluded, so that volume and
    surface describe the same root system as the skeleton traits. Mask components
    that are not connected to the cleaned roots (isolated fragments, which are not
    part of the skeleton tree) are excluded first, so that their voxels cannot be
    assigned to a distant root."""
    from scipy import ndimage
    vs = np.asarray(voxel_size, float)
    ref = [s['coords'] for s in all_segments]
    flag = [np.zeros(len(s['coords']), bool) for s in all_segments]
    keep_ids = set(id(s) for s in root_segments)
    for k, s in enumerate(all_segments):
        if id(s) in keep_ids:
            flag[k][:] = True
    ref.append(np.asarray(primary_path, int).reshape(-1, 3))
    flag.append(np.ones(len(ref[-1]), bool))
    ref = np.vstack(ref).astype(np.int64); flag = np.concatenate(flag)
    # a voxel shared by a kept and a removed segment counts as kept
    key = np.ravel_multi_index(ref.T, BW.shape)
    uk, inv = np.unique(key, return_inverse=True)
    uflag = np.zeros(len(uk), bool); np.logical_or.at(uflag, inv, flag)
    uvox = np.column_stack(np.unravel_index(uk, BW.shape))
    tree = cKDTree(uvox * vs)
    lab, _ = ndimage.label(BW, structure=np.ones((3, 3, 3), bool))
    kept_vox = uvox[uflag]
    keep_lab = np.unique(lab[kept_vox[:, 0], kept_vox[:, 1], kept_vox[:, 2]])
    keep_lab = keep_lab[keep_lab > 0]
    fg = np.argwhere(np.isin(lab, keep_lab))
    out = np.zeros_like(BW, dtype=bool)
    for i in range(0, len(fg), chunk):
        F = fg[i:i + chunk]
        _, j = tree.query(F * vs)
        F = F[uflag[j]]
        out[F[:, 0], F[:, 1], F[:, 2]] = True
    return out


def compute_all_traits(segments, primary_path, base_rcp, BW, edt, voxel_size, skv_full):
    vs = np.asarray(voxel_size, float)
    base = np.asarray(base_rcp, float)
    baseY = base[0] * vs[0]
    basephys = base * vs
    T = {}

    # ---- primary ----
    PPv = np.asarray(primary_path, int)
    PP = PPv.astype(float) * vs
    PRL = _len_phys(PP)
    T['PRL'] = PRL
    T['PIVOT_RETURN'] = pivot_return_mm(PPv, vs)   # quality control, not a trait

    # ---- per segment ----
    orders = np.array([s['order'] for s in segments])
    seglen = np.array([s['length_mm'] for s in segments], float)
    segdiam = np.array([_seg_diam_mm(s, edt) for s in segments], float)
    lat = orders >= 2

    T['TRL'] = float(seglen.sum())
    T['LTRL'] = float(seglen[lat].sum()) if lat.any() else 0.0
    T['MLRL'] = float(seglen[lat].mean()) if lat.any() else 0.0
    T['NLR'] = int(lat.sum())
    # distribution of lateral lengths into classes (mm)
    # Keys must match the pipeline output columns (COLS in run_pipeline.py).
    ll = seglen[lat]
    T['NLR_short_<5'] = int((ll < 5).sum())
    T['NLR_medium_5_15'] = int(((ll >= 5) & (ll < 15)).sum())
    T['NLR_long_>15'] = int((ll >= 15).sum())

    # ---- skeleton point cloud (depth, width, hull) ----
    A = skv_full.astype(float) * vs
    depth = A[:, 0] - baseY
    # "Top" correction: ignore material ABOVE the collar (negative depth) for the
    # depth traits. Those points (upward laterals, noise, or hypocotyl near the
    # collar) are not part of the root system below the collar and would bias MD
    # and the depth distribution. A small margin (-2mm) allows for collar position
    # uncertainty.
    below_collar = depth >= -2.0
    depth_valid = depth[below_collar]
    T['MD'] = float(depth_valid.max()) if len(depth_valid) else 0.0
    sd = np.sort(depth_valid)
    T['D50'] = float(np.percentile(sd, 50)) if len(sd) else 0.0
    T['D95'] = float(np.percentile(sd, 95)) if len(sd) else 0.0
    WX = float(A[:, 1].max() - A[:, 1].min()) if len(A) else 0.0
    WZ = float(A[:, 2].max() - A[:, 2].min()) if len(A) else 0.0
    T['WX'] = WX
    T['WZ'] = WZ
    T['MW'] = max(WX, WZ)
    T['WDR'] = (T['MW'] / T['MD']) if T['MD'] > 0 else np.nan

    def width_at(frac):
        if T['MD'] <= 0:
            return 0.0
        lo, hi = frac * T['MD'] - 2.0, frac * T['MD'] + 2.0
        m = (depth >= lo) & (depth <= hi)
        if m.sum() < 2:
            return 0.0
        return float(max(np.ptp(A[m, 1]), np.ptp(A[m, 2])))
    T['W25'] = width_at(0.25)
    T['W50'] = width_at(0.50)
    T['W75'] = width_at(0.75)

    # ---- angles of order-2 roots only (system framework) ----
    # ANGO2 (mean order-2 angle to vertical) and ANGO2_sd (its dispersion) describe the
    # inclination and regularity of the framework. Earlier variants to vertical
    # (ANGsys, ACRL, ANGO2_init) were dropped as redundant with ANGO2 (see the angle
    # correlation analysis); the insertion angle ANGI below adds an independent axis.
    angs = np.array([_angle_vs_vertical(s['coords'].astype(float) * vs) for s in segments])
    ok = ~np.isnan(angs)
    ord2 = np.array([s['order'] == 2 for s in segments])
    m_o2 = ord2 & ok
    if m_o2.any():
        T['ANGO2'] = float(np.mean(angs[m_o2]))
        T['ANGO2_sd'] = float(np.std(angs[m_o2]))
    else:
        T['ANGO2'] = np.nan
        T['ANGO2_sd'] = np.nan

    # ---- ANGI: mean insertion angle of order-2 laterals relative to their parent ----
    seg_by_id = {s['seg_id']: s for s in segments}
    ins = []
    for s in segments:
        if s['order'] >= 2 and s.get('parent_seg', -1) in seg_by_id:
            a = _insertion_angle(s, seg_by_id[s['parent_seg']], vs)
            if not np.isnan(a):
                ins.append(a)
    T['ANGI'] = float(np.mean(ins)) if ins else np.nan


    # ---- convex hull ----
    try:
        T['CHV'] = float(ConvexHull(A).volume) if len(A) >= 4 else np.nan
    except Exception:
        T['CHV'] = np.nan

    # ---- volume and surface from the mask (pass the mask restricted to the
    # cleaned roots, see root_mask, so that CI and RLV combine consistent quantities) ----
    vvol = float(np.prod(vs))
    TRV = float(BW.sum()) * vvol
    T['TRV'] = TRV
    try:
        verts, faces, _, _ = measure.marching_cubes(BW.astype(np.uint8), level=0.5, spacing=tuple(vs))
        T['TRSA'] = float(measure.mesh_surface_area(verts, faces))
    except Exception:
        T['TRSA'] = np.nan
    if T['CHV'] and T['CHV'] > 0 and TRV > 0:
        _ic = TRV / T['CHV']
        T['CI'] = _ic if _ic <= 1.0 else np.nan      # CI > 1 = degenerate convex hull
    else:
        T['CI'] = np.nan
    T['RLV'] = (T['TRL'] / TRV) if TRV > 0 else np.nan  # mm/mm3

    # ---- topology on the segment tree (robust to skeleton pollution) ----
    seg_ids = [s['seg_id'] for s in segments]
    parents_lat = set(s['parent_seg'] for s in segments if s['order'] >= 2 and s['parent_seg'] >= 0)
    is_parent = set(s['parent_seg'] for s in segments if s['parent_seg'] >= 0)
    T['NT'] = int(sum(1 for sid in seg_ids if sid not in is_parent))   # leaf segments = apices
    T['NBP'] = int(len(parents_lat))                                    # branch sites carrying a lateral
    T['MaxO'] = int(orders.max()) if len(orders) else 0
    T['LRD'] = (T['NLR'] / (PRL / 10.0)) if PRL > 0 else np.nan          # laterals per cm of pivot

    # ---- NCR: roots starting from the collar (proxy: proximal end <=4 mm from base, order <=2) ----
    nb = 0
    for s in segments:
        C = s['coords'].astype(float) * vs
        dmin = min(np.linalg.norm(C[0] - basephys), np.linalg.norm(C[-1] - basephys))
        if dmin <= 4.0 and s['order'] <= 2:
            nb += 1
    T['NCR'] = max(1, nb)

    # ---- IBD: spacing of lateral insertion points along the pivot ----
    if PP.shape[0] >= 2 and lat.any():
        arc = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(PP, axis=0), axis=1))])
        ptree = cKDTree(PP)
        pos = []
        for s in [s for s, fl in zip(segments, lat) if fl]:
            C = s['coords'].astype(float) * vs
            d0, i0 = ptree.query(C[0])
            d1, i1 = ptree.query(C[-1])
            idx = i0 if d0 <= d1 else i1
            pos.append(arc[idx])
        pos = np.sort(np.array(pos))
        T['IBD'] = float(np.mean(np.diff(pos))) if len(pos) >= 2 else np.nan
    else:
        T['IBD'] = np.nan

    # ---- diameters ----
    skv = skv_full.astype(int)
    dia_pp = 2.0 * edt[PPv[:, 0], PPv[:, 1], PPv[:, 2]]
    T['PRD'] = float(np.mean(dia_pp)) if len(dia_pp) else np.nan
    T['MLD'] = float(np.mean(segdiam[lat])) if lat.any() else np.nan
    T['DMAX'] = 2.0 * float(edt[skv[:, 0], skv[:, 1], skv[:, 2]].max()) if len(skv) else np.nan
    mu = np.mean(segdiam) if len(segdiam) else np.nan
    T['DD_cv'] = float(np.std(segdiam) / mu) if (mu and mu > 0) else np.nan
    # pivot taper: proximal 25% vs distal 25% diameter, relative, per cm
    if PP.shape[0] >= 4:
        arc = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(PP, axis=0), axis=1))])
        L = arc[-1]
        prox = dia_pp[arc <= 0.25 * L]
        dist = dia_pp[arc >= 0.75 * L]
        Dp = np.mean(prox) if len(prox) else np.nan
        Dd = np.mean(dist) if len(dist) else np.nan
        T['TAPER'] = float((Dp - Dd) / Dp / (L / 10.0)) if (Dp and Dp > 0 and L > 0) else np.nan
    else:
        T['TAPER'] = np.nan

    # ---- mean tortuosity of laterals (length / straight-line distance) ----
    tors = []
    for s in [s for s, fl in zip(segments, lat) if fl]:
        C = s['coords'].astype(float) * vs
        straight = np.linalg.norm(C[-1] - C[0])
        if straight > 1e-6:
            tors.append(s['length_mm'] / straight)
    T['TOR'] = float(np.mean(tors)) if tors else np.nan

    return T
