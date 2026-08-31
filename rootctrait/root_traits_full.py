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


def compute_all_traits(segments, primary_path, base_rcp, BW, edt, voxel_size, skv_full):
    vs = np.asarray(voxel_size, float)
    base = np.asarray(base_rcp, float)
    baseY = base[0] * vs[0]
    basephys = base * vs
    T = {}

    # ---- primary ----
    PPv = np.asarray(primary_path, int)
    # "Hook" correction: a primary root plunges, it does not go back up. If the
    # pivot tip goes back up (skeletonization artifact, or the pivot "jumping" onto
    # a neighbouring upward root), cut the pivot at its deepest point. This avoids
    # inflating LRP with a non-biological upward portion.
    if len(PPv) > 2:
        deepest = int(np.argmax(PPv[:, 0]))
        if deepest < len(PPv) - 1:
            # only cut if the rise is significant (>3mm), otherwise tolerate the
            # small natural undulations of the pivot
            rise_mm = (PPv[deepest, 0] - PPv[-1, 0]) * vs[0]
            if rise_mm > 3.0:
                PPv = PPv[:deepest + 1]
    PP = PPv.astype(float) * vs
    LRP = _len_phys(PP)
    T['LRP'] = LRP

    # ---- per segment ----
    orders = np.array([s['order'] for s in segments])
    seglen = np.array([s['length_mm'] for s in segments], float)
    segdiam = np.array([_seg_diam_mm(s, edt) for s in segments], float)
    lat = orders >= 2

    T['TRL'] = float(seglen.sum())
    T['LTRL'] = float(seglen[lat].sum()) if lat.any() else 0.0
    T['MLRL'] = float(seglen[lat].mean()) if lat.any() else 0.0
    T['NRL'] = int(lat.sum())
    # distribution of lateral lengths into classes (mm)
    # Keys must match the pipeline output columns (NRL_court_<5, NRL_moyen_5_15,
    # NRL_long_>15). They were previously English (short/medium), which left the
    # court and moyen columns empty.
    ll = seglen[lat]
    T['NRL_court_<5'] = int((ll < 5).sum())
    T['NRL_moyen_5_15'] = int(((ll >= 5) & (ll < 15)).sum())
    T['NRL_long_>15'] = int((ll >= 15).sum())

    # ---- skeleton point cloud (depth, width, hull) ----
    A = skv_full.astype(float) * vs
    depth = A[:, 0] - baseY
    # "Top" correction: ignore material ABOVE the collar (negative depth) for the
    # depth traits. Those points (upward laterals, noise, or hypocotyl near the
    # collar) are not part of the root system below the collar and would bias PM
    # and the depth distribution. A small margin (-2mm) allows for collar position
    # uncertainty.
    below_collar = depth >= -2.0
    depth_valid = depth[below_collar]
    T['PM'] = float(depth_valid.max()) if len(depth_valid) else 0.0
    sd = np.sort(depth_valid)
    T['D50'] = float(np.percentile(sd, 50)) if len(sd) else 0.0
    T['D95'] = float(np.percentile(sd, 95)) if len(sd) else 0.0
    WX = float(A[:, 1].max() - A[:, 1].min()) if len(A) else 0.0
    WZ = float(A[:, 2].max() - A[:, 2].min()) if len(A) else 0.0
    T['WX'] = WX
    T['WZ'] = WZ
    T['LM'] = max(WX, WZ)
    T['RLP'] = (T['LM'] / T['PM']) if T['PM'] > 0 else np.nan

    def width_at(frac):
        if T['PM'] <= 0:
            return 0.0
        lo, hi = frac * T['PM'] - 2.0, frac * T['PM'] + 2.0
        m = (depth >= lo) & (depth <= hi)
        if m.sum() < 2:
            return 0.0
        return float(max(np.ptp(A[m, 1]), np.ptp(A[m, 2])))
    T['W25'] = width_at(0.25)
    T['W50'] = width_at(0.50)
    T['W75'] = width_at(0.75)

    # ---- angles ----
    angs = np.array([_angle_vs_vertical(s['coords'].astype(float) * vs) for s in segments])
    ok = ~np.isnan(angs)
    T['ANGsys'] = float(np.average(angs[ok], weights=seglen[ok])) if ok.any() else np.nan
    m_lat = lat & ok
    T['ACRL'] = float(np.mean(angs[m_lat])) if m_lat.any() else np.nan

    # ---- angles of order-2 roots only (system framework) ----
    # Targets the overall shape of the system (compact vs spreading) using only the
    # order-2 laterals (that branch off the pivot), excluding orders 3+ and the noise
    # mixed into ACRL. Three variants for comparison:
    #   ANGO2      : overall mean angle (start-to-end vector) vs vertical
    #   ANGO2_sd   : standard deviation of these angles (framework regularity)
    #   ANGO2_init : mean angle of the initial portion (4 mm) vs vertical
    ord2 = np.array([s['order'] == 2 for s in segments])
    m_o2 = ord2 & ok
    if m_o2.any():
        T['ANGO2'] = float(np.mean(angs[m_o2]))
        T['ANGO2_sd'] = float(np.std(angs[m_o2]))
        angs_init = np.array([_angle_initial_vs_vertical(s['coords'], vs)
                              for s in segments])
        m_o2i = ord2 & ~np.isnan(angs_init)
        T['ANGO2_init'] = float(np.mean(angs_init[m_o2i])) if m_o2i.any() else np.nan
    else:
        T['ANGO2'] = np.nan
        T['ANGO2_sd'] = np.nan
        T['ANGO2_init'] = np.nan


    # ---- convex hull ----
    try:
        T['CHV'] = float(ConvexHull(A).volume) if len(A) >= 4 else np.nan
    except Exception:
        T['CHV'] = np.nan

    # ---- volume and surface from the mask ----
    vvol = float(np.prod(vs))
    VRT = float(BW.sum()) * vvol
    T['VRT'] = VRT
    try:
        verts, faces, _, _ = measure.marching_cubes(BW.astype(np.uint8), level=0.5, spacing=tuple(vs))
        T['SRT'] = float(measure.mesh_surface_area(verts, faces))
    except Exception:
        T['SRT'] = np.nan
    if T['CHV'] and T['CHV'] > 0 and VRT > 0:
        _ic = VRT / T['CHV']
        T['IC'] = _ic if _ic <= 1.0 else np.nan      # IC > 1 = degenerate convex hull
    else:
        T['IC'] = np.nan
    T['SRL'] = (T['TRL'] / VRT) if VRT > 0 else np.nan  # mm/mm3

    # ---- topology on the segment tree (robust to skeleton pollution) ----
    seg_ids = [s['seg_id'] for s in segments]
    parents_lat = set(s['parent_seg'] for s in segments if s['order'] >= 2 and s['parent_seg'] >= 0)
    is_parent = set(s['parent_seg'] for s in segments if s['parent_seg'] >= 0)
    T['NT'] = int(sum(1 for sid in seg_ids if sid not in is_parent))   # leaf segments = apices
    T['NBP'] = int(len(parents_lat))                                    # branch sites carrying a lateral
    T['MaxO'] = int(orders.max()) if len(orders) else 0
    T['DR'] = (T['NRL'] / (LRP / 10.0)) if LRP > 0 else np.nan          # laterals per cm of pivot

    # ---- NTR: roots starting from the collar (proxy: proximal end <=4 mm from base, order <=2) ----
    nb = 0
    for s in segments:
        C = s['coords'].astype(float) * vs
        dmin = min(np.linalg.norm(C[0] - basephys), np.linalg.norm(C[-1] - basephys))
        if dmin <= 4.0 and s['order'] <= 2:
            nb += 1
    T['NTR'] = max(1, nb)

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
    T['DRP'] = float(np.mean(dia_pp)) if len(dia_pp) else np.nan
    T['DRS'] = float(np.mean(segdiam[lat])) if lat.any() else np.nan
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
