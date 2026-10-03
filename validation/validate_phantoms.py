"""validate_phantoms.py

Validation of RootCTrait against synthetic phantoms of KNOWN geometry.

There is no ground truth for real 3D roots, so we build synthetic root systems
whose total length, number of laterals and branching angles are known exactly,
run the full pipeline on them, and compare measured vs true. This quantifies the
recovery accuracy independently of any segmentation step.

Run: python validation/validate_phantoms.py
Writes validation/phantom_results.csv
"""
import os, sys, csv
import numpy as np
import scipy.ndimage as ndimage
from skimage.morphology import skeletonize

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from rootctrait.graph_extraction import prune_skeleton
from rootctrait.root_decomposition import decompose_root_system
from rootctrait.decontamination import decontaminate, keep_base_component
from rootctrait.root_traits_full import compute_all_traits
from rootctrait.detection_hypocotyle import collar_and_hypocotyl

VOXEL = np.array([0.39, 0.39, 0.2])
SHAPE = (200, 140, 140)


def _capsule(V, p0, p1, r):
    p0 = np.array(p0, float); p1 = np.array(p1, float)
    lo = np.maximum(np.floor(np.minimum(p0, p1) - r - 1).astype(int), 0)
    hi = np.minimum(np.ceil(np.maximum(p0, p1) + r + 1).astype(int), np.array(V.shape) - 1)
    ys, xs, zs = np.mgrid[lo[0]:hi[0]+1, lo[1]:hi[1]+1, lo[2]:hi[2]+1]
    P = np.stack([ys, xs, zs], -1).astype(float)
    d = p1 - p0; L2 = (d*d).sum()
    t = np.clip(((P - p0) * d).sum(-1) / (L2 + 1e-9), 0, 1)
    proj = p0 + t[..., None] * d
    V[lo[0]:hi[0]+1, lo[1]:hi[1]+1, lo[2]:hi[2]+1] |= np.sqrt(((P - proj)**2).sum(-1)) <= r


def _phys_len(p0, p1):
    return float(np.linalg.norm((np.array(p1) - np.array(p0)) * VOXEL))


def _angle_vert(p0, p1):
    d = (np.array(p1) - np.array(p0)) * VOXEL
    return float(np.degrees(np.arccos(np.clip(abs(d[0]) / (np.linalg.norm(d) + 1e-9), 0, 1))))


def make_phantom(n_lat, lat_len_vox, spread_vox, pivot_len_vox=150, seed=0):
    """Straight pivot + n_lat straight order-2 laterals. Returns (volume, truths)."""
    rng = np.random.default_rng(seed)
    V = np.zeros(SHAPE, bool)
    cx = cz = SHAPE[1] // 2
    top = 10
    piv0 = (top, cx, cz); piv1 = (top + pivot_len_vox, cx, cz)
    _capsule(V, piv0, piv1, 2.4)
    LRP_true = _phys_len(piv0, piv1)
    TRL_true = LRP_true
    angles = []
    ys = np.linspace(top + 20, top + pivot_len_vox - 20, n_lat).astype(int)
    for i, y in enumerate(ys):
        ang = rng.uniform(-1, 1)
        ex = cx + int(spread_vox * (1 if i % 2 else -1))
        ez = cz + int(spread_vox * 0.4 * ang)
        ey = y + lat_len_vox
        p0 = (y, cx, cz); p1 = (ey, ex, ez)
        _capsule(V, p0, p1, 1.2)
        TRL_true += _phys_len(p0, p1)
        angles.append(_angle_vert(p0, p1))
    truths = {'LRP': LRP_true, 'TRL': TRL_true, 'NRL': n_lat,
              'ANGO2': float(np.mean(angles))}
    return V.astype(np.uint8), truths


def measure(V):
    """Run the full pipeline (run_pipeline.analyse_mask) with the default settings."""
    import run_pipeline as rp
    ctx = dict(voxel_size=tuple(VOXEL), prune_vox=5, min_seg_len_mm=2.0, bc_min=3, lin_max=0.7,
               len_max=15, drop_orphans=True, dens_max=35, dens_len_max=6.0, rescue_min_mm=5.0)
    return rp.analyse_mask(V.astype(bool), ctx)['T']


def main():
    phantoms = [
        ("P1_4lat_30", dict(n_lat=4, lat_len_vox=35, spread_vox=22, seed=1)),
        ("P2_6lat_45", dict(n_lat=6, lat_len_vox=40, spread_vox=35, seed=2)),
        ("P3_8lat_45", dict(n_lat=8, lat_len_vox=38, spread_vox=34, seed=3)),
        ("P4_6lat_60", dict(n_lat=6, lat_len_vox=45, spread_vox=52, seed=4)),
        ("P5_10lat", dict(n_lat=10, lat_len_vox=36, spread_vox=32, seed=5)),
        ("P6_5lat_long", dict(n_lat=5, lat_len_vox=55, spread_vox=30, pivot_len_vox=170, seed=6)),
    ]
    rows = []
    print(f"{'phantom':>13} {'trait':>6} {'true':>8} {'measured':>9} {'err %':>7}")
    for name, kw in phantoms:
        V, truth = make_phantom(**kw)
        M = measure(V)
        for tr in ['LRP', 'TRL', 'NRL', 'ANGO2']:
            t, m = truth[tr], M[tr]
            err = 100 * (m - t) / t if t else float('nan')
            rows.append([name, tr, round(t, 2), round(m, 2), round(err, 1)])
            print(f"{name:>13} {tr:>6} {t:8.2f} {m:9.2f} {err:7.1f}")
        print()
    with open(os.path.join(ROOT, "validation", "phantom_results.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["phantom", "trait", "true", "measured", "err_pct"]); w.writerows(rows)
    # aggregate absolute error per trait
    print("=== mean absolute error per trait ===")
    for tr in ['LRP', 'TRL', 'NRL', 'ANGO2']:
        errs = [abs(r[4]) for r in rows if r[1] == tr]
        print(f"  {tr:>6}: {np.mean(errs):.1f}% mean abs error")


if __name__ == "__main__":
    main()
