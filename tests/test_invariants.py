"""Invariant tests for the RootCTrait pipeline.

Runs the full pipeline on the bundled synthetic example (example/roots/sample_S1.npy)
and checks the internal-consistency invariants that must hold for any sample. Run
with pytest (`pytest -q`) or directly (`python tests/test_invariants.py`).
"""
import os
import sys
import numpy as np
import scipy.ndimage as ndimage
from skimage.morphology import skeletonize

# make the package importable when run from the repo root or from tests/
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from rootctrait.io_volume import load_volume
from rootctrait.graph_extraction import prune_skeleton
from rootctrait.root_decomposition import decompose_root_system
from rootctrait.decontamination import decontaminate, keep_base_component
from rootctrait.root_traits_full import compute_all_traits
from rootctrait.detection_hypocotyle import collar_and_hypocotyl

VOXEL = (0.39, 0.39, 0.2)
EXAMPLE = os.path.join(ROOT, "example", "roots", "sample_S1.npy")


def detect_base(sk, edt):
    sc = np.argwhere(sk)
    thr = np.percentile(sc[:, 0], 20)
    top = sc[sc[:, 0] <= thr]
    rad = edt[top[:, 0], top[:, 1], top[:, 2]]
    return top[np.argmax(rad)]


def _traits_on_example():
    V = load_volume(EXAMPLE)
    BW = V > 0.5 * np.max(V)
    c = np.argwhere(BW); mn = c.min(0); mx = c.max(0); m = 6
    BW = BW[tuple(slice(max(0, mn[i] - m), mx[i] + m) for i in range(3))]
    edt = ndimage.distance_transform_edt(BW, sampling=VOXEL)
    sk = prune_skeleton(skeletonize(BW).astype(bool), 5)
    base = detect_base(sk, edt)
    segs, prim, vox, _, _ = decompose_root_system(sk, base, list(VOXEL),
                                                  dist_map=edt, min_seg_len_mm=2.0)
    kept, removed, _ = decontaminate(segs, vox, VOXEL, base=base, drop_orphans=True)
    base2, hypo_ids, gain, _ = collar_and_hypocotyl(sk, edt, kept, base, VOXEL)
    roots = [s for s in kept if s['seg_id'] not in hypo_ids]
    roots, _ = keep_base_component(roots, base)
    skv = np.vstack([s['coords'] for s in roots]) if roots else vox
    return compute_all_traits(roots, np.vstack([base2, prim]), base2, BW, edt, VOXEL, skv)


T = _traits_on_example()


def test_example_produces_traits():
    assert T is not None and 'TRL' in T and T['TRL'] > 0


def test_trl_ge_lrp():
    assert T['TRL'] >= T['LRP'] - 1e-6


def test_ltrl_le_trl():
    assert T['LTRL'] <= T['TRL'] + 1e-6


def test_compactness_le_one():
    assert T['IC'] <= 1.0 + 1e-6


def test_nrl_equals_sum_of_classes():
    s = T['NRL_court_<5'] + T['NRL_moyen_5_15'] + T['NRL_long_>15']
    assert s == T['NRL']


def test_angles_in_range():
    for k in ('ANGO2', 'ANGO2_sd', 'ANGO2_init', 'ANGsys'):
        v = T.get(k)
        if v is not None and not (isinstance(v, float) and np.isnan(v)):
            assert 0.0 <= v <= 90.0, f"{k}={v} out of [0,90]"


def test_core_traits_not_nan():
    for k in ('LRP', 'TRL', 'NRL', 'PM'):
        v = T[k]
        assert not (isinstance(v, float) and np.isnan(v)), f"{k} is NaN"


if __name__ == "__main__":
    tests = [f for n, f in sorted(globals().items()) if n.startswith('test_')]
    ok = 0
    for f in tests:
        try:
            f(); print(f"PASS  {f.__name__}"); ok += 1
        except AssertionError as e:
            print(f"FAIL  {f.__name__}: {e}")
    print(f"\n{ok}/{len(tests)} tests passed")
    sys.exit(0 if ok == len(tests) else 1)
