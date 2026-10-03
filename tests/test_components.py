"""Unit tests for the connectivity, pivot return, root mask and reattachment helpers."""
import os
import sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from scipy.sparse.csgraph import connected_components
from rootctrait.decontamination import segment_adjacency, keep_base_component, reattach
from rootctrait.root_traits_full import pivot_return_mm, split_pivot_return, root_mask

VS = (0.39, 0.39, 0.2)


def _seg(i, pts, order=2):
    return {'seg_id': i, 'order': order, 'parent_seg': -1,
            'coords': np.asarray(pts, int), 'length_mm': float(len(pts))}


def test_corner_contact_connects():
    a = _seg(0, [[5, 5, 5], [6, 5, 5]], order=1)
    b = _seg(1, [[7, 6, 6], [8, 6, 6]])          # touches a by a voxel corner
    c = _seg(2, [[20, 20, 20], [21, 20, 20]])    # isolated
    n, lab = connected_components(segment_adjacency([a, b, c]), directed=False)
    assert n == 2 and lab[0] == lab[1] != lab[2]
    kept, orph = keep_base_component([a, b, c], base=(5, 5, 5))
    assert [s['seg_id'] for s in kept] == [0, 1] and [s['seg_id'] for s in orph] == [2]


def test_shared_voxel_connects():
    a = _seg(0, [[1, 1, 1], [2, 2, 2]], order=1)
    b = _seg(1, [[2, 2, 2], [9, 9, 9]])          # shares voxel (2, 2, 2)
    n, _ = connected_components(segment_adjacency([a, b]), directed=False)
    assert n == 1


def test_pivot_return_flag():
    down = [[d, 0, 0] for d in range(0, 40)]
    hook = down + [[39 - k, 1, 0] for k in range(1, 12)]   # 11 voxels up = 4.3 mm
    assert abs(pivot_return_mm(hook, VS) - 11 * VS[0]) < 1e-9
    d, u = split_pivot_return(hook, VS)
    assert d[-1][0] == 39 and len(u) == 12
    small = down + [[39 - k, 1, 0] for k in range(1, 5)]   # 4 voxels up = 1.6 mm
    d, u = split_pivot_return(small, VS)
    assert len(d) == len(small) and len(u) == 0
    assert pivot_return_mm(down, VS) == 0.0


def test_root_mask_excludes_removed_material():
    BW = np.zeros((30, 30, 30), bool)
    BW[2:28, 4:7, 4:7] = True        # root column
    BW[2:6, 15:25, 15:18] = True     # surface artifact
    root = _seg(0, [[d, 5, 5] for d in range(3, 27)], order=1)
    junk = _seg(1, [[3, x, 16] for x in range(16, 24)])
    M = root_mask(BW, [root, junk], [root], root['coords'], VS)
    assert M[2:28, 4:7, 4:7].all() and not M[2:6, 15:25, 15:18].any()


def test_root_mask_excludes_isolated_fragments():
    BW = np.zeros((30, 30, 30), bool)
    BW[2:28, 4:7, 4:7] = True        # root column
    BW[10:14, 20:24, 20:24] = True   # isolated fragment, not in the skeleton tree
    root = _seg(0, [[d, 5, 5] for d in range(3, 27)], order=1)
    M = root_mask(BW, [root], [root], root['coords'], VS)
    assert M[2:28, 4:7, 4:7].all() and not M[10:14, 20:24, 20:24].any()


def test_reattach_joins_rescued_fragment_through_removed_layer():
    kept = _seg(0, [[d, 5, 5] for d in range(0, 10)], order=1)
    layer = _seg(1, [[9, 5 + k, 5] for k in range(1, 6)])          # removed, links kept to the fragment
    frag = _seg(2, [[9 + k, 10, 5] for k in range(1, 8)])           # rescued fragment
    far = _seg(3, [[25, 25, 25], [26, 25, 25]])                     # cannot be joined
    B, joined, not_joined = reattach([kept, layer, frag, far], [kept], [frag, far], VS)
    assert [s['seg_id'] for s in joined] == [2] and [s['seg_id'] for s in not_joined] == [3]
    path = set(map(tuple, B))
    assert (9, 6, 5) in path and (9, 10, 5) in path or (10, 10, 5) in path
