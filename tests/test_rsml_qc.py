"""Tests for the RSML export and the QC ranking."""
import os
import sys
import numpy as np
import xml.etree.ElementTree as ET

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from rootctrait.rsml_export import build_rsml
from tools.qc_rank import rank_batch, robust_z

VS = (0.39, 0.39, 0.2)


def _seg(i, pts, order, parent):
    P = np.asarray(pts, int)
    return {'seg_id': i, 'order': order, 'parent_seg': parent, 'coords': P,
            'length_mm': float(np.linalg.norm(np.diff(P * VS, axis=0), axis=1).sum())}


def _tree():
    prim = [[d, 10, 10] for d in range(0, 40)]
    lat = [[10, 10 + k, 10] for k in range(0, 12)]               # order 2 from the primary
    sub = [[10 + k, 16, 10] for k in range(0, 8)]                # order 3 from the lateral
    roots = [_seg(0, prim, 1, -1), _seg(1, lat, 2, 0), _seg(2, sub, 3, 1)]
    edt = np.full((50, 30, 20), 0.4)
    return roots, np.asarray(prim), edt


def test_rsml_nesting_and_geometry():
    roots, prim, edt = _tree()
    r = build_rsml('S', roots, prim, edt, VS, offset=(5, 0, 0)).getroot()
    assert r.find('metadata/unit').text == 'mm'
    p = r.find('scene/plant/root')
    assert p.get('label') == 'primary'
    l2 = p.findall('root'); assert len(l2) == 1 and l2[0].get('label') == 'order2'
    l3 = l2[0].findall('root'); assert len(l3) == 1 and l3[0].get('label') == 'order3'
    pts = p.findall('geometry/polyline/point')
    assert len(pts) == len(prim)
    # z is depth in mm, in the frame of the input volume (crop offset added)
    assert abs(float(pts[0].get('z')) - 5 * VS[0]) < 1e-6
    samples = p.findall('functions/function/sample')
    assert len(samples) == len(prim) and abs(float(samples[0].get('value')) - 0.8) < 1e-6


def test_rsml_lateral_with_dropped_parent_hangs_on_nearest_root():
    roots, prim, edt = _tree()
    roots[2]['parent_seg'] = 99                                    # parent not kept
    r = build_rsml('S', roots, prim, edt, VS).getroot()
    l2 = r.find('scene/plant/root').findall('root')
    assert [e.get('label') for e in l2[0].findall('root')] == ['order3']


def test_robust_z_flags_outlier():
    z = robust_z([10, 11, 9, 10, 12, 10, 40])
    assert z[-1] > 3 and np.all(np.abs(z[:-1]) < 3)


def test_rank_puts_flagged_samples_first():
    rows = [{'name': f'S{i}', 'T': {'LRP': 60 + i, 'PM': 50, 'TRL': 300, 'NRL': 40, 'ANGO2': 45, 'DRP': 1.5,
                                    '%removed': 20, 'N_RESCUED': 1, 'COLLAR_RAISE': 2, 'HYPOCOTYL_LEN': 0,
                                    'PIVOT_RETURN': 0, 'LRP_PM': (60 + i) / 50}} for i in range(8)]
    rows[3]['T']['PIVOT_RETURN'] = 5.0
    rows[5]['T']['%removed'] = 90
    out = rank_batch(rows)
    assert {out[0]['ID'], out[1]['ID']} == {'S3', 'S5'}
    assert all(d['n_flags'] == 0 for d in out[2:])
