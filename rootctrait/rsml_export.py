"""Export of the cleaned root tree in the Root System Markup Language (RSML).

RSML (Lobet et al., 2015, Plant Physiology 167, 617-627) is the community format
for root architectures. One file is written per sample, with one plant whose
primary root is the primary path from the raised collar to its tip, and whose
laterals are nested in their parent root, as in the ordered tree of RootCTrait.

Coordinates are in millimetres in the frame of the input volume:
    x = second axis (X), y = third axis (Z), z = first axis (depth, downward).
Each root carries its polyline (from its insertion to its tip), a diameter
function sampled at every point (2 edt, mm) and its branching order.
The hypocotyl, the removed artifacts and the detached orphans are not exported.
"""
import datetime
import numpy as np
import xml.etree.ElementTree as ET


def _fmt(v):
    return f'{float(v):.3f}'


def _add_root(parent_el, rid, label, order, pts_vox, off, edt, vs):
    r = ET.SubElement(parent_el, 'root', ID=str(rid), label=label)
    props = ET.SubElement(r, 'properties')
    ET.SubElement(props, 'order', value=str(int(order)))
    poly = ET.SubElement(ET.SubElement(r, 'geometry'), 'polyline')
    P = (np.asarray(pts_vox, int) + off) * vs
    for p in P:
        ET.SubElement(poly, 'point', x=_fmt(p[1]), y=_fmt(p[2]), z=_fmt(p[0]))
    fun = ET.SubElement(ET.SubElement(r, 'functions'), 'function', name='diameter', domain='polyline')
    q = np.asarray(pts_vox, int)
    for d in 2.0 * edt[q[:, 0], q[:, 1], q[:, 2]]:
        ET.SubElement(fun, 'sample', value=_fmt(d))
    return r


def build_rsml(name, roots, primary_path, edt, voxel_size, offset=(0, 0, 0), software='RootCTrait'):
    """Return an ElementTree of the RSML document.

    roots         cleaned root segments (dicts with seg_id, order, parent_seg, coords),
                  coords in the frame of `edt` (the cropped volume)
    primary_path  primary path from the raised collar to the tip, same frame
    offset        position of the cropped volume in the input volume (voxels)
    """
    vs = np.asarray(voxel_size, float)
    off = np.asarray(offset, int)
    rsml = ET.Element('rsml')
    meta = ET.SubElement(rsml, 'metadata')
    ET.SubElement(meta, 'version').text = '1'
    ET.SubElement(meta, 'unit').text = 'mm'
    ET.SubElement(meta, 'resolution').text = '1'
    ET.SubElement(meta, 'last-modified').text = datetime.date.today().isoformat()
    ET.SubElement(meta, 'software').text = software
    ET.SubElement(meta, 'file-key').text = str(name)
    pdefs = ET.SubElement(meta, 'property-definitions')
    pd = ET.SubElement(pdefs, 'property-definition')
    ET.SubElement(pd, 'label').text = 'order'; ET.SubElement(pd, 'type').text = 'integer'
    fdefs = ET.SubElement(meta, 'function-definitions')
    fd = ET.SubElement(fdefs, 'function-definition')
    ET.SubElement(fd, 'label').text = 'diameter'; ET.SubElement(fd, 'type').text = 'float'
    ET.SubElement(fd, 'unit').text = 'mm'
    plant = ET.SubElement(ET.SubElement(rsml, 'scene'), 'plant', ID='1', label=str(name))

    # The primary root is the primary path (raised collar to tip). Order 1 segments
    # are parts of it, so their children hang on the primary root.
    prim_el = _add_root(plant, 'primary', 'primary', 1, primary_path, off, edt, vs)
    by_id = {s['seg_id']: s for s in roots}
    order1 = {s['seg_id'] for s in roots if s['order'] == 1}
    lat = [s for s in roots if s['order'] >= 2]

    # parent of every lateral: its parent segment when it was kept, otherwise the
    # nearest kept segment of lower order (the parent may have been dropped by the
    # length filter or the cleaning steps), otherwise the primary root.
    def parent_of(s):
        p = s.get('parent_seg', -1)
        if p in by_id and by_id[p]['order'] < s['order']:
            return p
        if p in by_id and by_id[p]['order'] == s['order']:   # same order chain cut by a filter
            return parent_of(by_id[p])
        x0 = s['coords'][0] * vs
        best, bd = None, np.inf
        for t in roots:
            if t['order'] < s['order']:
                d = np.min(np.linalg.norm(t['coords'] * vs - x0, axis=1))
                if d < bd:
                    best, bd = t['seg_id'], d
        if best is None:
            return None
        dprim = np.min(np.linalg.norm(np.asarray(primary_path) * vs - x0, axis=1))
        return None if dprim <= bd else best

    children = {}
    for s in lat:
        p = parent_of(s)
        key = 'primary' if (p is None or p in order1) else p
        children.setdefault(key, []).append(s)

    def add_children(el, key):
        for s in sorted(children.get(key, []), key=lambda t: t['seg_id']):
            c = _add_root(el, f"r{s['seg_id']}", f"order{s['order']}", s['order'], s['coords'], off, edt, vs)
            add_children(c, s['seg_id'])

    add_children(prim_el, 'primary')
    ET.indent(rsml, space='  ')
    return ET.ElementTree(rsml)


def write_rsml(path, name, roots, primary_path, edt, voxel_size, offset=(0, 0, 0), software='RootCTrait'):
    tree = build_rsml(name, roots, primary_path, edt, voxel_size, offset, software)
    tree.write(path, encoding='utf-8', xml_declaration=True)
    return path
