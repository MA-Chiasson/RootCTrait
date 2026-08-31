# root_decomposition.py
# Decomposition of the root system into individual segments with branching order.
#
# Principle:
#   - the skeleton tree is rooted at the collar (base_rcp);
#   - at each bifurcation, the branch that plunges deepest (Y max) keeps the
#     current order (main axis), the others take order + 1;
#   - the deepest axis from the collar is therefore order 1 (pivot), what grafts
#     onto it is order 2 (secondary), then 3 (tertiary), etc.
#
# Returns the list of all segments (one per "root"), each with its voxel
# coordinates, its order, and the id of its parent segment.

import numpy as np
from collections import deque
from scipy.sparse.csgraph import shortest_path
from .graph_extraction import build_skel_graph


def _select_primary_tip(voxels, rad, base_node, D, pred, deg, vs,
                        weights, crown_mm):
    """
    Choose the pivot tip via a multi-criteria score among the candidate axes
    (collar -> descending tip paths).

    Criteria, z-normalized then combined (weights = weights):
      proximal diameter below the collar (dominant), verticality, depth,
      straightness, and deep connections (branch points below the surface).
    Returns the node index of the chosen tip, or None.
    """
    b = voxels[base_node].astype(float)
    n = len(voxels)
    eps = [i for i in range(n)
           if i != base_node and np.isfinite(D[i]) and deg[i] == 1
           and (voxels[i, 0] - b[0]) * vs[0] >= 8.0]
    if not eps:
        eps = [i for i in range(n)
               if i != base_node and np.isfinite(D[i]) and deg[i] == 1]
    if not eps:
        return None

    def path_nodes(nd):
        out = []
        cu = nd
        while cu >= 0 and cu != base_node:
            out.append(cu)
            cu = pred[cu]
        out.append(base_node)
        return out[::-1]

    depth, diam, vert, straight, conn = [], [], [], [], []
    for e in eps:
        pth = path_nodes(e)
        P = voxels[pth].astype(float) * vs
        tip = voxels[e].astype(float)
        dep = (tip[0] - b[0]) * vs[0]
        cum = np.concatenate([[0.0],
                              np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))])
        r = rad[pth]
        prox = r[(cum >= crown_mm) & (cum <= crown_mm + 12.0)]
        if len(prox) == 0:
            prox = r[cum >= crown_mm] if np.any(cum >= crown_mm) else r
        vec = (tip - b) * vs
        nv = np.linalg.norm(vec)
        v_ = np.degrees(np.arccos(min(1.0, abs(vec[0]) / nv))) if nv > 0 else 90.0
        nd_depth = (voxels[pth, 0] - b[0]) * vs[0]
        depth.append(dep)
        diam.append(2.0 * float(np.mean(prox)))
        vert.append(v_)
        straight.append(nv / cum[-1] if cum[-1] > 0 else 0.0)
        conn.append(int(np.sum((deg[pth] >= 3) & (nd_depth > 15.0))))

    def z(x):
        x = np.asarray(x, float)
        s = x.std()
        return (x - x.mean()) / s if s > 1e-9 else x * 0.0

    wd, wv, wp, ws, wc = weights
    score = (wd * z(diam) + wv * z([-v for v in vert]) + wp * z(depth)
             + ws * z(straight) + wc * z(conn))
    return eps[int(np.argmax(score))]


def decompose_root_system(sk, base_rcp, voxel_size, dist_map=None,
                          min_seg_len_mm=2.0, crown_exclude_mm=0.0,
                          score_weights=(2.0, 1.5, 1.5, 0.5, 0.5),
                          crown_score_mm=4.0):
    graph, voxels, vol2node = build_skel_graph(sk, voxel_size=voxel_size)
    n = graph.shape[0]
    base_node = int(vol2node[base_rcp[0], base_rcp[1], base_rcp[2]])
    vs = np.array(voxel_size)

    # Spanning tree rooted at the collar (geodesic distances + predecessors)
    D, pred = shortest_path(graph, directed=False, indices=base_node,
                            return_predecessors=True)
    reachable = np.isfinite(D)
    reach_nodes = np.where(reachable)[0]

    # List of children of each node in the tree
    children = [[] for _ in range(n)]
    for node in range(n):
        if node == base_node or not reachable[node]:
            continue
        p = pred[node]
        if p >= 0:
            children[p].append(node)

    # Maximum depth reached in each node's subtree (used to order the laterals:
    # the child with the deepest subtree continues).
    sub_max = voxels[:, 0].astype(float).copy()
    for nd in reach_nodes[np.argsort(-D[reach_nodes])]:
        for c in children[nd]:
            if sub_max[c] > sub_max[nd]:
                sub_max[nd] = sub_max[c]

    # --- Pivot choice -----------------------------------------------------
    # Default (without a distance map): old rule = deepest tip. With dist_map:
    # multi-criteria score (diameter, verticality, depth, straightness, deep
    # connections).
    deg = np.asarray((graph != 0).sum(1)).ravel()
    primary_tip = None
    if dist_map is not None:
        rad = dist_map[voxels[:, 0], voxels[:, 1], voxels[:, 2]].astype(float)
        primary_tip = _select_primary_tip(voxels, rad, base_node, D, pred, deg,
                                           vs, score_weights, crown_score_mm)
    if primary_tip is None:
        primary_tip = int(reach_nodes[np.argmax(voxels[reach_nodes, 0])])

    # Pivot path (order 1): from the collar to the chosen tip
    path_nodes = []
    cur = primary_tip
    while cur != base_node and cur >= 0:
        path_nodes.append(cur)
        cur = pred[cur]
    path_nodes.append(base_node)
    path_nodes = path_nodes[::-1]
    primary_path = voxels[path_nodes]
    on_primary = set(path_nodes)

    # Order assignment: order 1 follows exactly the chosen pivot; at each branch
    # on the pivot, the pivot's child stays order 1 and the others start at
    # order 2; within a lateral, the child with the deepest subtree continues
    # the order, the others take +1.
    branch_order = np.full(n, -1, dtype=int)
    branch_order[base_node] = 1
    q = deque([base_node])
    while q:
        nd = q.popleft()
        o = branch_order[nd]
        ch = children[nd]
        if not ch:
            continue
        if nd in on_primary:
            for c in ch:
                branch_order[c] = o if c in on_primary else o + 1
                q.append(c)
        else:
            ch_sorted = sorted(ch, key=lambda c: sub_max[c], reverse=True)
            for i, c in enumerate(ch_sorted):
                branch_order[c] = o if i == 0 else o + 1
                q.append(c)

    # Building the segments: a segment = connected chain of nodes of the same
    # order. A node starts a new segment if its parent has a different order
    # (or if it is the collar).
    seg_id_of = np.full(n, -1, dtype=int)
    segments  = []  # each segment: dict(order, nodes, parent_seg)
    # We traverse from the collar to the leaves so the parent already exists
    bfs = deque([base_node])
    visited = np.zeros(n, dtype=bool)
    visited[base_node] = True
    # root of the first segment (order 1 starts from the collar)
    while bfs:
        nd = bfs.popleft()
        for c in children[nd]:
            if visited[c]:
                continue
            visited[c] = True
            same_order = (branch_order[c] == branch_order[nd])
            if same_order and seg_id_of[nd] >= 0:
                sid = seg_id_of[nd]
                segments[sid]['nodes'].append(c)
            else:
                parent_seg = seg_id_of[nd]
                sid = len(segments)
                segments.append({'order': int(branch_order[c]),
                                 'nodes': [nd, c],
                                 'parent_seg': int(parent_seg)})
            seg_id_of[c] = sid
            bfs.append(c)

    # Conversion to coordinates + length, length filter + collar zone filter
    base_xyz = voxels[base_node] * vs   # collar in mm
    out = []
    for sid, seg in enumerate(segments):
        coords = voxels[seg['nodes']]
        P = coords * vs
        length_mm = float(np.sum(np.sqrt(np.sum(np.diff(P, axis=0) ** 2, axis=1)))) \
            if len(coords) >= 2 else 0.0
        if length_mm < min_seg_len_mm:
            continue
        # excludes segments whose closest point to the collar is within the
        # collar zone (unreliable skeleton node at the collar/hypocotyl)
        if crown_exclude_mm > 0.0:
            dmin = float(np.min(np.linalg.norm(P - base_xyz, axis=1)))
            if dmin < crown_exclude_mm:
                continue
        out.append({'seg_id': sid,
                    'order': seg['order'],
                    'parent_seg': seg['parent_seg'],
                    'coords': coords,
                    'length_mm': length_mm})

    # primary_path was already computed above from the tip chosen by the score
    # (or the deepest tip as a fallback).
    return out, primary_path, voxels, branch_order, base_node


def flag_parallel_pollution(segments, voxel_size, rpar=6.0, min_parallel=3,
                            cos_par=0.9, cos_offset=0.5):
    """
    Flag the segments belonging to a parallel sheet (artifact signature).

    A segment is flagged if it has at least `min_parallel` neighbors (centroids
    within `rpar` mm) that are both parallel (|cos angle| > cos_par) and laterally
    offset (|cos of the offset| < cos_offset), i.e. side by side and not in the
    continuation. True roots, which diverge, have few; a sheet of parallel
    strands has many.

    Fixed geometric rule, identical for all samples (compatible with a GWAS
    study). Returns a boolean mask the size of `segments`.
    """
    import numpy as np
    from scipy.spatial import cKDTree
    n = len(segments)
    if n == 0:
        return np.zeros(0, dtype=bool)
    vs = np.array(voxel_size)
    cents = np.array([(s['coords'] * vs).mean(0) for s in segments])
    dirs = np.zeros((n, 3))
    for i, s in enumerate(segments):
        C = s['coords'] * vs
        X = C - C.mean(0)
        d = np.linalg.svd(X, full_matrices=False)[2][0] if len(C) >= 2 \
            else np.array([1.0, 0.0, 0.0])
        dirs[i] = d / (np.linalg.norm(d) + 1e-9)
    tree = cKDTree(cents)
    flag = np.zeros(n, dtype=bool)
    nbrs = tree.query_ball_point(cents, rpar)
    for i in range(n):
        nb = np.array([j for j in nbrs[i] if j != i])
        if len(nb) < min_parallel:
            continue
        cs = np.abs(dirs[nb] @ dirs[i])
        o = cents[nb] - cents[i]
        no = np.linalg.norm(o, axis=1)
        no[no < 1e-6] = 1e9
        oc = np.abs((o / no[:, None]) @ dirs[i])
        if int(np.sum((cs > cos_par) & (oc < cos_offset))) >= min_parallel:
            flag[i] = True
    return flag


def order_summary(segments):
    """Count the number of roots and the total length per order."""
    summary = {}
    for s in segments:
        o = s['order']
        if o not in summary:
            summary[o] = {'count': 0, 'total_length_mm': 0.0}
        summary[o]['count'] += 1
        summary[o]['total_length_mm'] += s['length_mm']
    return dict(sorted(summary.items()))


def save_interactive_html(paths, orders, base_rcp, voxel_size, out_path):
    """
    Save an interactive 3D figure (rotation, zoom, hover) in HTML format.
    One trace per order, clickable in the legend to show/hide it.
    Returns True if the figure was created, False if plotly is not installed.
    """
    try:
        import plotly.graph_objects as go
    except ImportError:
        print("[VIZ] plotly not installed: interactive figure skipped "
              "(pip install plotly)")
        return False

    from collections import defaultdict
    vs = voxel_size
    palette = {1: 'red', 2: '#1f77b4', 3: '#2ca02c', 4: '#ff7f0e', 5: '#9467bd'}

    by_order = defaultdict(lambda: {'x': [], 'y': [], 'z': []})
    for p, o in zip(paths, orders):
        d = by_order[o]
        d['x'] += list(p[:, 1] * vs[1]) + [None]   # X
        d['y'] += list(p[:, 2] * vs[2]) + [None]   # Z
        d['z'] += list(p[:, 0] * vs[0]) + [None]   # Y (depth)

    fig = go.Figure()
    for o in sorted(by_order):
        d = by_order[o]
        name = 'Pivot (order 1)' if o == 1 else f'order {o}'
        fig.add_trace(go.Scatter3d(
            x=d['x'], y=d['y'], z=d['z'], mode='lines',
            line=dict(color=palette.get(o, 'gray'), width=7 if o == 1 else 3),
            name=name))

    fig.add_trace(go.Scatter3d(
        x=[base_rcp[1] * vs[1]], y=[base_rcp[2] * vs[2]], z=[base_rcp[0] * vs[0]],
        mode='markers', marker=dict(color='green', size=6), name='Collar'))

    fig.update_layout(
        title=f'{len(paths)} roots  -  free rotation with the mouse',
        scene=dict(xaxis_title='X (mm)', yaxis_title='Z (mm)',
                   zaxis=dict(title='Y depth (mm)', autorange='reversed')))
    fig.write_html(out_path, include_plotlyjs=True)
    return True
