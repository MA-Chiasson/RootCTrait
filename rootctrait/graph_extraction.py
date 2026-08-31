# graph_extraction.py
# Root Phenotyping Pipeline
# Author: [Mana Eskandari], [2026]

import numpy as np
import scipy.ndimage as ndimage
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path


def get_endpoints(sk):
    kernel = np.ones((3, 3, 3), dtype=int)
    neighbor_count = ndimage.convolve(sk.astype(int), kernel, mode='constant', cval=0)
    return (neighbor_count == 2) & sk


def get_branchpoints(sk):
    kernel = np.ones((3, 3, 3), dtype=int)
    neighbor_count = ndimage.convolve(sk.astype(int), kernel, mode='constant', cval=0)
    return (neighbor_count > 3) & sk


def build_skel_graph(sk, voxel_size=None):
    if voxel_size is None:
        voxel_size = [1.0, 1.0, 1.0]

    z, y, x   = np.nonzero(sk)
    num_nodes = len(z)

    vol2node          = np.full(sk.shape, -1, dtype=int)
    vol2node[z, y, x] = np.arange(num_nodes)
    voxels            = np.column_stack((z, y, x))

    dz, dy, dx = np.meshgrid([-1, 0, 1], [-1, 0, 1], [-1, 0, 1], indexing='ij')
    offsets     = np.column_stack((dz.ravel(), dy.ravel(), dx.ravel()))
    offsets     = offsets[np.any(offsets != 0, axis=1)]
    dist_offsets = np.sqrt(np.sum((offsets * voxel_size)**2, axis=1))

    edges_from, edges_to, weights = [], [], []

    for c, offset in enumerate(offsets):
        nz = z + offset[0]
        ny = y + offset[1]
        nx = x + offset[2]

        valid = ((nz >= 0) & (nz < sk.shape[0]) &
                 (ny >= 0) & (ny < sk.shape[1]) &
                 (nx >= 0) & (nx < sk.shape[2]))

        vz, vy, vx     = nz[valid], ny[valid], nx[valid]
        neighbor_nodes  = vol2node[vz, vy, vx]
        valid_neighbors = neighbor_nodes >= 0

        from_nodes = np.arange(num_nodes)[valid][valid_neighbors]
        to_nodes   = neighbor_nodes[valid_neighbors]

        edges_from.extend(from_nodes)
        edges_to.extend(to_nodes)
        weights.extend(np.full(len(from_nodes), dist_offsets[c]))

    graph = csr_matrix((weights, (edges_from, edges_to)), shape=(num_nodes, num_nodes))
    return graph, voxels, vol2node


def get_path_from_predecessors(pred, start_node, end_node):
    path = []
    curr = end_node
    while curr != start_node and curr >= 0:
        path.append(curr)
        curr = pred[curr]
    if curr == start_node:
        path.append(start_node)
        return path[::-1]
    return []


def prune_skeleton(sk, min_branch_length_vox=5):
    sk_pruned = sk.copy()

    while True:
        graph, voxels, vol2node = build_skel_graph(sk_pruned, voxel_size=[1.0, 1.0, 1.0])
        ep_mask = get_endpoints(sk_pruned)
        bp_mask = get_branchpoints(sk_pruned)

        ep_coords = np.argwhere(ep_mask)
        bp_coords = np.argwhere(bp_mask)

        if len(ep_coords) == 0 or len(bp_coords) == 0:
            break

        ep_nodes     = vol2node[ep_coords[:, 0], ep_coords[:, 1], ep_coords[:, 2]]
        bp_nodes     = vol2node[bp_coords[:, 0], bp_coords[:, 1], bp_coords[:, 2]]
        bp_nodes_set = set(bp_nodes)

        num_nodes    = graph.shape[0]
        v_edges_from = np.full(len(bp_nodes), num_nodes)
        v_edges_to   = bp_nodes
        v_weights    = np.zeros(len(bp_nodes))

        hop_weights = np.ones_like(graph.data)
        edges_from  = np.concatenate((graph.nonzero()[0], v_edges_from, v_edges_to))
        edges_to    = np.concatenate((graph.nonzero()[1], v_edges_to, v_edges_from))
        all_weights = np.concatenate((hop_weights, v_weights, v_weights))

        graph_bp    = csr_matrix((all_weights, (edges_from, edges_to)),
                                  shape=(num_nodes + 1, num_nodes + 1))
        D_bp, pred  = shortest_path(csgraph=graph_bp, directed=False,
                                     indices=num_nodes, return_predecessors=True)

        pruned_any     = False
        nodes_to_remove = []
        for epn in ep_nodes:
            if D_bp[epn] < min_branch_length_vox:
                curr       = pred[epn] if epn != num_nodes else -1
                path_nodes = [epn]
                while curr != num_nodes and curr >= 0 and curr not in bp_nodes_set:
                    path_nodes.append(curr)
                    curr = pred[curr]
                nodes_to_remove.extend(path_nodes)
                pruned_any = True

        if not pruned_any:
            break

        for n in set(nodes_to_remove):
            c = voxels[n]
            sk_pruned[c[0], c[1], c[2]] = False

    return sk_pruned
