'''
NOTE
====

Trees are represented with a Pyhton dictionary of adjacency lists. Each key-value in the dictionary has one node for the
key, and a pair (parent-info, children-info) for the value; parent-info is a pair (parent, weight), which gives the
parent of the node and teh weight of the connecting edge; children-info is a list of pairs (child, weight), giving
all the children of the node along with the weights of the respective edges. If there is a root node in the tree,
its parent-info point to itself, with weight 0.
'''

from collections import deque


def add_node(tree, node1, node2, weight):
    def add_monodir_edge(tree, node1, node2, weight):
        adjs = tree.get(node1, [])
        adjs.append((node2, weight))
        tree[node1] = adjs

    add_monodir_edge(tree, node1, node2, weight)
    add_monodir_edge(tree, node2, node1, weight)


def dist_between_leaves(tree):
    def set_distance(distance, node1, node2, value):
        def set_monodir_distance(distance, node1, node2, value):
            if distance.get(node1) is None:
                distance[node1] = {}
            distance[node1][node2] = value

        set_monodir_distance(distance, node1, node2, value)
        if node1 != node2:
            set_monodir_distance(distance, node2, node1, value)

    # The matrix of distances between each pair of nodes is a dict, where distance[a][b] is the distance between node a
    # and node b; the matrix is square, symmetric, with all 0s on the diagonal.
    distance = {}

    # Select any of the tree nodes to start with
    starting_node = next(iter(tree))
    # Distance to itself is 0
    set_distance(distance, starting_node, starting_node, 0)

    processed_nodes = set()
    pending_nodes = deque([starting_node])
    while pending_nodes:
        current_node = pending_nodes.popleft()
        processed_nodes.add(current_node)
        for adj_node, adj_weight in tree[current_node]:
            if adj_node in processed_nodes:
                continue
            set_distance(distance, adj_node, adj_node, 0)
            # Set of nodes for which distances to/from current_node have been computed already
            computed_nodes = set(distance) - set([adj_node])
            for other_node in computed_nodes:
                set_distance(distance,
                             adj_node,
                             other_node,
                             adj_weight if other_node == current_node else distance[current_node][
                                                                               other_node] + adj_weight)
            pending_nodes.append(adj_node)

    for node1 in set(distance):
        if len(tree[node1]) > 1:
            del distance[node1]
        else:
            for node2 in set(distance[node1]):
                if len(tree[node2]) > 1:
                    del distance[node1][node2]

    return distance
