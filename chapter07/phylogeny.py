'''
NOTE
====

Trees are represented with a Pyhton dictionary of adjacency lists. Each key-value in the dictionary has one node for the
key, and a list of pairs (adj_node, adj_weight) for the value; adj_node is a node adjacent to the one in the key, and
adj_weight is the weight of the connecting edge. The tree may or may not have a root node. Note that, if node 'a' is
adjacent to node 'b', then 'b' is adjacent to 'a'.

Matrices are represented with a Python dictionary of dictionaries; given a matrix 'm', m[r][c] indicates its element
of row 'r' and column 'c'.

'''

from collections import deque


def add_node(tree, node1, node2, weight):
    """
    Add a weighted edge to a tree. Any end-point of the edge (node) that is not already in the tree, is added to the
    tree too.
    :param tree: The tree, a dictionary of dictionaries with its adjacency lists.
    :param node1: The first edge end-point, a tree node.
    :param node2: The second edge endpoint, a tree node.
    :param weight: The weight of the edge connecting the two nodes.
    """

    def add_monodir_edge(tree, node1, node2, weight):
        adjs = tree.get(node1, [])
        adjs.append((node2, weight))
        tree[node1] = adjs

    add_monodir_edge(tree, node1, node2, weight)
    add_monodir_edge(tree, node2, node1, weight)


def dist_between_leaves(tree):
    """
    Returns a matrix with the distance between every pair of leaf nodes in a weighted tree.
    :param tree: The tree, a dictionary of dictionaries with its adjacency lists.
    :return: The matrix, as a dictionary of dictionaries.
    """

    def set_distance(distance, node1, node2, value):
        """
        Set the value in a matrix, at given row and column coordinates, keeping the matrix symmetric; setting element
        [row][column] also sets the element [column][row] to the same value.
        :param distance: The matrix, a dictionary of dictionaries.
        :param node1: The row coordinate.
        :param node2: The column coordinate.
        :param value: The value to be set.
        """

        def set_monodir_distance(distance, node1, node2, value):
            if distance.get(node1) is None:
                distance[node1] = {}
            distance[node1][node2] = value

        set_monodir_distance(distance, node1, node2, value)
        if node1 != node2:
            set_monodir_distance(distance, node2, node1, value)

    ''' The matrix of distances between each pair of nodes is a dict, where distance[a][b] is the distance between node 
    a and node b; the matrix is square, symmetric, with all 0s on the diagonal. The matrix will grow one row and one
    column at a time, as the tree is traversed. '''
    distance = {}

    # Select any of the tree nodes, from where to start the tree traversal
    starting_node = next(iter(tree))
    # Distance of the starting node to itself is 0
    set_distance(distance, starting_node, starting_node, 0)

    # Keep track of nodes already traversed, to ensure they are not traversed again
    traversed_nodes = set()
    # Use a queue to schedule nodes for traversal, breadth-first
    pending_nodes = deque([starting_node])
    while pending_nodes:
        current_node = pending_nodes.popleft()
        traversed_nodes.add(current_node)
        ''' For every node adjacent to the current one, which hasn't been traversed yet, determine its distance
        from all nodes already in the matrix, and add it to the matrix (add rows/columns as necessary).'''
        for adj_node, adj_weight in tree[current_node]:
            if adj_node in traversed_nodes:
                continue
            set_distance(distance, adj_node, adj_node, 0)
            # Set of nodes for which distances to/from current_node have been computed already
            computed_nodes = set(distance) - {adj_node}
            for other_node in computed_nodes:
                set_distance(distance,
                             adj_node,
                             other_node,
                             adj_weight if other_node == current_node else distance[current_node][
                                                                               other_node] + adj_weight)
            # Queue the node for traversal
            pending_nodes.append(adj_node)

    # Remove from the distances matrix all rows/columns related to internal (non-leaf) nodes
    for node1 in set(distance):
        if len(tree[node1]) > 1:
            del distance[node1]
        else:
            for node2 in set(distance[node1]):
                if len(tree[node2]) > 1:
                    del distance[node1][node2]

    return distance


def limb_length(leaf, distance_matrix):
    """
    Returns the limb length of a given leaf in a tree with given distance matrix.
    :param leaf: The given leaf, a non-negative integer.
    :param distance_matrix: The distance matrix, a dictionary of dictionaries; it must be additive.
    :return: The limb length, an integer.
    """
    n = len(distance_matrix)
    length = min(
        [(distance_matrix[i][leaf] + distance_matrix[leaf][k] - distance_matrix[i][k]) // 2 for i in range(0, n) for k
         in range(0, n) if i != leaf != k])
    return length
