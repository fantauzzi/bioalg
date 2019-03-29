'''
NOTE
====

Trees and graphs are represented with a Pyhton dictionary of adjacency lists. Each key-value in the dictionary has one
vertex for the key, and a list of pairs (adj_vertex, adj_weight) for the value; adj_vertex is a vertex adjacent to the
one in the key, and adj_weight is the weight of the connecting edge. When it is a tree, it may or may not have a root
node. Note that, if node 'a' is adjacent to node 'b', then 'b' is adjacent to 'a', graphs and trees are not directed.

Matrices are represented with a Python dictionary of dictionaries; given a matrix 'm', m[r][c] indicates its element
of row 'r' and column 'c'.

'''

from itertools import product
from collections import deque, namedtuple

BinTreeAdj = namedtuple('BinTreeAdj', ['left', 'right', 'parent'])  # Adjacency list for a tree node


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


def hamming_distance(s1, s2):
    assert len(s1) == len(s2)
    dist = sum([char1 != char2 for char1, char2 in zip(s1, s2)])
    return dist


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
        [(distance_matrix[i][leaf] + distance_matrix[leaf][k] - distance_matrix[i][k]) / 2 for i in range(0, n) for
         k
         in range(0, n) if i != leaf != k])
    assert length == int(length)
    return int(length)


def path_in_a_tree(tree, v1, v2):
    """
    Returns the path in a weighted tree from one vertex to another, along with the overall weights from the starting vertex to each vertex in the path.
    :param tree: The adjacency lists for the tree, a dictionary of lists.
    :param v1: The vertex where the path shall start.
    :param v2: The vertex where the path shall end.
    :return: The path from v1 to v2, including v1 and v2, a list of pairs; in each pair, the first element is the traversed vertex, and the second element is the sum of the edge weights from v1 to that vertex, a number.
    """
    # Do a depth-first exploration of the tree starting from v1, until you reach v2

    pending = [(v1, 0)]  # Stack of nodes to be processed, each with its distance from v1
    unvisited = set(tree)
    current = None
    precedent = {}  # Track from which vertex any vertex is reached, and the cumulative weight to get there from v1
    while current != v2 and unvisited:
        current, current_dist = pending.pop()
        unvisited -= {current}
        for next_v, next_dist in tree[current]:
            if next_v in unvisited:
                pending.append((next_v, current_dist + next_dist))
                precedent[next_v] = current, current_dist
    assert current == v2
    path = [(current, current_dist)]
    # Backtrack from v2 to v1 using precedent
    while current != v1:
        current, dist = precedent[current]
        path.append((current, dist))
    path.reverse()
    return path


def break_edge_with_node(tree, node1, node2, new_node, dist_from_node1, dist_from_node2):
    """
    Insert a node along an edge in a weighted tree; the edge is removed, and the inserted node is connected to the previous edge endpoints.
    :param tree: The weighted tree adjacency lists, a dictionary of lists.
    :param node1: The first node delimiting the edge.
    :param node2: The second node delimiting the edge.
    :param new_node: The node to be inserted.
    :param dist_from_node1: The weight of the new edge between the first node and the inserted node.
    :param dist_from_node2: The weight of the new edge between the second node and the inserted node.
    """

    def remove_adjacent(tree, node1, node2):
        new_adj = []
        for adj_node, adj_dist in tree[node1]:
            if adj_node != node2:
                new_adj.append((adj_node, adj_dist))
        tree[node1] = new_adj

    remove_adjacent(tree, node1, node2)
    remove_adjacent(tree, node2, node1)
    add_node(tree, node1, new_node, dist_from_node1)
    add_node(tree, node2, new_node, dist_from_node2)


def additive_phylogeny(d):
    """
    Returns a simple tree fitting a given distance matrix, solving the distance based phylogeny problem.
    :param d: The distance matrix, must be additive; a dictionary of dictionaries.
    :return:The adjacency lists of the tree fitting the matrix, a dictionary of lists.
    """

    def remove_row_and_col(d, to_be_removed):
        """
        Remove a given row and column from a distance matrix.
        :param d: The distance matrix, a dictionary of dictionaries.
        :param to_be_removed: The key of the row and column to be removed. Both the row and the column with the given
        key will be removed.
        """
        del d[to_be_removed]
        for row in d:
            del d[row][to_be_removed]

    ''' Proceed by recursion. Remove one row and one column from the distance matrix, fit a tree to the smaller matrix,
    then add to the tree the node corresponding to the row and column previously removed. '''

    n = len(d)
    if n == 2:  # Recursion base case, a 2x2 distance matrix, the fitting simple tree is trivial
        tree = {0: [(1, d[0][1])],
                1: [(0, d[1][0])]}
        return tree

    # Determine the length of the limb for the node to be inserted in the tree
    limb = limb_length(n - 1, d)
    # Compute the "balded" distance matrix, based on the limb length
    for j in range(0, n - 1):
        d[j][n - 1] -= limb
        d[n - 1][j] = d[j][n - 1]

    ''' Find two leaves i and k such that they satisfy the equality of the limb length theorem. A new vertex will be
    added along the path from i to k '''
    for row, col in product(range(0, n), repeat=2):
        if d[row][col] == d[row][n - 1] + d[n - 1][col]:
            i, k = row, col
            break
    else:
        assert False

    # Determine the distance x from vertex i to the node where the limb for the new node should be connected.
    x = d[i][n - 1]
    # Remove the last row and last column from the distance matrix
    remove_row_and_col(d, n - 1)
    # Fit a simple tree on the smaller (trimmed) matrix
    tree = additive_phylogeny(d)
    # Determine a (potentially new) node in the path from i to k at distance x from i where to connect the limb for n-1
    selected_node = None
    path = path_in_a_tree(tree, i, k)
    for i_node in range(0, len(path)):
        node, dist_from_i = path[i_node]
        if dist_from_i == x:
            selected_node = node
            assert len(tree[selected_node]) > 1
            break
        if i_node < len(path) - 1:
            next_node, next_dist_from_i = path[i_node + 1]
            if next_dist_from_i > x:
                selected_node = min([-1, min(tree) - 1])
                break_edge_with_node(tree=tree,
                                     node1=node,
                                     node2=next_node,
                                     new_node=selected_node,
                                     dist_from_node1=x - dist_from_i,
                                     dist_from_node2=next_dist_from_i - x)
                break
    assert selected_node is not None
    ''' Add the new node to the fitting tree, along with its limb (a tree edge), and, if not already in the tree,
    the node where to attach the limb. 
    '''
    add_node(tree, selected_node, n - 1, limb)

    return tree


def upgma(d):
    """
    Returns the ultrametric tree produced by UPGMA (Unweighted Pair Group Method with Arithmetic Mean) for a given distance matrix.
    :param d: The distance matrix; a dictionary of dictionaries.
    :return: The adjacency lists for the tree, a dictionary of lists.
    """
    n = len(d)
    # For each cluster, keep track of its cardinality. Start with one cluster per leaf node, each cluster with cardinality one.
    clusters = {vertex: 1 for vertex in range(0, n)}
    # Start with a (non-connected) graph with one vertex per leaf (cluster), and no edges.
    graph = {vertex: [] for vertex in clusters}
    # Initialise the age for every vertex in the graph
    age = {vertex: 0 for vertex in graph}
    # Continue until you have clustered everything in one cluster only (i.e. you have built the root of the tree)
    while len(clusters) > 1:
        current_min = float('inf')  # Find i, j with i!=k that mininimises d[i][k].
        min_idxs = None
        for i in d:
            for j in d:
                if i != j and d[i][j] < current_min:
                    current_min = d[i][j]
                    min_idxs = i, j
        assert min_idxs is not None
        i, j = min_idxs
        assert i != j
        # Get the new cluster number. They get numbered n, n+1, n+2, ...
        new_cluster = max(clusters) + 1
        # Set the age for the new cluster
        age[new_cluster] = d[i][j] / 2
        # Add it as a vertex to the graph, connected to vertices (clusters) i and j
        add_node(graph, i, new_cluster, age[new_cluster] - age[i])
        add_node(graph, j, new_cluster, age[new_cluster] - age[j])
        # Update the matrix d, by removing rows and columns for i and j, and inserting a new row and column for new_cluster
        new_row = {col: (d[i][col] * clusters[i] + d[j][col] * clusters[j]) / (clusters[i] + clusters[j]) for col in
                   d
                   if col != i and col != j}
        new_row[new_cluster] = 0
        del d[i]
        del d[j]
        for row in d:
            del d[row][i]
            del d[row][j]
        d[new_cluster] = new_row
        for row in d:
            d[row][new_cluster] = new_row[row]
        # Add the new cluster to 'clusters' along with its cardinality
        clusters[new_cluster] = clusters[i] + clusters[j]
        # Clusters that have been clusterd together are moved from 'clusters'
        del clusters[i]
        del clusters[j]

    return graph


def total_distance(d):
    """
    Returns the total distance vector for a given distance matrix.
    :param d: The distance matrix, a dictionary of dictionaries.
    :return: The total distance, a dictionary, where item [i] is the sum of items on the i-th row of the matrix.
    """
    dist = {i: sum([d[i][k] for k in d if k != i]) for i in d}
    return dist


def neighbor_joining_matrix(d):
    """
    Returns the neighbor-joining matrix obtained from a given distance matrix.
    :param d: The distance matrix, a dictionary of dictionaries.
    :return: The neighbor-joining matrix, a dictionary of dictionaries.
    """
    n = len(d)
    total_dist = total_distance(d)
    d_star = {i: {j: (n - 2) * d[i][j] - total_dist[i] - total_dist[j] if j != i else 0 for j in d} for i in d}
    return d_star


def neighbor_joining(d):
    """
    Returns the tree fitting a given distance matrix via neighbor joining. The obtained tree has no root.
    :param d: The distance matrix; a dictionary of dictionaries.
    :return: The adjacency lists for the tree, a dictionary of lists.
    """
    n = len(d)
    # For each cluster, keep track of its cardinality. Start with one cluster per leaf node, each cluster with cardinality one.
    clusters = {vertex: 1 for vertex in range(0, n)}
    # Start with a (non-connected) graph with one vertex per leaf (cluster), and no edges.
    graph = {vertex: [] for vertex in clusters}
    # Continue until you have clustered everything in two clusters only (obtaining two trees)
    while len(clusters) > 2:
        # Get the neighbor joining matrix and its total distance vector.
        d_star = neighbor_joining_matrix(d)
        tot_dist = total_distance(d)
        # Find i, j with i!=k that minimises d[i][k].
        current_min = float('inf')
        min_idxs = None
        for i in d_star:
            for j in d_star:
                if i != j and d_star[i][j] < current_min:
                    current_min = d_star[i][j]
                    min_idxs = i, j
        assert min_idxs is not None
        i, j = min_idxs
        assert i != j
        # Get the new cluster number. They get numbered n, n+1, n+2, ...
        new_cluster = max(clusters) + 1
        # Add it as a vertex to the graph, connected to vertices (clusters) i and j
        delta_ij = (tot_dist[i] - tot_dist[j]) / (n - 2)
        limb_i = (d[i][j] + delta_ij) / 2
        limb_j = (d[i][j] - delta_ij) / 2
        add_node(graph, i, new_cluster, limb_i)
        add_node(graph, j, new_cluster, limb_j)
        # Update the matrix d, by removing rows and columns for i and j, and inserting a new row and column for new_cluster
        new_row = {col: (d[col][i] + d[col][j] - d[i][j]) / 2 for col in d if col != i and col != j}
        new_row[new_cluster] = 0
        del d[i]
        del d[j]
        for row in d:
            del d[row][i]
            del d[row][j]
        d[new_cluster] = new_row
        for row in d:
            d[row][new_cluster] = new_row[row]
        # Update n, as matrix d is now smaller by 1 row and 1 column
        n = len(d)
        # Add the new cluster to 'clusters' along with its cardinality
        clusters[new_cluster] = clusters[i] + clusters[j]
        # Clusters that have been clusterd together are moved from 'clusters'
        del clusters[i]
        del clusters[j]

    # Connect the two obtained clusters with an edge with appropriate weight. There is no root in the resulting tree.
    assert len(d) == 2
    v1, v2 = list(d.keys())
    add_node(graph, v1, v2, weight=d[v1][v2])
    return graph


def symbols_hamming_dist(symbol1, symbol2):
    return 0 if symbol1 == symbol2 else 1


def make_small_parsiomony_tree(strings, alphabet):
    n_leaves = len(strings)
    assert n_leaves % 2 == 0
    n_nodes = 2 * n_leaves - 1

    # Determine the set of symbols used in 'strings'.
    # alphabet = tuple(sorted(set.union(*[set(item) for item in strings])))
    # Compute the score for the leaves of the tree
    s = {node: {symbol: 0 if strings[node] == symbol else float('inf') for symbol in alphabet} for node in
         range(0, n_leaves)}

    tree = {}  # Adjacency lists for the tree, it associates the node number with a BinTreeAdj.
    ''' Nodes in the tree are numbered in decreasing order from the root going down to the leaves. The root is numbered
    n_nodes - 1 and the leaves from n_leaves-1 down to 0. '''
    next_for_insertion = n_nodes - 1  # Tracks the numbering of the next node to be inserted in the tree.
    ''' First node to be processed and inserted in the tree is the root. Pairs of the queue contain the number for the
    node and its parent number; parent of the root node is set to None.'''
    to_be_inserted = deque([(next_for_insertion, None)])
    next_for_insertion -= 1
    while to_be_inserted:
        current_node, parent = to_be_inserted.popleft()
        if current_node < n_leaves:  # If the node is a leaf...
            left, right = None, None
        else:  # If it is an internal node...
            left = next_for_insertion - 1
            right = next_for_insertion
            next_for_insertion -= 2
            to_be_inserted.extend([(right, current_node), (left, current_node)])
        tree[current_node] = BinTreeAdj(left=left, right=right, parent=parent)

    ''' Now that the tree has been built, queue its leaves for processing. They are queue in ascending order to
    facilitate debugging. '''
    to_be_processed = deque(sorted({tree[leaf].parent for leaf in range(0, n_leaves)}))
    while to_be_processed:
        current_node = to_be_processed.popleft()
        # Compute the score for 'current_node' and store it in 's'.
        daughter = tree[current_node].left
        son = tree[current_node].right
        s[current_node] = {k: min([s[daughter][i] + symbols_hamming_dist(i, k) for i in alphabet])
                              + min([s[son][j] + symbols_hamming_dist(j, k) for j in alphabet])
                           for k in alphabet}
        # If score has been computed already for 'current_node' and its sibling, then queue their parent for processing.
        parent = tree[current_node].parent
        if parent is None:  # The root node has no parent
            continue
        son_of_parent = tree[parent].left
        sibling = son_of_parent if son_of_parent != current_node else tree[parent].right
        if s.get(sibling) is not None:
            to_be_processed.append(parent)
    return tree, s


def small_parsimony(strings):
    def lowest_score_symbol(scores):
        min_score = float('inf')
        min_symbol = None
        for symbol, score in scores.items():
            if score < min_score:
                min_score = score
                min_symbol = symbol
        assert min_symbol is not None
        return min_symbol, min_score

    def append_to_best_string(symbol, best_string, node):
        the_list = best_string.get(node, '')
        the_list = the_list + symbol
        best_string[node] = the_list

    n_leaves = len(strings)
    assert n_leaves % 2 == 0
    n_nodes = 2 * n_leaves - 1
    root = n_nodes - 1
    # A tuple with the alphabet used by strings, in lexicographic order
    alphabet = tuple(sorted(set.union(*[set(item) for item in strings])))
    master_tree = None
    total_parsimony = 0
    best_string = {}

    for pos in range(0, len(strings[0])):
        characters = tuple(item[pos] for item in strings)
        tree, s = make_small_parsiomony_tree(characters, alphabet)
        if master_tree is None:
            master_tree = tree
        assert tree == master_tree

        root_symbol, parsimony = lowest_score_symbol(s[root])
        total_parsimony += parsimony
        to_be_processed = deque([(root, root_symbol)])
        append_to_best_string(root_symbol, best_string, root)

        while to_be_processed:
            current_node, current_node_symbol = to_be_processed.popleft()
            for child in (tree[current_node].left, tree[current_node].right):
                if child is None:
                    continue
                best_symbol = None
                lowest_score = float('inf')
                for symbol in alphabet:
                    score = s[child][symbol] + symbols_hamming_dist(current_node_symbol, symbol)
                    if score < lowest_score:
                        lowest_score = score
                        best_symbol = symbol
                assert best_symbol is not None
                append_to_best_string(best_symbol, best_string, child)
                to_be_processed.append((child, best_symbol))

    return tree, best_string, total_parsimony
