import numpy as np

""" 
NOTE
====

The code uses adjacency lists to represent weighted directed acyclic graphs (DAGs). Adjacency list are a Python dictionary. Keys in the dictionary are vertices, and values are lists of pairs. Every pair is given by a vertex adjacent to the vertex of the key, and the weight of the conencting edge. Vertices are named with strings. Weights are non-negative integer numbers. If a vertex has no outgoing edges, but only incoming edges, then it may or may not appear among the keys.
"""


def dp_change(money, coins):
    """
    Returns the smallest number of coins that allow to change a given amount of money using available coin denominations.
    :param money: The given amount of money, an integer number.
    :param coins: The available coin denominations, a sequence of integer numbers.
    :return: The smalles number of coins to change the money, an integer.
    """
    min_num_coins = [0]
    for m in range(1, money + 1):
        min_num_coins.append(float('inf'))
        for i in range(0, len(coins)):
            if m >= coins[i]:
                if min_num_coins[m - coins[i]] + 1 < min_num_coins[m]:
                    min_num_coins[m] = min_num_coins[m - coins[i]] + 1
    return min_num_coins[money]


def manhattan_tourist(down, right):
    """
    Returns the longest path for the Manhattan tourist problem.
    :param down: A matrix of the distances to travel downward in the map, a numpy array of integer numbers.
    :param right: A matrix of the distances to travel right in the map, a numpy array of integer numbers.
    :return: The length of the longest path, an integer.
    """
    n, _ = down.shape
    _, m = right.shape
    s = np.zeros((n + 1, m + 1), dtype=np.int)
    for i in range(1, n + 1):
        s[i, 0] = s[i - 1, 0] + down[i - 1, 0]
    for j in range(1, m + 1):
        s[0, j] = s[0, j - 1] + right[0, j - 1]

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s[i, j] = max(s[i - 1, j] + down[i - 1, j], s[i, j - 1] + right[i][j - 1])
    return s[-1, -1]


def topological_ordering(adj):
    """
    Returns the vertices of a weighted DAG in topological order.
    :param adj: The adjacency lists of the weighted DAG.
    :return: The topological order, a list of vertices (string).
    """

    ''' Determine the number of incoming edges for each vertex in the graph. Relevant down below, will need to know which vertices don't have any incoming edge. '''
    fan_in = {}  # Associates each vertex with its count of incoming vertices
    for vertex, adjs in adj.items():
        if fan_in.get(vertex) is None:
            fan_in[vertex] = 0
        for adj_vertex, _ in adjs:
            current_count = fan_in.get(adj_vertex, 0)
            fan_in[adj_vertex] = current_count + 1

    topo_order = []
    pending = [vertex for vertex, count in fan_in.items() if count == 0]  # All vertices with no incoming edges
    while pending:
        current = pending.pop()
        topo_order.append(current)
        adjs = adj.get(current)
        if adjs is None:
            continue
        for vertex, _ in adjs:
            vertex_fan_in = fan_in[vertex]
            vertex_fan_in -= 1
            fan_in[vertex] = vertex_fan_in
            if vertex_fan_in == 0:
                pending.append(vertex)

    # If the graph is a DAG, then all edges have been traversed. Verify that it is a DAG.
    max_fan_in = max(fan_in.values())
    min_fan_in = min(fan_in.values())
    assert max_fan_in == min_fan_in == 0

    return topo_order


def dag_longest_path(adj, source, sink):
    """
    Returns the longest path, and its length, between two given nodes in a weighted DAG.
    :param adj: The adjacency lists for the weighted DAG.
    :param source: The starting point for the longest path.
    :param sink: The end point for the longest path.
    :return: A pair, with the length of the longest path as first item and the path as second item; the path is a sequence of vertices.
    """

    ''' Traverse the vertices in topological order, from source to sink. Note that not all vertices in the graph will necessarily be traversed. '''
    topo_order = topological_ordering(adj)
    start_idx = topo_order.index(source)
    ''' Backtracking information, will be used to reconstruct the longest path found from source to sink, walking it backward. '''
    backtrack = {}
    for current in topo_order[start_idx:]:
        if current == sink:
            break
        ''' Check all vertices adjacent to the one being traversed, and fill in/update their backtracking information as needed.'''
        adjs = adj.get(current, None)
        if adjs is None:
            continue
        ''' Get the max distance currently known from source to the current vertex. If the current vertex is the source, then the max distance is 0; if it is not the source and the max distance is not known yet, then it must be a vertex with no incoming edges, and the max distance is set to -infinity.'''
        _, current_max_dist = backtrack.get(current, (None, None))
        if current_max_dist is None:
            current_max_dist = 0 if current == source else float('-inf')
        for vertex, distance in adjs:
            _, vertex_max_dist = backtrack.get(vertex, (None, None))
            if vertex_max_dist is None or vertex_max_dist < current_max_dist + distance:
                backtrack[vertex] = (current, current_max_dist + distance)

    ''' Reconstruct the longest path from source to sink by walking the backtracking information backward, from sink to source. '''
    longest_path = [sink]
    current = sink
    while current != source:
        current, _ = backtrack.get(current, (None, None))
        if current is None:
            longest_path = None
            break
        longest_path.append(current)
    if longest_path is not None:
        longest_path.reverse()

    # Get the length of the longest path
    _, path_length = backtrack.get(sink, (None, None))

    return path_length, longest_path


def longest_common_string(string1, string2):
    """
    Returns the longest common string between two given strings.
    :param string1: The first given string, a string.
    :param string2: The second given string, a string.
    :return: The longest common string, a string. If multiple such strings exist, then only one of them is returned.
    """

    def vertex_name(row, col):
        """
        Converts a pair of integer numbers into a string, and returns the string. It is given by the two numbers separated by a comma.
        :param row: The first integer.
        :param col: The second integer.
        :return: The string.
        """
        name = str(row) + ',' + str(col)
        return name

    def vertex_from_name(name):
        """
        Given a string representing two integer numbers separated by a comma, returns the two numbers.
        :param name:  The given string.
        :return: A pair of two integer numbers.
        """
        r, c = name.split(',')
        row, col = int(r), int(c)
        return row, col

    # Convert the two strings to be matched into their alignment graph
    adj = {}  # The aligment graph will go here
    for row_i, row_item in enumerate(string1):
        adj[vertex_name(row_i, len(string2))] = [(vertex_name(row_i + 1, len(string2)), 0)]
        for col_i, col_item in enumerate(string2):
            adj[vertex_name(row_i, col_i)] = [(vertex_name(row_i, col_i + 1), 0),
                                              (vertex_name(row_i + 1, col_i), 0),
                                              (vertex_name(row_i + 1, col_i + 1), 1 if row_item == col_item else 0)]

    for col_i in range(0, len(string2)):
        adj[vertex_name(len(string1), col_i)] = [(vertex_name(len(string1), col_i + 1), 0)]

    # Determine the longest path along the alignment graph
    length, longest_path = dag_longest_path(adj, source=vertex_name(0, 0), sink=vertex_name(len(string1), len(string2)))

    # Convert the longest path into the longest common string
    longest_string = []
    for i in range(0, len(longest_path) - 1):
        vertex1 = longest_path[i]
        vertex2 = longest_path[i + 1]
        v1_r, v1_c = vertex_from_name(vertex1)
        v2_r, v2_c = vertex_from_name(vertex2)
        if v2_c - v1_c == v2_r - v1_r == 1 and string1[v1_r] == string2[v1_c]:
            longest_string.append(string1[v1_r])

    longest_string = ''.join(longest_string)

    # Verify that the length of the longest common string is the same as the length of the longest path
    assert length == len(longest_string)

    return longest_string
