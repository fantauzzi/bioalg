import numpy as np
from pandas import DataFrame

""" 
NOTE
====

The code uses adjacency lists to represent weighted directed acyclic graphs (DAGs). Adjacency lists are a Python dictionary. Keys in the dictionary are vertices, and values are lists of pairs. Every pair gives a vertex adjacent to the vertex of the key, and the weight of the conencting edge. Vertices are named with strings. Weights are non-negative integer numbers. If a vertex has no outgoing edges, but only incoming edges, then it may or may not appear among the keys.
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
    :return: A pair (path, distances), where path is a list of vertices (strings) giving the longest path from source to sink included, and distances is a list of numbers, giving the respective distances from the source to vertices in path.
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
    distances = []
    current = sink
    while current != source:
        current, distance = backtrack.get(current, (None, None))
        if current is None:
            longest_path = None
            distances = None
            break
        longest_path.append(current)
        distances.append(distance)
    if longest_path is not None:
        longest_path.reverse()
        distances.reverse()

    # Get the length of the longest path
    _, path_length = backtrack.get(sink, (None, None))
    assert distances[-1] == path_length

    return longest_path, distances


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


def alignment_graph_from_strings(string1, string2, scoring_matrix, alphabet, sigma, local=False):
    """
    Returns the alignment graph (a weighted DAG) for two given strings.
    :param string1: The first string.
    :param string2: The second string.
    :param scoring_matrix: The scoring matrix to be used, a sequence of sequences of numbers.
    :param alphabet: The alphabet used by the strings, a sequence of strings whose items order is matched by the order of rows and columns in the scoring_matrix.
    :param sigma: The penalty to be used for insertions and deletions, a number, usually positive; if negative, it becomes a reward for indels.
    :param local: A boolean, indicating whether "free rides" should be allowed in the aligment graph, i.e. if a graph for local alignment is required.
    :return: The adjacency lists of the weighted DAG, a dictionary.
    """

    scoring_matrix = DataFrame(scoring_matrix, columns=alphabet, index=alphabet)
    # Convert the two strings to be matched into their alignment graph
    adj = {}  # The aligment graph will go here
    for row_i, row_item in enumerate(string1):
        adj[vertex_name(row_i, len(string2))] = [(vertex_name(row_i + 1, len(string2)), -sigma)]
        for col_i, col_item in enumerate(string2):
            adj[vertex_name(row_i, col_i)] = [(vertex_name(row_i, col_i + 1), -sigma),
                                              (vertex_name(row_i + 1, col_i), -sigma),
                                              (vertex_name(row_i + 1, col_i + 1),
                                               scoring_matrix[string2[col_i]][string1[row_i]])]
            if local:
                # (row_i, col_i) -> sink
                adj[vertex_name(row_i, col_i)].append((vertex_name(len(string1), len(string2)), 0))
                # source -> (row_i, col_i)
                if (row_i, col_i) != (0, 0):  # Avoid to add a loop around source
                    adj[vertex_name(0, 0)].append((vertex_name(row_i, col_i), 0))

    for col_i in range(0, len(string2)):
        assert adj.get(vertex_name(len(string1), col_i)) is None
        adj[vertex_name(len(string1), col_i)] = [(vertex_name(len(string1), col_i + 1), -sigma)]
        if local:
            # (bottom row, col_i) -> sink
            adj[vertex_name(len(string1), col_i)].append((vertex_name(len(string1), len(string2)), 0))
            # source -> (bottom_row, col_i)
            adj[vertex_name(0, 0)].append((vertex_name(len(string1), col_i), 0))

    if local:
        for row_i in range(0, len(string1)):
            # source -> (row_i, starboard side)
            adj[vertex_name(0, 0)].append((vertex_name(row_i, len(string2)), 0))
            # (row_1, starboard side) -> sink
            adj[vertex_name(row_i, len(string2))].append((vertex_name(len(string1), len(string2)), 0))

    return adj


def longest_common_string(string1, string2):
    """
    Returns the longest common string between two given strings.
    :param string1: The first given string, a string.
    :param string2: The second given string, a string.
    :return: The longest common string, a string. If multiple such strings exist, then only one of them is returned.
    """

    # Determine the alphabet of the two strings, and store it in a list. The list is sorted to facilitate debugging.
    alphabet = list(sorted(set(string1 + string2)))

    # Determine the scoring matrix to be used, which is just an identity matrix
    scoring_matrix = [[1 if row_item == col_item else 0 for col_item in alphabet] for row_item in alphabet]
    adj = alignment_graph_from_strings(string1, string2, scoring_matrix, alphabet, sigma=0, local=False)

    # Determine the longest path along the alignment graph
    longest_path, lengths = dag_longest_path(adj, source=vertex_name(0, 0),
                                             sink=vertex_name(len(string1), len(string2)))
    length = lengths[-1]

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


def best_alignment(string1, string2, scoring_matrix, alphabet, sigma, local=False):
    """
    Returns the best global or local alignment between two strings, based on a scoring matrix.
    :param string1: The first string.
    :param string2: The second string.
    :param scoring_matrix: The scoring matrix, a sequence of sequences of integer numbers.
    :param alphabet: All the characters that could appear in either string, a sequence of strings.
    :param sigma: The penalty to be applied to insertions and deletion; an integer, typically positive; if negative, it becomes a reward for indels.
    :param local: True if a local alignment is sought, False for a global alignment.
    :return:the score of the alignment and a pair of strings, with the calculated alignment; insertions and deletions are indicated in the strings by a '-'.
    """
    blank = '-'

    adj = alignment_graph_from_strings(string1, string2, scoring_matrix, alphabet, sigma, local)

    # Determine the longest path along the alignment graph
    best_path, scores = dag_longest_path(adj, source=vertex_name(0, 0), sink=vertex_name(len(string1), len(string2)))

    # Convert the longest path into the alignment of the two given strings
    aligned1, aligned2 = [], []
    for i in range(0, len(best_path) - 1):
        vertex1 = best_path[i]
        vertex2 = best_path[i + 1]
        v1_r, v1_c = vertex_from_name(vertex1)
        v2_r, v2_c = vertex_from_name(vertex2)
        if local and ((scores[i] == 0 and (v1_r, v1_c) == (0, 0)) or (
                scores[i] == scores[-2] and (v2_r, v2_c) == (len(string1), len(string2)))):
            continue
        if v2_r == v1_r and v2_c == v1_c + 1:
            aligned1.append(blank)
            aligned2.append(string2[v1_c])
        elif v2_r == v1_r + 1 and v2_c == v1_c:
            aligned2.append(blank)
            aligned1.append(string1[v1_r])
        else:
            assert local or (v2_r == v1_r + 1 and v2_c == v1_c + 1)
            aligned1.append(string1[v1_r])
            aligned2.append(string2[v1_c])

    aligned1 = ''.join(aligned1)
    aligned2 = ''.join(aligned2)

    return scores[-1], (aligned1, aligned2)


def get_pam250():
    """
    Returns the alphabet and matrix of a PAM 250 alignment matrix.
    :return: A pair, the first element is the alphabet, a sequence of strings, and the second element is the matrix, a sequence of sequences of numbers. Order of the alphabet elements matches the order of rows and columns in the matrix.
    """
    alphabet = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    pam250 = [[2, -2, 0, 0, -3, 1, -1, -1, -1, -2, -1, 0, 1, 0, -2, 1, 1, 0, -6, -3],
              [-2, 12, -5, -5, -4, -3, -3, -2, -5, -6, -5, -4, -3, -5, -4, 0, -2, -2, -8, 0],
              [0, -5, 4, 3, -6, 1, 1, -2, 0, -4, -3, 2, -1, 2, -1, 0, 0, -2, -7, -4],
              [0, -5, 3, 4, -5, 0, 1, -2, 0, -3, -2, 1, -1, 2, -1, 0, 0, -2, -7, -4],
              [-3, -4, -6, -5, 9, -5, -2, 1, -5, 2, 0, -3, -5, -5, -4, -3, -3, -1, 0, 7],
              [1, -3, 1, 0, -5, 5, -2, -3, -2, -4, -3, 0, 0, -1, -3, 1, 0, -1, -7, -5],
              [-1, -3, 1, 1, -2, -2, 6, -2, 0, -2, -2, 2, 0, 3, 2, -1, -1, -2, -3, 0],
              [-1, -2, -2, -2, 1, -3, -2, 5, -2, 2, 2, -2, -2, -2, -2, -1, 0, 4, -5, -1],
              [-1, -5, 0, 0, -5, -2, 0, -2, 5, -3, 0, 1, -1, 1, 3, 0, 0, -2, -3, -4],
              [-2, -6, -4, -3, 2, -4, -2, 2, -3, 6, 4, -3, -3, -2, -3, -3, -2, 2, -2, -1],
              [-1, -5, -3, -2, 0, -3, -2, 2, 0, 4, 6, -2, -2, -1, 0, -2, -1, 2, -4, -2],
              [0, -4, 2, 1, -3, 0, 2, -2, 1, -3, -2, 2, 0, 1, 0, 1, 0, -2, -4, -2],
              [1, -3, -1, -1, -5, 0, 0, -2, -1, -3, -2, 0, 6, 0, 0, 1, 0, -1, -6, -5],
              [0, -5, 2, 2, -5, -1, 3, -2, 1, -2, -1, 1, 0, 4, 1, -1, -1, -2, -5, -4],
              [-2, -4, -1, -1, -4, -3, 2, -2, 3, -3, 0, 0, 0, 1, 6, 0, -1, -2, 2, -4],
              [1, 0, 0, 0, -3, 1, -1, -1, 0, -3, -2, 1, 1, -1, 0, 2, 1, -1, -2, -3],
              [1, -2, 0, 0, -3, 0, -1, 0, 0, -2, -1, 0, 0, -1, -1, 1, 3, 0, -5, -3],
              [0, -2, -2, -2, -1, -1, -2, 4, -2, 2, 2, -2, -1, -2, -2, -1, 0, 4, -6, -2],
              [-6, -8, -7, -7, 0, -7, -3, -5, -3, -2, -4, -4, -6, -5, 2, -2, -5, -6, 17, 0],
              [-3, 0, -4, -4, 7, -5, 0, -1, -4, -1, -2, -2, -5, -4, -4, -3, -3, -2, 0, 10]]
    return alphabet, pam250


def get_blosum62():
    """
    Returns the alphabet and matrix of a BLOSUM 62 alignment matrix.
    :return: A pair, the first element is the alphabet, a sequence of strings, and the second element is the matrix, a sequence of sequences of numbers. Order of the alphabet elements matches the order of rows and columns in the matrix.
    """
    alphabet = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    blosum62 = [[4, 0, -2, -1, -2, 0, -2, -1, -1, -1, -1, -2, -1, -1, -1, 1, 0, 0, -3, -2],
                      [0, 9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2],
                      [-2, -3, 6, 2, -3, -1, -1, -3, -1, -4, -3, 1, -1, 0, -2, 0, -1, -3, -4, -3],
                      [-1, -4, 2, 5, -3, -2, 0, -3, 1, -3, -2, 0, -1, 2, 0, 0, -1, -2, -3, -2],
                      [-2, -2, -3, -3, 6, -3, -1, 0, -3, 0, 0, -3, -4, -3, -3, -2, -2, -1, 1, 3],
                      [0, -3, -1, -2, -3, 6, -2, -4, -2, -4, -3, 0, -2, -2, -2, 0, -2, -3, -2, -3],
                      [-2, -3, -1, 0, -1, -2, 8, -3, -1, -3, -2, 1, -2, 0, 0, -1, -2, -3, -2, 2],
                      [-1, -1, -3, -3, 0, -4, -3, 4, -3, 2, 1, -3, -3, -3, -3, -2, -1, 3, -3, -1],
                      [-1, -3, -1, 1, -3, -2, -1, -3, 5, -2, -1, 0, -1, 1, 2, 0, -1, -2, -3, -2],
                      [-1, -1, -4, -3, 0, -4, -3, 2, -2, 4, 2, -3, -3, -2, -2, -2, -1, 1, -2, -1],
                      [-1, -1, -3, -2, 0, -3, -2, 1, -1, 2, 5, -2, -2, 0, -1, -1, -1, 1, -1, -1],
                      [2, -3, 1, 0, -3, 0, 1, -3, 0, -3, -2, 6, -2, 0, 0, 1, 0, -3, -4, -2],
                      [-1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2, 7, -1, -2, -1, -1, -2, -4, -3],
                      [-1, -3, 0, 2, -3, -2, 0, -3, 1, -2, 0, 0, -1, 5, 1, 0, -1, -2, -2, -1],
                      [-1, -3, -2, 0, -3, -2, 0, -3, 2, -2, -1, 0, -2, 1, 5, -1, -1, -3, -3, -2],
                      [1, -1, 0, 0, -2, 0, -1, -2, 0, -2, -1, 1, -1, 0, -1, 4, 1, -2, -3, -2],
                      [0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1, 0, -1, -1, -1, 1, 5, 0, -2, -2],
                      [0, -1, -3, -2, -1, -3, -3, 3, -2, 1, 1, -3, -2, -2, -3, -2, 0, 4, -3, -1],
                      [-3, -2, -4, -3, 1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11, 2],
                      [-2, -2, -3, -2, 3, -3, 2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1, 2, 7]]
    return alphabet, blosum62


def score_and_check(input1, input2, solution1, solution2, scoring_matrix, alphabet, sigma):
    """
    Verifies if a given solution is admissible as the best local alignment between two strings, and computes and returns the related score.
    :param input1: The first string.
    :param input2: The second string.
    :param solution1: The candidate solution for the first string; insertions are marked with a '-'.
    :param solution2: The candidate solution for the second string; deletions are marked with a '-'.
    :param scoring_matrix: The scoring matrix to be used, a sequence of sequences.
    :param alphabet: The alphabet used by the strings, does not include the '-' character, used for indels; a sequence of strings.
    :param sigma: The penalty for insertions and deletions; typically a positive number, if negative, it becomes a reward for indels.
    :return: A pair with a boolean and a number; the boolean is True if and only if solution1 and solution2 are possible local alignments for input1 and input2, the number is the alignment score. If the boolean is False, the score is meaningless.
    """
    scoring_matrix = DataFrame(scoring_matrix, columns=alphabet, index=alphabet)
    s1, s2 = 0, 0
    i1, i2 = 0, 0
    score = 0
    while s1 < len(solution1) and s2 < len(solution2):
        if solution1[s1] == '-' or solution2[s2] == '-':
            score -= sigma
            s1 += 1
            s2 += 1
            if solution1[s1] == '-':
                assert solution2[s2] != '-'
                i2 += 1
            else:
                assert solution1[s1] != '-'
                i1 += 1
            continue
        score += scoring_matrix[solution2[s2]][solution1[s1]]
        s1 += 1
        s2 += 1

    solution1_compacted = solution1.replace('-', '')
    solution2_compacted = solution2.replace('-', '')

    return score, solution1_compacted in input1 and solution2_compacted in input2


def edit_distance(string1, string2):
    """
    Returns the edit distance between two strings, i.e. the number of insertions, deletions and replacements to turn one string into another. The operation is commutative.
    :param string1: The first string.
    :param string2: The second string.
    :return:The edit distance, an integer.
    """
    # Sigma is positive because it is a penalty (if negative, it would be a reward)
    sigma = 1

    # Determine the alphabet of the two strings, and store it in a list. The list is sorted to facilitate debugging
    alphabet = list(sorted(set(string1 + string2)))

    # Determine the scoring matrix to be used
    scoring_matrix = [[-1 if row_item != col_item else 0 for col_item in alphabet] for row_item in alphabet]
    adj = alignment_graph_from_strings(string1, string2, scoring_matrix, alphabet, sigma=sigma, local=False)

    # Determine the longest path along the alignment graph
    longest_path, lengths = dag_longest_path(adj, source=vertex_name(0, 0),
                                             sink=vertex_name(len(string1), len(string2)))
    length = lengths[-1]

    return -length