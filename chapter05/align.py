import numpy as np

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


def vertex_name(row, col, layer=None):  # TODO fix the documentation
    """
    Converts a pair of integer numbers into a string, and returns the string. It is given by the two numbers separated by a comma.
    :param row: The first integer.
    :param col: The second integer.
    :return: The string.
    """
    name = str(row) + ',' + str(col)
    if layer is not None:
        name = name + ',' + str(layer)
    return name


def vertex_from_name(name):
    """
    Given a string representing two integer numbers separated by a comma, returns the two numbers.
    :param name:  The given string.
    :return: A pair of two integer numbers.
    """
    items = name.split(',')
    if len(items) == 2:
        r, c = items
        row, col = int(r), int(c)
        return row, col
    assert len(items) == 3
    r, c, l = items
    row, col, layer = int(r), int(c), int(l)
    return row, col, layer


def scoring_matrix_as_dict(alphabet, scoring_matrix):
    """
    Returns a scoring matrix as a dictionary of dictionaries.
    :param alphabet: The alphabet used by the scoring matrix, and to access its elements; a sequence of strings.
    :param scoring_matrix: The scoring matrix, a sequence of sequences of numbers.
    :return: The scoring matrix, organised in a dictionary of dictionaries. E.g., if scoring_m is the return value, the scoring_m['A']['B'] accesses the score for the pair 'A', 'B'.
    """
    as_dict = {}
    for row_i in range(0, len(alphabet)):
        row_dict = {}
        for col_i in range(0, len(alphabet)):
            row_dict[alphabet[col_i]] = scoring_matrix[row_i][col_i]
        as_dict[alphabet[row_i]] = row_dict
    return as_dict


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
    scoring_matrix = scoring_matrix_as_dict(alphabet, scoring_matrix)

    # TODO Making of the alignment graph could be vastly simplified by the approach used in three_way_alignment()
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


def alignment_from_longest_path(best_path, scores, string1, string2, local=False):
    """
    Returns the alignment between two strings, given a path in their alignment graph.
    :param best_path: The path in the alignment graph; a sequence of pairs, where each pair corresponds to a vertex; the first element of each pair being the row index, and the second element the column index.
    :param scores: A sequence of numbers; the i-th item is the overall score to traverse the graph along the path from source to the i-th vertex of the path. Because the last vertex of the path is the sink of the graph, the last item in the sequence is the overall score for the path.
    :param string1: The first string.
    :param string2: The second string.
    :param local: A boolean, indicating if local alignment is sought; if False, then global alignment is sought.
    :return: A pair of strings, the alignment strings corresponding respectively to string1 and string2; any insertion and deletion in the strings is represented by a "-" character.
    """
    blank = '-'
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
    return aligned1, aligned2


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
    adj = alignment_graph_from_strings(string1, string2, scoring_matrix, alphabet, sigma, local)

    # Determine the longest path along the alignment graph
    best_path, scores = dag_longest_path(adj, source=vertex_name(0, 0), sink=vertex_name(len(string1), len(string2)))

    aligned1, aligned2 = alignment_from_longest_path(best_path, scores, string1, string2, local)

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
                [-2, -3, 1, 0, -3, 0, 1, -3, 0, -3, -2, 6, -2, 0, 0, 1, 0, -3, -4, -2],
                [-1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2, 7, -1, -2, -1, -1, -2, -4, -3],
                [-1, -3, 0, 2, -3, -2, 0, -3, 1, -2, 0, 0, -1, 5, 1, 0, -1, -2, -2, -1],
                [-1, -3, -2, 0, -3, -2, 0, -3, 2, -2, -1, 0, -2, 1, 5, -1, -1, -3, -3, -2],
                [1, -1, 0, 0, -2, 0, -1, -2, 0, -2, -1, 1, -1, 0, -1, 4, 1, -2, -3, -2],
                [0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1, 0, -1, -1, -1, 1, 5, 0, -2, -2],
                [0, -1, -3, -2, -1, -3, -3, 3, -2, 1, 1, -3, -2, -2, -3, -2, 0, 4, -3, -1],
                [-3, -2, -4, -3, 1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11, 2],
                [-2, -2, -3, -2, 3, -3, 2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1, 2, 7]]
    return alphabet, blosum62


def score_and_check(input1, input2, solution1, solution2, scoring_matrix, alphabet, sigma, local=False):
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
    scoring_matrix = scoring_matrix_as_dict(alphabet, scoring_matrix)
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

    if local:
        return score, solution1_compacted in input1 and solution2_compacted in input2

    return score, solution1_compacted == input1 and solution2_compacted == input2


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


def fit_align(long, short):
    """
    Returns a highest score fitting alignment of one string into another, where matches are scored 1, and both mismatches and indels are penalised by 1.
    :param long: The string to be aligned against.
    :param short: The string to be aligned against the long string.
    :return: A pair, the first element is the score of the alignment, the second element is a pair of strings, with the alignment of long and short respectively.
    """
    assert len(long) >= len(short)

    # Determine the alphabet of the two strings, and store it in a list. The list is sorted to facilitate debugging
    alphabet = list(sorted(set(long + short)))

    # Set the scoring matrix
    scoring_matrix = [[-1 if row_item != col_item else 1 for col_item in alphabet] for row_item in alphabet]

    # Compute the alignment graph for the two given strings and the scoring matrix
    adj = alignment_graph_from_strings(long, short, scoring_matrix, alphabet, sigma=1, local=False)

    # Add free taxi-rides (local aligment) for long only

    for row_i in range(1, len(long) + 1):
        # source -> (row_i, 0)
        adj[vertex_name(0, 0)].append((vertex_name(row_i, 0), 0))

    for row_i in range(0, len(long)):
        # (row_1, starboard side) -> sink
        adj[vertex_name(row_i, len(short))].append((vertex_name(len(long), len(short)), 0))

    # Find the longest path in the DAG
    best_path, scores = dag_longest_path(adj, source=vertex_name(0, 0), sink=vertex_name(len(long), len(short)))

    # Compute the alignment for the two strings, given the longest path
    aligned1, aligned2 = alignment_from_longest_path(best_path, scores, long, short, local=True)

    return scores[-1], (aligned1, aligned2)


def overlap_align(string1, string2):
    """
    Returns the highest scoring overlap alignment between two strings, and their score, where mismatches and deletions are penalised by 2, and matches are scored 1. The operation is not commutative.
    :param string1: The first string.
    :param string2: The second string.
    :return: A pair, the first element is the score, the second element is a pair of strings, giving the alignment of a suffix of string1 and a prefix of string2.
    """
    match_score = 1
    mismatch_score = -2
    deletion_score = -2

    # Determine the alphabet of the two strings, and store it in a list. The list is sorted to facilitate debugging
    alphabet = list(sorted(set(string1 + string2)))

    # Set the scoring matrix
    scoring_matrix = [[mismatch_score if row_item != col_item else match_score for col_item in alphabet] for row_item in
                      alphabet]

    # Compute the alignment graph for the two given strings and the scoring matrix
    adj = alignment_graph_from_strings(string1, string2, scoring_matrix, alphabet, sigma=-deletion_score, local=False)

    # Add free taxi-rides (local aligment) to allow skipping a prefix of string1 and a postfix of string2

    for row_i in range(1, len(string1) + 1):
        # source -> (row_i, 0)
        adj[vertex_name(0, 0)].append((vertex_name(row_i, 0), 0))

    for col_i in range(0, len(string2)):
        # (bottom side, col_i) -> sink
        adj[vertex_name(len(string1), col_i)].append((vertex_name(len(string1), len(string2)), 0))

    # Find the longest path in the DAG
    best_path, scores = dag_longest_path(adj, source=vertex_name(0, 0), sink=vertex_name(len(string1), len(string2)))

    # Compute the alignment for the two strings, given the longest path
    aligned1, aligned2 = alignment_from_longest_path(best_path, scores, string1, string2, local=True)

    return scores[-1], (aligned1, aligned2)


def align_with_gap_penalties(ammino1, ammino2, alphabet, scoring_matrix, gap_open_penalty, gap_ext_penalty):
    """
    Returns the highest scoring alignment between two strings, using affine gap penalties, and its score.
    :param ammino1: The first string.
    :param ammino2: The second string.
    :param alphabet: The strings and scoring matrix alphabet, a sequence of strings.
    :param scoring_matrix: The scoring matrix, a sequence of sequences of strings.
    :param gap_open_penalty: The gap opening penalty, a number, usually positive; if negative, it rewards opening a gap.
    :param gap_ext_penalty: The gap extension penalty, a number, usually positive; if negative, it rewards extending a gap.
    :return: The score of the highest scoring alignment, a number, followed by a pair of strings with the alignment. A '-' character is used in the alignment to indicate and isertion or deletion.
    """
    scoring_matrix = scoring_matrix_as_dict(alphabet, scoring_matrix)
    # Convert the two strings to be matched into their alignment graph
    adj = {}  # The aligment graph will go here
    for row_i, row_item in enumerate(ammino1):
        # Fill in the starboard side, for row row_i
        adj[vertex_name(row_i, len(ammino2), 0)] = [(vertex_name(row_i + 1, len(ammino2), 0), -gap_ext_penalty),  # Down
                                                    (vertex_name(row_i, len(ammino2), 1), 0)]  # To layer 1
        adj[vertex_name(row_i, len(ammino2), 2)] = [(vertex_name(row_i, len(ammino2), 1), 0)]  # To layer 1
        adj[vertex_name(row_i, len(ammino2), 1)] = [
            (vertex_name(row_i + 1, len(ammino2), 0), -gap_open_penalty)]  # To layer 0, -sigma

        # Fill in the rest of row row_i
        for col_i, col_item in enumerate(ammino2):
            adj[vertex_name(row_i, col_i, 0)] = [(vertex_name(row_i + 1, col_i, 0), -gap_ext_penalty),  # Down, -epsilon
                                                 (vertex_name(row_i, col_i, 1), 0)]  # To layer 1, cost 0
            adj[vertex_name(row_i, col_i, 1)] = [
                (vertex_name(row_i + 1, col_i + 1, 1), scoring_matrix[ammino2[col_i]][ammino1[row_i]]),  # Down-right
                (vertex_name(row_i + 1, col_i, 0), -gap_open_penalty),  # To layer 0, cost -sigma
                (vertex_name(row_i, col_i + 1, 2), -gap_open_penalty)]  # To layer 2, cost -sigma
            adj[vertex_name(row_i, col_i, 2)] = [(vertex_name(row_i, col_i + 1, 2), -gap_ext_penalty),
                                                 # right, -epsilon
                                                 (vertex_name(row_i, col_i, 1), 0)]  # To layer 1, cost 0

    # Complete the bottom
    for col_i in range(0, len(ammino2)):
        assert adj.get(vertex_name(len(ammino1), col_i, 0)) is None
        adj[vertex_name(len(ammino1), col_i, 0)] = [(vertex_name(len(ammino1), col_i, 1), 0)]  # To layer 1, cost 0
        assert adj.get(vertex_name(len(ammino1), col_i, 2)) is None
        adj[vertex_name(len(ammino1), col_i, 2)] = [(vertex_name(len(ammino1), col_i + 1, 2), -gap_ext_penalty),
                                                    # Right, -epsilon
                                                    (vertex_name(len(ammino1), col_i, 1), 0)]  # To layer 1, cost 0

    # Complete the bottom-right corner
    assert adj.get(vertex_name(len(ammino1), len(ammino2), 0)) is None
    assert adj.get(vertex_name(len(ammino1), len(ammino2), 2)) is None
    adj[vertex_name(len(ammino1), len(ammino2), 0)] = [(vertex_name(len(ammino1), len(ammino2), 1), 0)]
    adj[vertex_name(len(ammino1), len(ammino2), 2)] = [(vertex_name(len(ammino1), len(ammino2), 1), 0)]

    # Find the longest path in the DAG
    best_path, scores = dag_longest_path(adj, source=vertex_name(0, 0, 1),
                                         sink=vertex_name(len(ammino1), len(ammino2), 1))

    # Construct the aligned strings based on the longest path

    blank = '-'
    aligned1, aligned2 = [], []
    for i in range(0, len(best_path) - 1):
        vertex1 = best_path[i]
        vertex2 = best_path[i + 1]
        v1_r, v1_c, v1_l = vertex_from_name(vertex1)
        v2_r, v2_c, v2_l = vertex_from_name(vertex2)
        # Match and mismatch
        if v1_l == v2_l == 1:
            aligned1.append(ammino1[v1_r])
            aligned2.append(ammino2[v1_c])
        # Deletion
        elif v1_l == v2_l == 0 or (v1_l == 1 and v2_l == 0):
            aligned1.append(ammino1[v1_r])
            aligned2.append(blank)
        # Insertion
        elif v1_l == v2_l == 2 or (v1_l == 1 and v2_l == 2):
            aligned1.append(blank)
            aligned2.append(ammino2[v1_c])

    aligned1 = ''.join(aligned1)
    aligned2 = ''.join(aligned2)

    best_score = scores[-1]

    return best_score, (aligned1, aligned2)


def middle_edge(string1, string2, alphabet, scoring_matrix, sigma):
    """
    Returns the middle edge in the alignment graph between two strings, and its score.
    :param string1: The first string.
    :param string2: The second string.
    :param alphabet: The strings alphabet.
    :param scoring_matrix: The scoring matrix to be used in the alignment graph, a sequence of sequences of numbers.
    :param sigma: The penalty for insertions and deletions, a positive number.
    :return: A pair (score, ((v1_row, v1_col), (v2_row, v2_col)) where score is the middle edge score, and the next element of the pair is the middle edge; (v1_row, v1_col) is the vertex from where the edge originates, and (v2_row, v2_col) is the vertex where the edge goes into.
    """

    def split(string):
        if len(string) >= 2:
            col = len(string) // 2
            return string[:col], string[-(len(string) - col):][::-1]
        if len(string) == 1:
            return '', string[0]
        return '', ''

    def last_2columns_score(string1, string2, alphabet, scoring_matrix, sigma):
        scoring_matrix = scoring_matrix_as_dict(alphabet, scoring_matrix)

        prev_scores = None
        col_scores = [i * (-sigma) for i in range(0, len(string1) + 1)]
        for col_i in range(1, len(string2) + 1):
            prev_scores = col_scores
            col_scores = [prev_scores[0] - sigma]
            for row in range(1, len(string1) + 1):
                col_scores.append(max((col_scores[-1] - sigma,
                                       prev_scores[row] - sigma,
                                       prev_scores[row - 1] + scoring_matrix[string2[col_i - 1]][string1[row - 1]]
                                       ))
                                  )

        return prev_scores, col_scores

    # Select the middle column of the alignment graph
    col = len(string2) // 2
    string2_l, string2_r = split(string2)

    # Get the scores from the source to the middle column
    _, scores_l = last_2columns_score(string1,
                                      string2_l,
                                      alphabet,
                                      scoring_matrix,
                                      sigma)

    # Get the scores from the sink to the middle column
    prev_scores_r, scores_r = last_2columns_score(string1[::-1],
                                                  string2_r,
                                                  alphabet,
                                                  scoring_matrix,
                                                  sigma)
    scores_r.reverse()
    if prev_scores_r is not None:
        prev_scores_r.reverse()

    scoring_matrix = scoring_matrix_as_dict(alphabet, scoring_matrix)

    # Find the edge with the highest score originating from the middle column

    highest_score = float('-inf')
    vertex1 = None  # The vertex from where the edge originates (in the middle column)
    vertex2 = None  # The vertex into which the edge goes
    for row in range(0, len(string1) + 1):
        # Check the edge going to the right, from a vertex in the middle column to a vertex in the next column on the same row
        if col < len(
                string2):  # You cannot go to the next column to the right if you are already on the rightmost columns
            score_right = scores_l[row] - sigma + prev_scores_r[row]
            if score_right > highest_score:
                highest_score = score_right
                vertex1 = (row, col)
                vertex2 = (row, col + 1)
        # Check the diagonal edge and the edge going straight down
        if row < len(string1) and col < len(
                string2):  # You cannot go diagonal if already on the bottom row or at the rightmost column
            score_diagonal = scores_l[row] + scoring_matrix[string2[col]][string1[row]] + prev_scores_r[row + 1]
            if score_diagonal > highest_score:
                highest_score = score_diagonal
                vertex1 = (row, col)
                vertex2 = (row + 1, col + 1)
        if row < len(string1):  # You cannot go one row lower if you are already on the bottom row
            score_vertical = scores_l[row] - sigma + scores_r[row + 1]
            if score_vertical > highest_score:
                highest_score = score_vertical
                vertex1 = (row, col)
                vertex2 = (row + 1, col)
    return highest_score, (vertex1, vertex2)


def linear_space_alignment(string1, string2, alphabet, scoring_matrix, sigma):
    """
    Returns the score of the highest score global alignment between two strings, and the path through their alignment graph.
    :param string1:The first string.
    :param string2: The second string.
    :param alphabet: The strings alphabet.
    :param scoring_matrix: The scoring matrix, a sequence of sequences of strings.
    :param sigma: The penalty for insertions and deletions, a positive number.
    :return: The score of the best global alignment between the two strings, a number, followed by the corresponding path in the alignment graph, a string; the string provides the path from source to sink, and may contain the following characters: '-' to indicate an edge going to the right in the graph, '|' to indicate an adge going down, and '\' to indicate an edge going diagonal, right and down.
    """
    # Find the middle edge, and the cost of the associated longest path
    middle_score, middle = middle_edge(string1, string2, alphabet, scoring_matrix, sigma)

    (v1_row, v1_col), (v2_row, v2_col) = middle

    if v1_row == v2_row and v2_col == v1_col + 1:
        middle_path = '-'
    elif v2_row == v1_row + 1 and v2_col == v1_col + 1:
        middle_path = '\\'
    elif v2_row == v1_row + 1 and v2_col == v1_col:
        middle_path = '|'
    else:
        assert False

    # If the upper-left sub-problem is non empty, call recursively on the upper-left sub-problem
    left_score, left_path = 0, ''
    if string2[:v1_col] or string2[:v1_col]:
        left_score, left_path = linear_space_alignment(string1[:v1_row], string2[:v1_col], alphabet, scoring_matrix,
                                                       sigma)
    # If the lower-right sub-problem is non empty, call recursively on the lower-right sub-problem
    right_score, right_path = 0, ''
    if string1[v2_row:] or string2[v2_col:]:
        right_score, right_path = linear_space_alignment(string1[v2_row:], string2[v2_col:], alphabet, scoring_matrix,
                                                         sigma)

    # Assemble the results
    overall_path = ''.join((left_path, middle_path, right_path))

    return middle_score, overall_path


def aligned_strings_from_path(s1, s2, path):
    """
    Returns the aligned strings from two given strings and a path through their alignment graph.
    :param s1: The first string.
    :param s2: The second string.
    :param path: The path, a string, as returned by linear_space_alignment(): a concatenation of characters indicating the movement of the path through the graph, '-' indicates an edge going right, '|' indicates and edge going down, and '\' indicates a diagonal edge going right and down.
    :return: A pair of aligned strings; the '-' character in the strings represent insertions and deletions.
    """
    blank = '-'
    aligned1, aligned2 = [], []
    pos_in_s1, pos_in_s2 = 0, 0
    for item in path:
        if item == '-':  # Insertion
            aligned2.append(s2[pos_in_s2])
            pos_in_s2 += 1
            aligned1.append(blank)
        elif item == '|':  # Deletion
            aligned1.append(s1[pos_in_s1])
            pos_in_s1 += 1
            aligned2.append(blank)
        else:  # Match and mis-match
            assert item == '\\'
            aligned1.append(s1[pos_in_s1])
            pos_in_s1 += 1
            aligned2.append(s2[pos_in_s2])
            pos_in_s2 += 1

    # Check that all characters from both input strings have been used
    assert pos_in_s1 == len(s1)
    assert pos_in_s2 == len(s2)

    aligned1 = ''.join(aligned1)
    aligned2 = ''.join(aligned2)

    return aligned1, aligned2


def three_way_alignment(string1, string2, string3):
    """
    Returns the alignment between three strings with the highest score, and its score. The score of a column of the aligned strings is 1 if all of the column's symbol are the same, 0 otherwise; the alignment score is given by the sum of the column scores across all columns.
    :param string1: The first string.
    :param string2: The second string.
    :param string3: The third string.
    :return: The alignment score, a number, followed by a tuple, with the three aligned strings; the '-' character in the strings indicate a blank.
    """
    blank = '-'

    # Make the alignment graph
    adj = {}
    increments = ((0, 0, 1),
                  (0, 1, 0),
                  (1, 0, 0),
                  (0, 1, 1),
                  (1, 1, 0),
                  (1, 0, 1),
                  (1, 1, 1))
    for row_i in range(0, len(string1) + 1):
        for col_i in range(0, len(string2) + 1):
            for depth_i in range(0, len(string3) + 1):
                adjs = []
                for incr in increments:
                    adj_vertex = tuple(map(sum, zip((row_i, col_i, depth_i), incr)))
                    if adj_vertex[0] > len(string1) or adj_vertex[1] > len(string2) or adj_vertex[2] > len(string3):
                        continue
                    weight = 1 if (incr == (1, 1, 1)) and (string1[row_i] == string2[col_i] == string3[depth_i]) else 0
                    adjs.append((vertex_name(*adj_vertex), weight))
                current_vertex = vertex_name(row_i, col_i, depth_i)
                assert adj.get(current_vertex) is None
                if adjs:
                    adj[current_vertex] = adjs

    # Find the longest path in the alignment graph
    path, scores = dag_longest_path(adj,
                                    vertex_name(0, 0, 0),
                                    vertex_name(len(string1), len(string2), len(string3)))

    # Build the strings alignment based on the longest path
    prev_r, prev_c, prev_d = 0, 0, 0
    aligned1, aligned2, aligned3 = [], [], []
    for v_name in path[1:]:
        r, c, d = vertex_from_name(v_name)
        aligned1.append(string1[prev_r] if r > prev_r else blank)
        aligned2.append(string2[prev_c] if c > prev_c else blank)
        aligned3.append(string3[prev_d] if d > prev_d else blank)
        prev_r, prev_c, prev_d = r, c, d

    aligned1 = ''.join(aligned1)
    aligned2 = ''.join(aligned2)
    aligned3 = ''.join(aligned3)
    assert len(aligned1) == len(aligned2) == len(aligned3)

    return scores[-1], (aligned1, aligned2, aligned3)
