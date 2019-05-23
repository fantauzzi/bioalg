import matplotlib.pyplot as plt
from copy import deepcopy
from collections import namedtuple
from itertools import product, chain
import networkx as nx
# import networkx.drawing.nx_pylab as nxp
from math import log, exp
from numpy import isclose
from functools import reduce
import operator

HMM = namedtuple('HMM', ['alphabet', 'states', 'transition', 'emission'])
"""
A Hidden Markov Model, described by:
- alphabet: the alphabet of symbols that can be emitted, a sequence of strings;
- states: the hidden states, a sequence of strings;
- transition: the state transition probabilities, a dictionary of dictionaries where transition[l][k] is the probability
    of moving from state l to state k; l and k are strings belonging to 'states';
- emission: the emission probabilities, a dictionary of dictionaries where emission[k][b] is the probability of emitting
    symbol b when in state k; b is a string belonging to 'alphabet' and k a string belonging to 'states'.
"""


def make_graph(emissions, model):
    """
    Builds and return a directed graph corresponding to a HMM and a sequence of emissions, such that the longest path in the graph from source to sink traverses the sequence of states most likely to produce the given sequence of emissions.
    :param emissions: The sequence of emissions, a sequence of strings.
    :param model: The HMM, a HMM named tuple.
    :return: The built graph.
    """
    n = len(emissions)
    source = ('__source', -1)
    sink = ('__sink', n)
    graph = nx.DiGraph()

    # Add vertices to the graph (except source and sink)
    graph.add_nodes_from([(state, i_emission) for state in model.states for i_emission in range(0, n)])

    # Now add source and sink
    graph.add_node(source)
    graph.add_node(sink)

    # Add edges to the graph
    graph.add_weighted_edges_from([((state1, i_emission), (state2, i_emission + 1), weight)
                                   for i_emission in range(0, n - 1)
                                   for state1, state2 in product(model.states, model.states)
                                   for weight in (model.transition[state1][state2] * model.emission[state2][
            emissions[i_emission + 1]],)
                                   if weight != 0])
    # These are the edges going out of the source
    graph.add_weighted_edges_from(
        [(source, (state2, 0), 1 / len(model.states) * model.emission[state2][emissions[0]]) for state2 in
         model.states])
    # Finally the edges going into the sink
    graph.add_weighted_edges_from([((state1, n - 1), sink, 1) for state1 in model.states])

    return graph


def viterbi(emissions, model):
    """
    Returns the succession of hidden states in a Hidden Markov Model (HMM) with the highest probability to have produced a given succession of emissions.
    :param emissions: The emissions, a string.
    :param model: The HMM, a HMM name tuple.
    :return: The succession of hidden states, a string of the same length as 'emissions'
    """
    n = len(emissions)
    source = ('__source', -1)  # TODO replace with some sort of global
    sink = ('__sink', n)
    graph = make_graph(emissions, model)

    # Set the length of the longest path from source to the source (which is 0, of course)
    graph.nodes[source]['s'] = .0
    # Set the predecessor for source (which is None)
    graph.nodes[source]['previous'] = None
    # Iterate along the emission steps
    for i_emission in range(0, n + 1):
        # Compute the length of the longest path to source to every vertex at the current i_emission emission step
        states = model.states if i_emission < n else [sink[0]]
        for state in states:
            previous_vertices = [source] if i_emission == 0 else [(state, i_emission - 1) for state in model.states]
            highest_s = float('-inf')
            for prev_vertex in previous_vertices:
                s = graph.nodes[prev_vertex]['s'] + log(graph.edges[(prev_vertex, (state, i_emission))]['weight'])
                if s > highest_s:
                    highest_s = s
                    graph.nodes[(state, i_emission)]['s'] = highest_s
                    graph.nodes[(state, i_emission)]['previous'] = prev_vertex[0]

    # Backtrack from sink to source to determine the longest path found
    longest_path = []
    current_vertex = sink
    previous_state = ''
    while current_vertex != source:
        longest_path.append(previous_state)
        previous_state = graph.nodes[current_vertex]['previous']
        current_vertex = (previous_state, current_vertex[1] - 1)
    longest_path = ''.join(reversed(longest_path))

    return longest_path


class Counter(dict):
    """
    A dictionary that, when retrieving missing values, default them to 0. Retrieving a missing value does not add it to the dictionary.
    """

    def __missing__(self, key):
        return 0


def flatten(seq_of_seq):
    """
    Returns a sequence obtained by flattening a sequence of sequences. E.g. [[1, 2, 3], [4, 5]] becomes [1, 2, 3, 4, 5]
    :param seq_of_seq: The sequence of sequences.
    :return: The sequence after flattening, a list.
    """
    return list(chain(*seq_of_seq))


def make_profile_HMM(theta, sigma, alphabet, alignment):
    """
    Returns a profile HMM with pseudocounts from a multiple alignment.
    :param theta: The column removal threshold, a number between 0 and 1; columns in the alignment with a proportion of blanks '-' greater than or equal to the threshold are ignored.
    :param sigma: The pseudocount, a positive number; it is added to all allowed state transitions (corresponding to edges in a HMM graph) and to all emission probabilities for states that can emit (matches and insertions).
    :param alphabet: The alphabet of symbols making the alignment, a string. It must not contain the special symbol '-', reserved for blanks.
    :param alignment: The multiple alignment, a sequence of strings. Symbol '-' is used to indicate a blank.
    :return: The profile HMM, an HMM named tuple..
    """

    def threshold_alignment(alignment, theta):
        """
        Returns the set of column indices in an alignment that meet or exceed a column removal threshold.
        :param alignment: The alignment, a sequence of strings. Symbol '-' is used for blanks.
        :param theta: The removal threshold, a number between 0 and 1.
        :return: The set of column indices, numbered starting from 0, where the proportion of blanks is greater than, or equal to, the threshold theta.
        """
        n_cols = len(alignment[0])
        n_rows = len(alignment)
        alignment_T = list(map(list, zip(*alignment)))  # Transpose alignment
        # Count occurrences of '-' in the transposed of alignment, going by rows
        count_spaces = [sum([item == '-' for item in row]) for row in alignment_T]
        # The return value
        thresholded = {col for col in range(0, n_cols) if count_spaces[col] / n_rows >= theta}
        return thresholded

    def validate_prob_matrix(matrix):
        """
        Checks that the sums by row of a given matrix are all equal to 1.
        :param matrix: The matrix, a dictionary of dictionaries of numbers; matrix[r][c] is the element of row 'r' and column 'c' in the matrix.
        :raises AssertionError: exception raised when the condition is not met.
        """
        for row in matrix.values():
            row_total = sum(row.values())
            assert isclose(row_total, 1) or isclose(row_total, 0)

    def normalise_prob_matrix(matrix):
        """
        Modifies the non-zero values in a numerical matrix, such that their sum by row is equal to 1, in every row. Zero values are left unchanged.
        :param matrix: The matrix to be modified, a dictionary of dictionaries of numbers; matrix[r][c] is the element of row 'r' and column 'c' in the matrix.
        """
        for row_key, row in matrix.items():
            total = sum(row.values())
            matrix[row_key] = {col_key: 0 if item == 0 else item / total for col_key, item in row.items()}

    # The source and sink of the HMM Viterbi graph
    source = ('S', None)
    sink = ('E', None)
    # The total number of rows and columns in the alignment
    n_rows = len(alignment)
    n_cols = len(alignment[0])
    # The set of indices of alignment columns that meet or exceed the removal threshold theta
    thresholded = threshold_alignment(alignment, theta)
    # The number of columns in alignment, net of thresholded columns
    n_unthresholded = n_cols - len(thresholded)
    # The set of vertices (states) in the Viterbi graph
    states = [source, ('I', 0)] + \
             flatten([(('M', pos), ('D', pos), ('I', pos)) for pos in range(1, n_unthresholded + 1)]) + \
             [sink]

    # Counters that will be used to compute transition and emission probabilities
    node_visits = Counter()  # node_visits[a] will be the number of sequences in alignment that go through node 'a'
    edge_visits = Counter()  # edge_visits[(a,b)] will be the number of sequences in alignment that go through edge 'a'-'b'
    # node_emissions[state][symbol] will be the number of times the sequences in alignment make 'state' emit 'symbol'
    node_emissions = {state: Counter() for state in states}

    ''' For every sequence in alignement, follow the sequence through the Viterbi graph, source to sink,  and increment
    the counters '''
    for sequence in alignment:
        prev_state = ('S', None)
        pos = 0  # Numbering of the states (vertices in the Viterbi graph)
        for col, symbol in enumerate(sequence):
            assert symbol in alphabet or symbol == '-'
            if col in thresholded:
                if symbol == '-':
                    continue
                curr_state = ('I', pos)
            else:
                pos += 1
                curr_state = ('D', pos) if symbol == '-' else ('M', pos)
            node_visits[curr_state] = node_visits[curr_state] + 1
            edge_visits[(prev_state, curr_state)] = edge_visits[(prev_state, curr_state)] + 1
            if curr_state[0] in 'MI':  # Only matches and insertions can emit symbols
                node_emissions[curr_state][symbol] = node_emissions[curr_state][symbol] + 1
            prev_state = curr_state
        edge_visits[(prev_state, sink)] = edge_visits[(prev_state, sink)] + 1
        node_visits[source] = n_rows
        node_visits[sink] = n_rows

    # Compute the transition probabilities matrix based on the counters
    transition = {state_row: {
        state_col: 0 if node_visits[state_row] == 0 else edge_visits[state_row, state_col] / node_visits[state_row]
        for state_col in states}
        for state_row in states}

    validate_prob_matrix(transition)

    # Add pseudocounts to the transition probabilities matrix if requested
    if sigma:
        for i in range(0, n_unthresholded + 1):
            if i == 0:
                for node in [source, ('I', 0)]:
                    transition[node]['I', 0] = transition[node]['I', 0] + sigma
                    transition[node]['M', 1] = transition[node]['M', 1] + sigma
                    transition[node]['D', 1] = transition[node]['D', 1] + sigma
            elif i == n_unthresholded:
                for node in [('I', i), sink]:
                    transition['M', i][node] = transition['M', i][node] + sigma
                    transition['D', i][node] = transition['D', i][node] + sigma
                    transition['I', i][node] = transition['I', i][node] + sigma
            else:
                for node1 in [('M', i), ('D', i), ('I', i)]:
                    for node2 in [('I', i), ('M', i + 1), ('D', i + 1)]:
                        transition[node1][node2] = transition[node1][node2] + sigma
        normalise_prob_matrix(transition)
        validate_prob_matrix(transition)

    # Compute the emission probabilities matrix based on the counters
    emission = {state: {symbol: 0 if total == 0 else node_emissions[state][symbol] / total
                        for total in (sum(node_emissions[state].values()),)
                        for symbol in alphabet}
                for state in states}

    validate_prob_matrix(emission)

    # Add pseudocounts to the emission probabilities matrix if requested
    if sigma:
        for row_key, row in emission.items():
            # Pseudo-counts only apply to states that can emit, that is Insertion ('I') and Match ('M')
            if row_key[0] in 'MI':
                emission[row_key] = {col_key: item + sigma for col_key, item in row.items()}
        normalise_prob_matrix(emission)

    the_hmm = HMM(alphabet=alphabet, states=states, transition=transition, emission=emission)
    return the_hmm


def make_profile_graph(emissions, model):
    """
    Builds and returns the Viterbi graph to score a text against a given profile HMM, and a topological order of its nodes.
    :param emissions: The text, a string.
    :param model: The profile HMM, an HMM named tuple.
    :return: A pair with the graph and its nodes in topological order. The graph is a nx.DiGraph; the topological order is a list of nodes. Weights of the graph edges are probabilities (not logs of probabilities).
    """
    n = len(emissions)
    n_rows = max([index for _, index in model.states if index is not None])
    graph = nx.DiGraph()  # The graph bein built
    # The HMM states except source and sink
    states = deepcopy(model.states)
    states.remove(('S', None))
    states.remove(('E', None))

    # Compile the list of nodes in topological order
    nodes = [('D', row, 0) for row in range(0, n_rows + 1)]
    nodes.extend([(state, row_i, col_i) for col_i in range(1, n + 1) for state, row_i in states])
    nodes.append(('E', n_rows, n + 1))

    # Time to add them to the graph
    graph.add_nodes_from(nodes)

    # Compile the list of edges
    edges = []
    for node in nodes:
        if node[0] == 'E':
            continue
        state = node[:2]
        if node[2] < n:
            node2 = ('I', node[1], node[2] + 1)
            state2 = node2[:2]
            weight = model.transition[state][state2] * model.emission[state2][emissions[node2[2] - 1]] if node[
                                                                                                              2] > 0 \
                else 1 / 3 if node[1] < n_rows else 1
            edges.append((node, node2, weight))
        if node[1] < n_rows:
            node2 = ('D', node[1] + 1, node[2])
            state2 = node2[:2]
            weight = model.transition[state][state2] if node[2] > 0 else 1 / 3
            edges.append((node, node2, weight))
        if node[2] < n and node[1] < n_rows:
            node2 = ('M', node[1] + 1, node[2] + 1)
            state2 = node2[:2]
            weight = model.transition[state][state2] * model.emission[state2][emissions[node2[2] - 1]] if node[
                                                                                                              2] > 0 else 1 / 3
            edges.append((node, node2, weight))

    # Now the edges going into the sink
    edges.extend([((state, n_rows, n), ('E', n_rows, n + 1), 1) for state in 'MDI'])

    # Add the edges to the graph
    graph.add_weighted_edges_from(edges)

    return graph, nodes


def draw_viterbi_profile_graph(graph):
    """
    Draws a Viterbi graph for a profile HMM using matplotlib, in its own window. Program execution resumes after the user has closed the window.
    :param graph: The Viterbi graph, an nx.DiGraph.
    """

    def viterbi_profile_layout(graph):
        u = 10
        v_offset = {'M': 0, 'D': u, 'I': 2 * u, 'S': 2 * u, 'E': 2 * u}
        h_offset = {'M': 0, 'D': u / 6, 'I': 2 * u / 6, 'S': u / 6, 'E': u / 6}

        n_rows = max([index for _, index, _ in graph.nodes() if index is not None])
        pos = {node: (node[2] * u + h_offset[node[0]], n_rows * 2 * u - (node[1] * 3 * u + v_offset[node[0]])) for node
               in
               graph.nodes()}
        return pos

    pos = viterbi_profile_layout(graph)
    nx.draw(graph, pos=pos, with_labels=True, font_weight='bold', node_color='orange')
    plt.show()


def align(emissions, theta, sigma, alphabet, alignment):
    """
    Returns the optimal hidden path in the HMM for a given alignment that emits a given text, and its score. The score is the natural logarithm of a probability, therefore it is a number between -infinity and 0 (-infinity not included).
    :param emissions: The text, a string.
    :param theta: The column removal threshold, a number between 0 and 1; columns in the alignment with a proportion of blanks '-' greater than or equal to the threshold are ignored.
    :param sigma: The pseudocount, a positive number; it is added to all allowed state transitions (corresponding to edges in a HMM graph) and to all emission probabilities for states that can emit (matches and insertions).
    :param alphabet: The alphabet of symbols making the alignment, a string. It must not contain the special symbol '-', reserved for blanks.
    :param alignment: The multiple alignment, a sequence of strings. Symbol '-' is used to indicate a blank.
    :return: The optimal hidden path and its score; a pair, the first element is a sequence of nodes (not including source and sink), the second a number.
    """
    # This is an HMM align model, which includes source, sink, and vertices like (<state>, <index>)
    model = make_profile_HMM(theta, sigma, alphabet, alignment)
    graph, ordered_nodes = make_profile_graph(emissions, model)

    n = len(emissions)
    n_rows = max([index for _, index, _ in graph.nodes() if index is not None])
    source = ('D', 0, 0)  # TODO replace with some sort of global
    sink = ('E', n_rows, n + 1)

    # Set the length of the longest path from source to the source (which is 0, of course)
    graph.nodes[source]['s'] = .0
    # Set the predecessor for source (which is None)
    graph.nodes[source]['previous'] = None
    # Follow the topological order
    for node in ordered_nodes:
        # For every node2 adjacent to node, determine if the longest path from source to node2 goes through node
        for node2 in graph[node].keys():
            s = graph.nodes[node]['s'] + log(graph[node][node2]['weight'])
            if s > graph.nodes[node2].get('s', float('-inf')):
                graph.nodes[node2]['s'] = s
                graph.nodes[node2]['previous'] = node

    # Backtrack from sink to source to determine the longest path found
    longest_path = []
    current_vertex = sink
    while current_vertex != source:
        previous_state = graph.nodes[current_vertex]['previous']
        longest_path.append(previous_state)
        current_vertex = previous_state
    longest_path = [(item1, item2) for item1, item2, _ in reversed(longest_path)][1:]

    return longest_path, graph.nodes[sink]['s']


def hidden_path_prob(path, transition_prob):
    prob = reduce(operator.mul, [transition_prob[i][j] for i, j in zip(path[:len(path) - 1], path[1:])], 1)
    prob *= .5
    return prob
