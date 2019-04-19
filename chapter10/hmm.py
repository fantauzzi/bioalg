from collections import namedtuple
from itertools import product, chain
import networkx as nx
from math import log
from numpy import isclose

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
    def threshold_alignment(alignment, theta):
        n_cols = len(alignment[0])
        n_rows = len(alignment)
        alignment_T = list(map(list, zip(*alignment)))  # Transpose alignment
        # Count occurrences of '-' in the transposed of alignment, going by rows
        count_spaces = [sum([item == '-' for item in row]) for row in alignment_T]
        # Set of column indices in alignment (row indices in alignment_T) where the count of '-' is above threshold
        thresholded = {col for col in range(0, n_cols) if count_spaces[col] / n_rows >= theta}
        return thresholded

    def validate_prob_matrix(matrix):
        for row in matrix.values():
            row_total = sum(row.values())
            assert isclose(row_total, 1) or isclose(row_total, 0)

    def normalise_prob_matrix(matrix):
        for row_key, row in matrix.items():
            total = sum(row.values())
            matrix[row_key] = {col_key: 0 if item == 0 else item / total for col_key, item in row.items()}

    source = ('S', None)
    sink = ('E', None)
    n_rows = len(alignment)
    n_cols = len(alignment[0])
    thresholded = threshold_alignment(alignment, theta)
    n_unthresholded = n_cols - len(thresholded)
    states = [source, ('I', 0)] + \
             flatten([(('M', pos), ('D', pos), ('I', pos)) for pos in range(1, n_cols - len(thresholded) + 1)]) + \
             [sink]

    node_visits = Counter()
    edge_visits = Counter()
    node_emissions = {state: Counter() for state in states}

    # For every sequence in alignement, follow the sequence through the HMM and increment counters:
    # - For every vertex: number of visits to the vertex; for every symbol, number of symbols emitted per symbol in that vertex
    # - For every edge: number of visits to the edge
    # Beware of source and sink

    for sequence in alignment:
        prev_state = ('S', None)
        pos = 0
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
            if curr_state[0] in 'MI':
                node_emissions[curr_state][symbol] = node_emissions[curr_state][symbol] + 1
            prev_state = curr_state
        edge_visits[(prev_state, sink)] = edge_visits[(prev_state, sink)] + 1
        node_visits[source] = n_rows
        node_visits[sink] = n_rows

    transition = {state_row: {
        state_col: 0 if node_visits[state_row] == 0 else edge_visits[state_row, state_col] / node_visits[state_row]
        for state_col in states}
        for state_row in states}

    validate_prob_matrix(transition)

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

    emission = {state: {symbol: 0 if total == 0 else node_emissions[state][symbol] / total
                        for total in (sum(node_emissions[state].values()),)
                        for symbol in alphabet}
                for state in states}

    validate_prob_matrix(emission)

    if sigma:
        for row_key, row in emission.items():
            # Pseudo-counts only apply to states that can emit, that is Insertion ('I') and Match ('M')
            if row_key[0] in 'MI':
                emission[row_key] = {col_key: item + sigma for col_key, item in row.items()}
        normalise_prob_matrix(emission)

    the_hmm = HMM(alphabet=alphabet, states=states, transition=transition, emission=emission)
    return the_hmm
