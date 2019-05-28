import networkx as nx


def get_amino_mass():
    """
    Returns an association between amminoacids and their masses.
    :return: A dictionary, associating the ammino acid one letter code in the key, a string, with its mass, an integer number.
    """
    amino_mass = {
        'G': 57,
        'A': 71,
        'S': 87,
        'P': 97,
        'V': 99,
        'T': 101,
        'C': 103,
        'I': 113,
        'L': 113,
        'N': 114,
        'D': 115,
        'K': 128,
        'Q': 128,
        'E': 129,
        'M': 131,
        'H': 137,
        'F': 147,
        'R': 156,
        'Y': 163,
        'W': 186
    }
    """amino_mass = {'X': 4,
                  'Z': 5
                  }"""
    return amino_mass


def get_reverse_amino_mass():
    """
    Returns the reverse association between amino acids and their masses.
    :return: A dictionary, associating an amino acid mass, an integer, with the one letter code of the amino acid with that mass.
    """
    amino_mass = get_amino_mass()
    reversed = {value: key for key, value in amino_mass.items()}
    # reversed[113] = 'I/L'
    # reversed[128] = 'K/Q'
    return reversed


def graph_from_spectrum(spectrum):
    """
    Returns the graph of a peptide spectrum.
    :param spectrum: The spectrum, a sequence of integer numbers in non-decreasing order; it may contain repeated values and it should not contain 0.
    :return: The graph of the spectrum, an nx.DiGraph with one node per spectrum value and edges labelled with attribute 'amino', a string. The graph contains a node for mass 0.
    """
    assert spectrum == sorted(spectrum)
    spectrum = [0] + spectrum
    graph = nx.DiGraph()
    graph.add_nodes_from(spectrum)
    reverse_amino_mass = get_reverse_amino_mass()
    max_mass = max(reverse_amino_mass)

    # Scan pair of nodes node1 < node2 in the graph such that node2-node1 does not exceed the maximum amino acid mass
    for i1, node1 in enumerate(spectrum):
        for node2 in spectrum[i1 + 1:]:
            if node2 - node1 > max_mass:
                break
            amino = reverse_amino_mass.get(node2 - node1)
            if amino is not None:
                graph.add_edge(node1, node2, amino=amino)
    return graph


def ideal_spectrum(peptide):
    """
    Returns the ideal spectrum of a given peptide.
    :param peptide: The peptide, a sequence of strings.
    :return: The peptide ideal spectrum, a list of non-decreasing integer numbers, beginning with 0.
    """

    def measure_mass(peptide):
        amino_mass = get_amino_mass()
        mass = sum([amino_mass[amino] for amino in peptide])
        return mass

    spectrum = []
    # Add to 'spectrum' the total mass of prefixes and postfixes in 'peptide'.
    for i in range(0, len(peptide) + 1):
        prefix = peptide[:i]
        postfix = peptide[i:]
        spectrum.append(measure_mass(prefix))
        # Make sure that the empty peptide and the one that matches the given peptide are processed only once
        if i != 0 and i != len(peptide):
            spectrum.append(measure_mass(postfix))
    return sorted(spectrum)


def peptide_from_ideal_spectrum(spectrum):
    """
    Returns a peptide corresponding to a given ideal spectrum.
    :param spectrum: The ideal spectrum, a sequence of integer numbers in non-decreasing order; it may contain repeated values and it should not contain 0.
    :return: The peptide decoded from the given spectrum, a string, if such a peptide exists; None otherwise.
    """
    spectrum = [0] + spectrum
    assert sorted(spectrum) == spectrum
    graph = graph_from_spectrum(spectrum)
    paths = nx.all_simple_paths(graph, 0, max(spectrum))
    for path in map(nx.utils.pairwise, paths):
        peptide = [graph.edges[edge]['amino'] for edge in path]
        ideal = ideal_spectrum(peptide)
        if ideal == spectrum:
            return ''.join(peptide)
    return None


def vector_from_peptide(peptide):
    """
    Returns the peptide vector for a given peptide.
    :param peptide: The given peptide, a string.
    :return: The peptide vector, a list of integer numbers (0 and 1).
    """
    amino_mass = get_amino_mass()
    peptide_vector = []
    for amino in peptide:
        mass = amino_mass[amino]
        peptide_vector.extend([0] * (mass - 1) + [1])
    return peptide_vector


def peptide_from_vector(vector):
    """
    Returns the peptide for a given peptide vector.
    :param vector: The peptide vector, a sequence of numbers, each number being either 0 or 1.
    :return: The peptide for the given vector, a string, if it exists; None otherwise.
    """
    reverse_amino_mass = get_reverse_amino_mass()
    peptide = []
    previous_mass = 0
    for i in range(0, len(vector)):
        if vector[i]:
            amino = reverse_amino_mass.get(i + 1 - previous_mass)
            if amino is None:
                return None
            peptide.append(amino)
            previous_mass = i + 1
    return ''.join(peptide)


def peptide_from_spectral_vector(spectrum):
    """
    Returns the peptide with maximum score against a given spectral vector.
    :param spectrum: The spectral vector, a sequence of integer numbers.
    :return: The peptide, a string.
    """
    spectrum = [0] + spectrum
    m = len(spectrum)
    # Make a DAG for 'spectrum', with nodes 0, 1, ... m
    dag = nx.DiGraph()
    dag.add_nodes_from(range(0, m))
    # Connect node i to j if j-i is the mass of any amino acid. Assign weight -spectrum[i] to all edges leaving from node i.
    # Weights are assigned with changed sign so that finding the longest path becomes findind the shortest, and we can apply Bellman-Ford
    reverse_amino_mass = get_reverse_amino_mass()
    max_mass = max(reverse_amino_mass)
    # Scan pair of nodes i < j in the graph such that j-i does not exceed the maximum amino acid mass
    for i in range(0, m):
        for j in range(i + 1, m):
            if j - i > max_mass:
                break
            amino = reverse_amino_mass.get(j - i)
            if amino is not None:
                dag.add_edge(i, j, amino=amino, weight=-spectrum[i])

    # Find the shortest path from source (0) to sink (m), which corresponds to the longest path before weights were changed of sign
    longest_path = nx.algorithms.shortest_paths.generic.shortest_path(dag,
                                                                      source=0,
                                                                      target=m - 1,
                                                                      weight='weight',
                                                                      method='bellman-ford')  # Because weights can be < 0
    path_edges = nx.utils.pairwise(longest_path)
    peptide = [dag[node1][node2]['amino'] for (node1, node2) in path_edges]
    return ''.join(peptide)


def identify_peptide_from_proteome(spectrum, proteome):
    """
    Returns a substring of a given proteome with maximum score against a spectral vector, and its score. The score is the dot product between the peptide vector of the substring and the spectrum, when they have the same length, -infinity otherwise.
    :param spectrum: The spectral vector, a sequence of 1 and 0.
    :param proteome: The given proteome, a string.
    :return: A pair with the proteome substring that maximises the score, a string, and its score, an integer.
    """

    def score_peptide(peptide, spectrum):
        peptide_v = vector_from_peptide(peptide)
        assert len(peptide_v) == len(spectrum)
        score = sum([digit * spectrum_item for digit, spectrum_item in zip(peptide_v, spectrum)])
        return score

    amino_mass = get_amino_mass()
    expected_mass = len(spectrum)
    m = len(proteome)
    best_score = float('-inf')
    best_peptide = None
    for i in range(0, m):
        mass = 0
        for j in range(i + 1, m):
            peptide = proteome[i:j]
            mass += amino_mass[proteome[j - 1]]
            if mass == expected_mass:
                score = score_peptide(peptide, spectrum)
                if score > best_score:
                    best_score = score
                    best_peptide = peptide
            elif mass > expected_mass:
                break

    return best_peptide, best_score


def psm_search(proteome, spectral_vectors, threshold):
    """
    Returns all substrings of a given proteome with the highest score against one vector in a collection of spectral vectors that meet a certain threshold.
    :param proteome: The given proteome, a string.
    :param spectral_vectors: The collection of spectral vectors, a sequence of sequences of 1s and 0s.
    :param threshold: The threshold, a number.
    :return: The highest scoring substrings whose score is at least equal to the threhsold, a list of strings.
    """
    res = []  # The resulting list of proteoms
    for spectrum in spectral_vectors:
        peptide, score = identify_peptide_from_proteome(spectrum, proteome)
        if score >= threshold:
            res.append(peptide)
    return res


def size_of_spectral_dict(spectrum, threshold, max_score):
    amino_mass = get_amino_mass()
    pre_computed_size = {}

    def size(i, t):
        if (i, t) == (0, 0):
            return 1  # Only the empty amino acid (empty string) has mass 0
        if t < 0 or i <= 0:
            return 0
        the_size = pre_computed_size.get((i, t))
        if the_size is None:
            the_size = sum([size(i - amino_mass[amino], t - spectrum[i - 1]) for amino in amino_mass.keys()])
            pre_computed_size[(i, t)] = the_size
        return the_size

    m = len(spectrum)
    totalsize = sum([size(m, t) for t in range(threshold, max_score + 1)])

    return totalsize


def prob_of_spectral_dict(spectrum, threshold, max_score):
    amino_mass = get_amino_mass()
    pre_computed_prob = {}

    def prob(i, t):
        if (i, t) == (0, 0):
            return 1  # Only the empty amino acid (empty string) has mass 0
        if t < 0 or i <= 0:
            return 0
        the_prob = pre_computed_prob.get((i, t))
        if the_prob is None:
            the_prob = sum([prob(i - amino_mass[amino], t - spectrum[i - 1]) for amino in amino_mass.keys()]) / len(
                amino_mass)
            pre_computed_prob[(i, t)] = the_prob
        return the_prob

    m = len(spectrum)
    totalprob = sum([prob(m, t) for t in range(threshold, max_score + 1)])

    return totalprob


def spectral_alignment(peptide, spectrum, k):
    amino_mass = get_amino_mass()
    m_plus_Delta = len(spectrum)

    pre_computed = {}
    diff = {}
    cumulative = 0
    for pept_i in range(1, len(peptide) + 1):
        cumulative += amino_mass[peptide[pept_i - 1]]
        diff[cumulative] = amino_mass[peptide[pept_i - 1]]

    def score(i, j, t):
        if (i, j, t) == (0, 0, 0):
            return 0, []
        if (i, j) == (0, 0) and t >= 1:
            return float('-inf'), None
        if i == 0 and j > 0 or j == 0 and i > 0:
            return float('-inf'), None
        best_score, best_path = pre_computed.get((i, j, t), (None, None))
        if best_score is None:
            best_score, best_path = score(i - diff[i], j - diff[i], t) if j - diff[i] >= 0 else (float('-inf'), None)
            for j_p in range(0, j):
                score2, path2 = score(i - diff[i], j_p, t - 1) if t >= 1 else (float('-inf'), None)
                if score2 > best_score:
                    best_score = score2
                    best_path = path2
            pre_computed[i, j, t] = spectrum[j - 1] + best_score, (best_path if best_path is None else best_path + [(i, j, t)])
            return spectrum[j - 1] + best_score, (best_path if best_path is None else best_path + [(i, j, t)])
        else:
            return best_score, best_path

    m = sum([amino_mass[amino] for amino in peptide])
    best_path_score = float('-inf')
    best_path = None
    for t in range(0, k + 1):
        the_score, path = score(m, m_plus_Delta, t)
        if the_score > best_path_score:
            best_path_score = the_score
            best_path = path

    best_path = [(0, 0, 0)] + best_path


    assert len(best_path) == len(peptide) + 1
    res = ''
    for pept_i in range(1, len(peptide) + 1):
        amino = peptide[pept_i - 1]
        i, j, _ = best_path[pept_i]
        i_prev, j_prev, _ = best_path[pept_i - 1]
        variant = j - (j_prev + (i - i_prev))
        res = (res + amino) if not variant else res + '{}({:+d})'.format(amino, variant)

    return res, best_path_score
