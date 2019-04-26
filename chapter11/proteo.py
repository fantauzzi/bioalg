import networkx as nx


def get_ammino_mass():
    """
    Returns an association between amminoacids and their masses.
    :return: A dictionary, associating the ammino acid one letter code in the key, a string, with its mass, an integer number.
    """
    ammino_mass = {
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
    return ammino_mass


def get_reverse_ammino_mass():
    """
    Returns the reverse association between amminoacids and their masses.
    :return: A dictionary, associating an amminoacid mass, an integer, with the one letter code of the amminoacid with that mass.
    """
    ammino_mass = get_ammino_mass()
    reversed = {value: key for key, value in ammino_mass.items()}
    # reversed[113] = 'I/L'
    # reversed[128] = 'K/Q'
    return reversed


def graph_from_spectrum(spectrum):
    """
    Returns the graph of a peptide spectrum.
    :param spectrum: The spectrum, a sequence of integer numbers in non-decreasing order; it may contain repeated values and it should not contain 0.
    :return: The graph of the spectrum, an nx.DiGraph with one node per spectrum value and edges labelled with attribute 'ammino', a string. The graph contains a node for mass 0.
    """
    assert spectrum == sorted(spectrum)
    spectrum = [0] + spectrum
    graph = nx.DiGraph()
    graph.add_nodes_from(spectrum)
    reverse_ammino_mass = get_reverse_ammino_mass()
    max_mass = max(reverse_ammino_mass)

    # Scan pair of nodes node1 < node2 in the graph such that node2-node1 does not exceed the maximum amminoacid mass
    for i1, node1 in enumerate(spectrum):
        for node2 in spectrum[i1 + 1:]:
            if node2 - node1 > max_mass:
                break
            ammino = reverse_ammino_mass.get(node2 - node1)
            if ammino is not None:
                graph.add_edge(node1, node2, ammino=ammino)
    return graph


def ideal_spectrum(peptide):
    """
    Returns the ideal spectrum of a given peptide.
    :param peptide: The peptide, a sequence of strings.
    :return: The peptide ideal spectrum, a list of non-decreasing integer numbers, beginning with 0.
    """
    def measure_mass(peptide):
        ammino_mass = get_ammino_mass()
        mass = sum([ammino_mass[ammino] for ammino in peptide])
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
        peptide = [graph.edges[edge]['ammino'] for edge in path]
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
    ammino_mass = get_ammino_mass()
    peptide_vector = []
    for ammino in peptide:
        mass = ammino_mass[ammino]
        peptide_vector.extend([0]*(mass-1)+[1])
    return peptide_vector


def peptide_from_vector(vector):
    """
    Returns the peptide for a given peptide vector.
    :param vector: The peptide vector, a sequence of numbers, each number being either 0 or 1.
    :return: The peptide for the given vector, a string, if it exists; None otherwise.
    """
    reverse_ammino_mass = get_reverse_ammino_mass()
    peptide = []
    previous_mass = 0
    for i in range(0, len(vector)):
        if vector[i]:
            ammino = reverse_ammino_mass.get(i+1-previous_mass)
            if ammino is None:
                return None
            peptide.append(ammino)
            previous_mass = i+1
    return ''.join(peptide)
