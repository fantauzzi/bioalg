from copy import deepcopy


def composition(k, text):
    """
    Produces the list of k-mers in a text, in lexicographic order.
    :param k: The k-mer size (number of nucleotides).
    :param text: The text, a string.
    :return: A list of strings, the k-mers in lexicographic order.
    """

    res = []
    for i in range(0, len(text) - k + 1):
        res.append(text[i:i + k])
    return sorted(res)


def path_to_genome(path):
    """
    Reconstructs a string from its genome path.
    :param path:The genome path, a sequence of k-mers (strings).
    :return:The string corresponding to the given genome path.
    """

    nucleotides = [item[-1] for item in path[1:]]
    genome = path[0] + ''.join(nucleotides)
    return genome


def overlap_graph(kmers):
    """
    Returns the overlap graph for a sequence of k-mers.
    :param kmers: The k-mers sequence.
    :return: The overlap graph, as a dictionary associating every k-mer (key in the dictionary, vertex in the graph) with its adjacency list.
    """

    # Keep count of how many times any given k-mer was in the input list of k-mers
    kmers_count = {}

    # Adjacency list for the overlap graph. Note that the same k-mer may recur multiple times in the input list of k-mers, but it appears as a key only once.
    graph_adj = {}

    for kmer in kmers:
        if kmers_count.get(kmer) is None:
            kmers_count[kmer] = 1
            graph_adj[kmer] = []
        else:
            kmers_count[kmer] += 1

    for kmer in graph_adj.keys():
        postfix = kmer[1:]
        candidates_next = [postfix + 'A', postfix + 'C', postfix + 'G', postfix + 'T']
        for candidate in candidates_next:
            if graph_adj.get(candidate) is not None:
                assert candidate not in graph_adj[kmer]
                graph_adj[kmer].append(candidate)

    return graph_adj, kmers_count


def print_overlap_graph(graph_adj, kmers_count):  # TODO who is calling this?
    """
    Pretty print for an overlap graph.
    :param graph_adj: The adjacency lists for the overlap graph, as a dictionary.
    :param kmers_count: ???
    """
    for kmer, adj in graph_adj.items():
        if len(adj) == 0:
            continue
        for _ in range(0, kmers_count[kmer]):
            print(kmer, ' -> ', sep='', end='')
            for i, adj_kmer in enumerate(adj):
                adj_kmer_count = kmers_count[adj_kmer]
                for j in range(0, adj_kmer_count):
                    print(adj_kmer, sep='', end='')
                    if j < adj_kmer_count - 1 or i < len(adj) - 1:
                        print(',', sep='', end='')
            print('')


def print_graph(adj):
    '''
    Pretty print for a graph.
    :param adj: The adjacency lists of the graph, in a dictionary.
    '''
    for vertex, adjs in sorted(adj.items()):
        print(vertex, sep='', end='')
        separator = ' -> '
        for item in sorted(adjs):
            print(separator, item, sep='', end='')
            separator = ','
        print()


def de_brujin_graph(k, text):
    '''
    Returns the De Brujin graph for a given text and k-mer length.
    :param k: The k-mer length.
    :param text: The text.
    :return: The adjacency lists of a De Brujin graph, a dictionary associating each vertex with a list of vertices.
    '''
    m = k - 1
    adj = {}
    for i in range(0, len(text) - m + 1):
        prefix = text[i:i + m]
        if i + 1 + m <= len(text):
            if adj.get(prefix) is None:
                adj[prefix] = [text[i + 1:i + 1 + m]]
            else:
                adj[prefix].append(text[i + 1:i + 1 + m])
    return adj


def eulerian_cycle(adj):
    '''
    Returns the Eulerian cycle for a given graph, assuming the graph contains such a cycle. It implements the Hierholzer algorithm.
    :param adj: The graph adjacency lists, a dictionary associating each vertex with a list of vertices.
    :return: The list of vertices in the Eulerian cycle; the list begins and ends with the same vertex.
    '''
    untraversed_adj = deepcopy(adj)

    # Choose any starting vertex
    current = next(iter(untraversed_adj.keys()))
    # Initialise the current path (used as a stack)
    curr_path = [current]
    # Initialise the cycle, to be built
    cycle = []

    while len(curr_path) > 0:
        succ_list = untraversed_adj[current]
        # If there are un-traversed edges outgoing from the current vertex...
        if len(succ_list) > 0:
            # ... then choose one and traverse it...
            curr_path.append(current)
            current = succ_list.pop()
        else:
            # ... otherwise add the current vertex to the cycle being built, and then backtrack along the current path
            cycle.append(current)
            current = curr_path.pop()

    res = list(reversed(cycle))
    return res


def de_brujin_graph_from_kmers(kmers):
    '''
    Returns a De Brujin graph from a sequence of k-mers.
    :param kmers: The sequence of k-mers.
    :return: The adjacency lists for the corresponding De Brujin graph, a dictionary.
    '''
    adj = {}
    for kmer in kmers:
        prefix = kmer[0:len(kmer) - 1]
        if adj.get(prefix) is None:
            adj[prefix] = [kmer[1:]]
        else:
            adj[prefix].append(kmer[1:])
    return adj


def parse_graph(text):
    '''
    Returns the oriented graph described by a given text. Useful to feed input from Stepik challenges.
    :param text: The text,
    :return:The grap adjacency lists, a dictionary.
    '''
    adj = {}
    for line in text:
        vertex, adj_list = str.split(line, ' -> ')
        vertex = str(int(vertex))
        adjs = str.split(adj_list, ',')
        adjs_corrected = []
        for item in adjs:
            try:
                item_corrected = str(int(item))
                adjs_corrected.append(item_corrected)
            except ValueError:
                adjs_corrected = None
                break
        if adjs_corrected is not None:
            adj[vertex] = adjs_corrected
    return adj


def print_cycle(cycle):
    '''
    Prints a graph cycle in a human readable way, e.g. Vertex0->Vertex1->Vertex0.
    :param cycle: The sequence of vertices in the cycle.
    '''
    for i, vertex in enumerate(cycle):
        if i > 0:
            print('->', sep='', end='')
        print(vertex, sep='', end='')


def is_eulerian_cycle(adj, cycle):
    '''
    Returns True if and only if a sequence of vertices is a Eulerian cycle in a graph.
    :param adj: The graph adjacency lists, a dictionary.
    :param cycle: The sequence of vertices from the graph.
    :return: True or False.
    '''
    for i in range(0, len(cycle) - 1):
        vertex1 = cycle[i]
        vertex2 = cycle[i + 1]
        adjs = adj.get(vertex1)
        if adjs is None:
            return False
        try:
            adjs.remove(vertex2)
        except ValueError:
            return False

    for _, adjs in adj.items():
        if len(adjs) > 0:
            return False

    return True


def eulerian_path(adj):
    '''
    Returns the Eulerian path for a given graph, assuming the graph contains such a path. It implements the Hierholzer algorithm, after augmenting the graph to turn the wanted path into a Eulerian cycle.
    :param adj: The graph adjacency lists, a dictionary associating each vertex with a list of vertices.
    :return: The list of vertices in the Eulerian path.
    '''

    ''' Find the two unbalanced vertices in the graph, i.e. vertices that each have a different number of outgoing and incoming edges. '''

    ''' First count and store, for every vertex, its number of outgoing and incoming edges. '''
    outbound_count = {}
    inbound_count = {}
    for vertex, adjs in adj.items():
        outbound_count[vertex] = len(adjs)
        if inbound_count.get(vertex) is None:
            inbound_count[vertex] = 0
        for adj_vertex in adjs:
            count = inbound_count.get(adj_vertex, 0)
            inbound_count[adj_vertex] = count + 1
            if outbound_count.get(adj_vertex) is None:
                outbound_count[adj_vertex] = 0

    ''' Next find the two vertices that have one more outgoing edge than incoming edges (pos_vertex),
    and one more incoming edge than outgoing edges (neg_vertex). '''
    pos_vertex, neg_vertex = None, None
    for vertex, outbound in outbound_count.items():
        inbound = inbound_count[vertex]
        if outbound > inbound:
            assert pos_vertex is None
            pos_vertex = vertex
            pos_amount = outbound - inbound
        elif outbound < inbound:
            assert neg_vertex is None
            neg_vertex = vertex
            neg_amount = outbound - inbound

    assert pos_vertex is not None and pos_amount == 1
    assert neg_vertex is not None and neg_amount == -1

    ''' Build an augmented graph, by adding one edge to the given graph, going from neg_vertex to pos_vertex. '''
    augmented_adj = deepcopy(adj)
    if augmented_adj.get(neg_vertex) is None:
        augmented_adj[neg_vertex] = [pos_vertex]
    else:
        augmented_adj[neg_vertex].append(pos_vertex)

    ''' Find the Eulerian cycle in the augmented graph. '''
    cycle = eulerian_cycle(augmented_adj)

    ''' Convert the Eulerian cycle to a path in the given graph. '''
    cycle.pop()
    cycle_cut_position = None
    for i in range(0, len(cycle)):
        if cycle[i] == neg_vertex and cycle[(i + 1) % len(cycle)] == pos_vertex:
            cycle_cut_position = (i + 1) % len(cycle)
            break
    assert cycle_cut_position is not None
    path = cycle[cycle_cut_position:] + cycle[:cycle_cut_position]

    return path


def is_k_universal(k, string):
    '''
    Returns whether a given binary string is k-universal.
    :param k: the k parameter, a positive integer.
    :param string: the binary string to be tested, a string.
    :return: True if the string is k-universal, False otherwise.
    '''

    def binary_kmer_to_int(kmer):
        total = 0
        for pos, char in enumerate(reversed(kmer)):
            assert char in ('0', '1')
            total += (char == '1') * (2 ** pos)
        return total

    found = [False] * (2 ** k)
    extended_string = string + string[0:k-1]
    for i in range(0, len(string)):
        kmer = extended_string[i:i + k]
        kmer_as_int = binary_kmer_to_int(kmer)
        if found[kmer_as_int]:
            return False
        found[kmer_as_int] = True

    res = sum(found) == len(found)
    return res


def make_k_universal(k):
    '''
    Returns a binary k-universal circular string.
    :param k: Value for parameter k.
    :return: The circular string, a Python string.
    '''

    # Generate all binary k-mers.
    kmers = [bin(number)[2:].rjust(k, '0') for number in range(0, 2**k)]

    # Build the De Brujin graph for the generated k-mers.

    adj = de_brujin_graph_from_kmers(kmers)

    # Find a Eulerian cycle in the De Brujin graph.

    cycle = eulerian_cycle(adj)

    # Spell the string along the cycle.

    res = ''.join(vertex[0] for vertex in cycle[0: len(cycle)-1])

    return res



def main_eulerian_cycle():
    """
    k = int(input())
    text = input()
    adj = de_brujin_graph(k, text)
    print_graph(adj)
    """

    text = []
    try:
        while True:
            item = input()
            text.append(item)
    except EOFError:
        pass
    adj = parse_graph(text)
    cycle = eulerian_cycle(adj)
    print_cycle(cycle)


def main_eulerian_path():
    """
    k = int(input())
    text = input()
    adj = de_brujin_graph(k, text)
    print_graph(adj)
    """

    text = []
    try:
        while True:
            item = input()
            text.append(item)
    except EOFError:
        pass
    adj = parse_graph(text)
    path = eulerian_path(adj)
    print_cycle(path)


def main_k_universal():
    k = int(input())
    res = make_k_universal(k)
    print(res)



if __name__ == '__main__':
    main_k_universal()
