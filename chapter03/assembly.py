from copy import deepcopy

""" 
NOTE
====
    
The code uses adjacency lists to represent directed graphs. Each adjacency list is a Python dictionary. Keys in the dictionary are vertices, and values are lists of vertices, that are adjacent to the vertex in the key. Every vertex is named with a string. There can be multiple edges, even with the same orientation, between two given vertices. If a vertex has no outgoing edges, but only incoming edges, then it does not appear among the keys.
"""


def composition(k, text):
    """
    Returns the list of k-mers from a text, in lexicographic order.
    :param k: The k-mer size (number of characters).
    :param text: The text, a string.
    :return: A list of strings, the k-mers in lexicographic order.
    """

    res = []
    for i in range(0, len(text) - k + 1):
        res.append(text[i:i + k])
    return sorted(res)


def path_to_genome(path):
    """
    Returns the string corresponding to the given genome path.
    :param path:The genome path, a sequence of k-mers (strings).
    :return:The string corresponding to the given genome path.
    """

    genome = path[0] + ''.join(vertex[-1] for vertex in path[1: len(path)])
    return genome


def gapped_path_to_genome(d, path):
    """
    Returns the path corresponding to a given sequence of gapped (k, d)-mers, provided such a path exists.
    :param d: The length of the gap within each gapped (k, d)-mer.
    :param path: The sequence of (k-d)-mers, each item is a pair of strings, representing two gapped k-mers.
    :return: The string of the genome path spelled by the given (k, d)-mers, if it exists, None otherwise.
    """

    pre_path, post_path = [item[0] for item in path], [item[1] for item in path]
    pre_genome = path_to_genome(pre_path)
    post_genome = path_to_genome(post_path)
    k = len(pre_path[0])
    assert len(pre_genome) == len(post_genome)
    return pre_genome[0:k + d] + post_genome if pre_genome[k + d:] == post_genome[:-k - d] else None


def print_graph(adj):
    """
    Pretty print for a graph.
    :param adj: The adjacency lists of the graph, in a dictionary.
    """
    for vertex, adjs in sorted(adj.items()):
        print(vertex, sep='', end='')
        separator = ' -> '
        for item in sorted(adjs):
            print(separator, item, sep='', end='')
            separator = ','
        print()


def de_brujin_graph(k, text):
    """
    Returns the De Brujin graph for a given text and k-mer length.
    :param k: The k-mer length.
    :param text: The text.
    :return: The adjacency lists of the De Brujin graph, a dictionary.
    """
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
    """
    Returns the Eulerian cycle for a given graph, assuming the graph contains such a cycle. It implements the Hierholzer algorithm. If the graph contains multiple Eulerian cycles, only one of them is returned.
    :param adj: The graph adjacency lists, a dictionary.
    :return: The list of vertices in the Eulerian cycle; it begins and ends with the same vertex.
    """
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
    """
    Returns a De Brujin graph from a sequence of k-mers.
    :param kmers: The sequence of k-mers.
    :return: The adjacency lists for the corresponding De Brujin graph, a dictionary.
    """
    adj = {}
    for kmer in kmers:
        prefix = kmer[0:len(kmer) - 1]
        if adj.get(prefix) is None:
            adj[prefix] = [kmer[1:]]
        else:
            adj[prefix].append(kmer[1:])
    return adj


def de_brujin_graph_from_paired_reads(reads):
    """
    Returns a De Brujin graph from a sequence of read pairs. Each read-pair is a (k-d) mer. Note that the result doesn't depend on the value for d, which therefore is not required by the function.
    :param reads: The sequence of read pairs; each element in the sequence is a (k-d)-mer expressed as a pair of strings.
    :return: The adjacency lists of the resulting graph, a dictionary.
    """
    adj = {}
    k = len(reads[0][0])
    for read in reads:
        prefix = (read[0][:k - 1], read[1][:k - 1])
        postfix = (read[0][1:], read[1][1:])
        if adj.get(prefix) is None:
            adj[prefix] = [postfix]
        else:
            adj[prefix].append(postfix)
    return adj


def parse_graph(text):
    """
    Returns the oriented graph described by a given text. Useful to feed input from Stepik challenges.
    :param text: The text, a sequence of lines; each line is an adjacency list like "Vertex1 -> Vertex2, Vertex3"
    :return:The grap adjacency lists, a dictionary.
    """
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
    """
    Prints a graph cycle or path in a human readable way, e.g. "Vertex0 -> Vertex1 -> Vertex0.
    :param cycle: The sequence of vertices in the cycle.
    """
    print(*cycle, sep=' -> ')


def is_eulerian_cycle(adj, cycle):
    """
    Returns True if and only if a given sequence of vertices is a Eulerian cycle in a graph.
    :param adj: The graph adjacency lists, a dictionary.
    :param cycle: The sequence of vertices from the graph.
    :return: True or False.
    """
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
    """
    Returns the Eulerian path for a given graph, assuming the graph contains such a path. The implementation augments the graph, turning the Eulerian path into a Eulerian cycle, finds the cycle via the Hierholzer algorithm, then converts it to a path before returning it.
    :param adj: The graph adjacency lists, a dictionary.
    :return: The list of vertices in the Eulerian path.
    """

    ''' Find the two unbalanced vertices in the graph, i.e. vertices that each have a different number of outgoing and incoming edges. They will be used to augment the graph'''

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
    """
    Returns whether a given binary string is k-universal.
    :param k: the k parameter, a positive integer.
    :param string: the binary string to be tested, a string.
    :return: True if the string is k-universal, False otherwise.
    """

    def binary_string_to_int(string):
        total = 0
        for pos, char in enumerate(reversed(string)):
            assert char in ('0', '1')
            total += (char == '1') * (2 ** pos)
        return total

    found = [False] * (2 ** k)
    extended_string = string + string[0:k - 1]
    for i in range(0, len(string)):
        kmer = extended_string[i:i + k]
        kmer_as_int = binary_string_to_int(kmer)
        if found[kmer_as_int]:
            return False
        found[kmer_as_int] = True

    res = sum(found) == len(found)
    return res


def make_k_universal_string(k):
    """
    Returns a binary k-universal circular string.
    :param k: Value for parameter k, a positive integer.
    :return: The circular string, a Python string.
    """

    # Generate all binary k-mers.
    kmers = [bin(number)[2:].rjust(k, '0') for number in range(0, 2 ** k)]

    # Build the De Brujin graph for the generated k-mers.
    adj = de_brujin_graph_from_kmers(kmers)

    # Find a Eulerian cycle in the De Brujin graph.
    cycle = eulerian_cycle(adj)

    # Spell the string along the cycle.
    res = ''.join(vertex[0] for vertex in cycle[: len(cycle) - 1])

    return res


def reconstruct_string_from_kmers(kmers):
    """
    Returns the string spelled by given k-mers.
    :param kmers: The k-mers, a sequence of strings.
    :return: The string.
    """
    k = len(kmers[0])

    # Build the De Brujin graph for the given k-mers.

    adj = de_brujin_graph_from_kmers(kmers)

    # Find a Eulerian path in the De Brujin graph.

    path = eulerian_path(adj)

    # Spell the string along the cycle.

    res = path_to_genome(path)

    return res


def reconstruct_string_from_paired_reads(d, reads):
    """
    Returns the string spelled by given paired reads, i.e. (k, d)-mers.
    :param d: The gap between reads in a pair.
    :param reads: A sequence of reads, each is a pair of (k-d)-mers (a pair of strings).
    :return: The string.
    """
    adj = de_brujin_graph_from_paired_reads(reads)

    path = eulerian_path(adj)

    ''' Need to pass d+1 as the gap length because the path is a sequence of vertices (through the De Brujin graph), which are (k-1)-mers, where the reads are k-mers. '''
    gen = gapped_path_to_genome(d + 1, path)

    return gen


def max_no_branch_paths(adj):
    """
    Returns all maximal non branching paths in a graph.
    :param adj: The graph adjacency lists, a dictionary.
    :return: A list of maximal non branching paths; each item in the list is a list of vertices (strings).
    """

    def compute_in_out(adj):
        fan_in = {}
        fan_out = {}
        for vertex, adjs in adj.items():
            # count_out = fan_out.get(vertex, 0)
            fan_out[vertex] = len(adj[vertex])
            fan_in.setdefault(vertex, 0)
            for vertex2 in adjs:
                count_in = fan_in.get(vertex2, 0)
                fan_in[vertex2] = count_in + 1
                fan_out.setdefault(vertex2, 0)
        return fan_in, fan_out

    # List of paths, to be built and returned
    paths = []
    # The in and out degree of every vertex in the graph
    fan_in, fan_out = compute_in_out(adj)
    ''' Set of vertices reachable through any of the branching paths; it is initialised to all vertices in the graph, then vertices are removed from it as discovered paths go across them. Vertices still in the set thereafter must be in non-connected cycles, and will be processed separately '''
    unprocessed_vertices = set(adj)
    # Process every vertex in the graph
    for vertex, adjs in adj.items():
        # If the vertex is a 1-in-1-out or has no outgoing egdes, then skip it
        fan_in_vertex = fan_in.get(vertex, 0)
        fan_out_vertex = fan_out.get(vertex, 0)
        if (fan_in_vertex == fan_out_vertex == 1) or fan_out_vertex == 0:
            continue
        # Build all the maximal non branching paths starting from the given vertex, and append them to paths
        for vertex2 in adjs:
            new_path = [vertex, vertex2]
            unprocessed_vertices.discard(vertex)
            unprocessed_vertices.discard(vertex2)
            # Follow each maximal non branching path until its end
            while fan_in[vertex2] == 1 and fan_out[vertex2] == 1:
                next_vertex = adj[vertex2][0]
                new_path.append(next_vertex)
                unprocessed_vertices.discard(next_vertex)
                vertex2 = next_vertex
            paths.append(new_path)

    # Now process all vertices that are in non-connected cycles
    while unprocessed_vertices:
        vertex = unprocessed_vertices.pop()
        assert fan_in[vertex] == fan_out[vertex] == 1
        next_vertex = adj[vertex][0]
        new_path = [vertex, next_vertex]
        unprocessed_vertices.discard(vertex)
        unprocessed_vertices.discard(next_vertex)
        # Follow each cycle
        while True:
            assert fan_in[next_vertex] == fan_out[next_vertex] == 1
            next_vertex = adj[next_vertex][0]
            new_path.append(next_vertex)
            # Make sure you don't go around the same cycle more than once
            if next_vertex not in unprocessed_vertices:
                break
            unprocessed_vertices.discard(next_vertex)
        paths.append(new_path)

    return paths


def contigs_from_kmers(kmers):
    """
    Returns all the contigs from the given kmers.
    :param kmers: A sequence of kmers.
    :return: A list of contigs; each element in the list is a string.
    """
    adj = de_brujin_graph_from_kmers(kmers)
    paths = max_no_branch_paths(adj)
    contigs = [path_to_genome(path) for path in paths]
    return contigs
