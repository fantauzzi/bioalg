from copy import deepcopy


def composition(k, text):
    """
    Produces the list of k-mers in a string, in lexicographic order.
    :param k: The k-mer size (number of nucleotides).
    :param text: The string.
    :return: A list of strings, the k-mers in lexicographic order.
    """

    res = []
    for i in range(0, len(text) - k + 1):
        res.append(text[i:i + k])
    return sorted(res)


def path_to_genome(path):
    """
    Reconstruct a string from its genome path.
    :param path:The genome path, a sequence of k-mers.
    :return:The string corresponding to the given genome path.
    """

    nucleotides = [item[-1] for item in path[1:]]
    genome = path[0] + ''.join(nucleotides)
    return genome


def overlap_graph(kmers):
    """
    Construct the overlap graph of a collection of k-mers.
    :param kmers: The given k-mers.
    :return: The overlap graph, as a dictionary associating every k-mer (key in the dictionary) with its adjacency list.
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


def print_overlap_graph(graph_adj, kmers_count):
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
    for vertex, adjs in sorted(adj.items()):
        print(vertex, sep='', end='')
        separator = ' -> '
        for item in sorted(adjs):
            print(separator, item, sep='', end='')
            separator = ','
        print()


def de_brujin_graph(k, text):
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
    adj = {}
    for kmer in kmers:
        prefix = kmer[0:len(kmer) - 1]
        if adj.get(prefix) is None:
            adj[prefix] = [kmer[1:]]
        else:
            adj[prefix].append(kmer[1:])
    return adj


def cycle_from_adj(adj, start):
    cycle_adj = deepcopy(adj)

    vertex, adj_list = start, cycle_adj[start]
    cycle = [vertex]
    while len(adj_list) > 0:
        vertex = adj_list.pop(0)
        adj_list = cycle_adj[vertex]
        cycle.append(vertex)

    return cycle


def parse_graph(text):
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
    for i, vertex in enumerate(cycle):
        if i > 0:
            print('->', sep='', end='')
        print(vertex, sep='', end='')


def is_eulerian_cycle(adj, cycle):
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

    augmented_adj = deepcopy(adj)
    if augmented_adj.get(neg_vertex) is None:
        augmented_adj[neg_vertex] = [pos_vertex]
    else:
        augmented_adj[neg_vertex].append(pos_vertex)

    cycle = eulerian_cycle(augmented_adj)
    cycle.pop()
    cycle_cut_position = None
    for i in range(0, len(cycle)):
        if cycle[i] == neg_vertex and cycle[(i + 1) % len(cycle)] == pos_vertex:
            cycle_cut_position = (i + 1) % len(cycle)
            break
    assert cycle_cut_position is not None
    path = cycle[cycle_cut_position:] + cycle[:cycle_cut_position]
    return path


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


if __name__ == '__main__':
    main_eulerian_path()
