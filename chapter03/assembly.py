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

    return cycle


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


def main():
    for k in (2, 3, 4):
        adj = de_brujin_graph(k, 'TAATGCCATGGGATGTT')
        print_graph(adj)
        print()


def main2():
    """
    k = int(input())
    text = input()
    adj = de_brujin_graph(k, text)
    print_graph(adj)
    """

    kmers = []
    try:
        while True:
            item = input()
            kmers.append(item)
    except EOFError:
        pass
    adj = de_brujin_graph_from_kmers(kmers)
    print_graph(adj)


if __name__ == '__main__':
    main()
