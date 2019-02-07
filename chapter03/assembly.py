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
    genome = path[0]+''.join(nucleotides)
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
        candidates_next=[postfix+'A', postfix+'C', postfix+'G', postfix+'T']
        for candidate in candidates_next:
            if graph_adj.get(candidate) is not None:
                assert candidate not in graph_adj[candidate]
                graph_adj[kmer].append(candidate)

    return graph_adj, kmers_count


def print_overlap_graph(graph_adj, kmers_count):
    for kmer, adj in graph_adj.items():
        if len(adj) == 0:
            continue
        print(kmer, sep='', end='')
        for adj_kmer in adj:
            adj_kmer_count = kmers_count[adj_kmer]
            for i in range(0, adj_kmer_count):
                print(' -> ', adj_kmer, sep='', end='')
        print('')


def main():
    kmers = []
    try:
        while True:
            item = input()
            kmers.append(item)
    except EOFError:
        pass
    adj, count = overlap_graph(kmers)
    print_overlap_graph(adj, count)



if __name__ == '__main__':
    main()
