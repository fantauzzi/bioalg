from copy import deepcopy
from itertools import chain, islice


def pretty_print(printme):
    """
    Prints the result of greedy_sorting() in a format accepted by the Stepik challenge.
    :param printme: The list of lists of numbers to be printed.
    """
    for line in printme:
        fmat = '{:+} ' * len(line)
        fmat = fmat.rstrip(' ')
        print(fmat.format(*line))


def pretty_print_2_break_seq(printme):
    for line in printme:
        for permutation in line:
            print('(', sep='', end='')
            fmat = '{:+} ' * len(permutation)
            fmat = fmat.rstrip(' ')
            print(fmat.format(*permutation), end='')
            print(')', sep='', end='')
        print()


def pretty_print_seq_of_seqs(printme):
    for item in printme:
        fmat = '{}, ' * len(item)
        fmat = fmat.rstrip(', ')
        print('(', fmat.format(*item), ')', sep='')


def greedy_sorting(p):
    """
    Returns the sequence of permutations leading from a given permutation to the indentity permutation, implementing greedy sorting.
    :param p: The given permutation, a list of integer numbers optionally with sign; positions in the permutation must be numbered starting from 1 (not 0).
    :return: A list of permutations that, applying greedy sorting, lead from the given one to the identity one, a list of lists of integer numbers. The list does not contain the given permutation, but ends with the identity permutation.
    """
    steps = []
    for i in range(0, len(p)):
        if abs(p[i]) != i + 1:
            try:
                pos = p.index(i + 1, i + 1)
            except ValueError:
                pos = p.index(-(i + 1), i + 1)
            p = p[:i] + [-item for item in p[i:pos + 1][::-1]] + p[pos + 1:]
            steps.append(deepcopy(p))
            assert abs(p[i]) == i + 1
        if p[i] == - (i + 1):
            p[i] = i + 1
            steps.append(deepcopy(p))

    return steps


def count_breakpoints(p):
    """
    Returns the number of breakpoints in a given permutation.
    :param p: The permutation, a list of integer numbers.
    :return: The number of breakpoints, an integer.
    """
    p = [0] + p + [len(p) + 1]
    count = sum([p[i + 1] - p[i] != 1 for i in range(0, len(p) - 1)])
    return count


def chromosome_to_cycle(chromosome):
    def flatten(seq_of_seq):
        return list(chain.from_iterable(seq_of_seq))

    cycle = flatten([(2 * item - 1, 2 * item) if item > 0 else (-2 * item, -2 * item - 1) for item in chromosome])
    return cycle


def pairwise(seq):
    assert len(seq) % 2 == 0
    return zip(islice(seq, 0, None, 2), islice(seq, 1, None, 2))


def cycle_to_chromosome(cycle):
    chromosome = [item2 // 2 if item1 < item2 else -item1 // 2 for item1, item2 in pairwise(cycle)]
    return chromosome


def colored_edges(genome):
    cycles = [chromosome_to_cycle(chromosome) for chromosome in genome]
    cycles = [cycle + [cycle[0]] for cycle in cycles]
    edges = [(block1, block2) for cycle in cycles for (block1, block2) in pairwise(cycle[1:])]
    return edges


def graph_to_genome(edges):
    def get_block(vertex):
        block = vertex // 2 if vertex % 2 == 0 else (vertex + 1) // 2
        return block

    genome = []
    cycle = []
    cycle_start = None
    for edge in edges:
        if cycle_start is None:
            cycle_start = edge[0]
        cycle.extend(edge)
        block_0, block_1 = get_block(edge[0]), get_block(edge[1])
        is_trivial_cycle = block_0 == block_1
        if ((edge[1] + 1 == cycle_start or edge[1] - 1 == cycle_start) and edge[0] != cycle_start) or is_trivial_cycle:
            cycle = [cycle[-1]] + cycle[0:len(cycle) - 1]
            genome.append(cycle_to_chromosome(cycle))
            cycle = []
            cycle_start = None
    return genome


def add_unidirectional_edge(adj, vertex1, vertex2, color):
    adjs = adj.get(vertex1, [])
    adjs.append((vertex2, color))
    adj[vertex1] = adjs


def add_bidirectional_edge(adj, vertex1, vertex2, color):
    assert vertex1 != vertex2
    add_unidirectional_edge(adj, vertex1, vertex2, color)
    add_unidirectional_edge(adj, vertex2, vertex1, color)


def find_adj_vertex_with_color(adj, vertex, color):
    for vertex2, color2 in adj[vertex]:
        if color2 == color:
            return vertex2
    assert False


def remove_unidirectional_edge(adj, vertex1, vertex2, color):
    adjs = adj[vertex1]
    new_adjs = []
    for adj_vertex, adj_color in adjs:
        if not (adj_vertex == vertex2 and adj_color == color):
            new_adjs.append((adj_vertex, adj_color))
    adj[vertex1] = new_adjs


def remove_bidirectional_edge(adj, vertex1, vertex2, color):
    assert vertex1 != vertex2
    remove_unidirectional_edge(adj, vertex1, vertex2, color)
    remove_unidirectional_edge(adj, vertex2, vertex1, color)


def DNA_complement(dna):
    complements = {'A': 'T',
                   'T': 'A',
                   'G': 'C',
                   'C': 'G'}
    complement = [complements[nucleotide] for nucleotide in reversed(dna)]
    as_string = ''.join(complement)
    return as_string


def find_shared_kmers(k, string1, string2):
    def add_position(kmers_position, kmer, pos):
        positions = kmers_position.get(kmer, [])
        positions.append(pos)
        kmers_position[kmer] = positions

    kmers_position = {}
    for i in range(0, len(string1) - k + 1):
        kmer = string1[i: i + k]
        add_position(kmers_position, kmer, i)
        kmer_rc = DNA_complement(kmer)
        if kmer != kmer_rc:
            add_position(kmers_position, kmer_rc, i)

    res = []
    for i in range(0, len(string2) - k + 1):
        kmer = string2[i: i + k]
        positions_in_1 = kmers_position.get(kmer)
        if positions_in_1 is None:
            continue
        for pos in positions_in_1:
            res.append((pos, i))

    return res


def fix_graph_cycles(graph):
    # Build the graph adjacency lists, colore edges first, black edges next
    adj = {}
    for edge in graph:
        add_bidirectional_edge(adj, edge[0], edge[1], 'colored')
        add_unidirectional_edge(adj, edge[0], edge[0] - 1 if edge[0] % 2 == 0 else edge[0] + 1, 'black')
        add_unidirectional_edge(adj, edge[1], edge[1] - 1 if edge[1] % 2 == 0 else edge[1] + 1, 'black')

    # Follow the cycles in the graph, and produce the colored edges in the right order and orientation
    unvisited = set(adj)

    edges = []
    # previous_edge_color = None
    cycles_count = 0
    while unvisited:
        starting_vertex = next(iter(unvisited))
        current_vertex = None
        cycles_count += 1
        previous_edge_color = None
        while current_vertex != starting_vertex:
            if current_vertex is None:
                current_vertex = starting_vertex
            unvisited.remove(current_vertex)
            current_edge_color = 'black' if previous_edge_color in ('None', 'colored') else 'colored'
            next_vertex = find_adj_vertex_with_color(adj, current_vertex, current_edge_color)
            if current_edge_color == 'colored':
                edges.append((current_vertex, next_vertex))
            current_vertex = next_vertex
            previous_edge_color = current_edge_color

    return edges, cycles_count


def two_break_on_genome_graph(edges, i1, i2, i3, i4):
    new_graph = []
    for edge in edges:
        if edge == (i3, i4):
            new_graph.append((i3, i1))
        elif edge == (i4, i3):
            new_graph.append((i1, i3))
        elif edge == (i1, i2):
            new_graph.append((i4, i2))
        elif edge == (i2, i1):
            new_graph.append((i2, i4))
        else:
            new_graph.append(edge)

    return new_graph


def two_break_on_genome(genome, i1, i2, i3, i4):
    edges = colored_edges(genome)
    edges = two_break_on_genome_graph(edges, i1, i2, i3, i4)
    edges, cycles = fix_graph_cycles(edges)
    genome = graph_to_genome(edges)
    return genome


def add_colored_edges_to_graph(adj, edges, color):
    for edge in edges:
        add_unidirectional_edge(adj, edge[0], edge[1], color)
        add_unidirectional_edge(adj, edge[1], edge[0], color)


def two_break_distance(genome1, genome2):
    red_edges = colored_edges(genome1)
    blue_edges = colored_edges(genome2)
    adj = {}
    add_colored_edges_to_graph(adj, red_edges, 'red')
    add_colored_edges_to_graph(adj, blue_edges, 'blue')
    # Follow the cycles in the graph, and produce the colored edges in the right order and orientation
    unvisited = set(adj)

    previous_edge_color = None
    cycles_count = 0
    while unvisited:
        starting_vertex = next(iter(unvisited))
        current_vertex = None
        cycles_count += 1
        while current_vertex != starting_vertex:
            if current_vertex is None:
                current_vertex = starting_vertex
            unvisited.remove(current_vertex)
            current_edge_color = 'red' if previous_edge_color in ('None', 'blue') else 'blue'
            next_vertex = find_adj_vertex_with_color(adj, current_vertex, current_edge_color)
            current_vertex = next_vertex
            previous_edge_color = current_edge_color

    assert len(adj) % 2 == 0
    dist = len(adj) // 2 - cycles_count
    return dist


def two_break_sorting(genome1, genome2):
    sequence = [deepcopy(genome1)]
    current_genome = genome1

    red_edges = colored_edges(genome1)
    blue_edges = colored_edges(genome2)

    # Fill in the adjacency lists of the breakpoint graph for genome1 and genome2
    adj = {}
    add_colored_edges_to_graph(adj, red_edges, 'red')
    add_colored_edges_to_graph(adj, blue_edges, 'blue')

    has_non_trivial_cycles = True
    # Iterate as long as the breakpoint graph contains non-trivial (red-blue) cycle(s), but anyway at least once
    while has_non_trivial_cycles:
        ''' Look for a non-trivial cycle, exploring cycles one at a time, until you have explored them all, or have
        found a non-trivial one '''
        has_non_trivial_cycles = False
        unvisited = set(adj)
        while unvisited and not has_non_trivial_cycles:
            vertex1 = next(iter(unvisited))
            unvisited.remove(vertex1)
            vertex2 = find_adj_vertex_with_color(adj, vertex1, 'red')
            unvisited.remove(vertex2)
            vertex3 = find_adj_vertex_with_color(adj, vertex2, 'blue')
            # If it is a trivial red-blue cycle, proceed looking for another cycle
            if vertex3 == vertex1:
                continue
            unvisited.remove(vertex3)
            has_non_trivial_cycles = True
            vertex4 = find_adj_vertex_with_color(adj, vertex3, 'red')
            unvisited.remove(vertex4)
        if has_non_trivial_cycles:
            # Update the breakpoint graph, based on the red-blue-red edges at the beginning of the non-trivial cycle just found
            remove_bidirectional_edge(adj, vertex1, vertex2, 'red')
            remove_bidirectional_edge(adj, vertex3, vertex4, 'red')
            add_bidirectional_edge(adj, vertex1, vertex4, 'red')
            add_bidirectional_edge(adj, vertex2, vertex3, 'red')
            # Do a 2-break on current_genome, and use the result to update the sequence of steps to be returned at the end
            current_genome = two_break_on_genome(current_genome, vertex1, vertex2, vertex4, vertex3)
            sequence.append(current_genome)

    return sequence
