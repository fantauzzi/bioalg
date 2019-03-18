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


def graph_to_genome(graph):
    genome = []
    cycle = []
    cycle_start = None
    for edge in graph:
        if cycle_start is None:
            cycle_start = edge[0]
        cycle.extend(edge)
        if (edge[1] + 1 == cycle_start or edge[1] - 1 == cycle_start) and edge[0] != cycle_start:
            cycle = [cycle[-1]] + cycle[0:len(cycle) - 1]
            genome.append(cycle_to_chromosome(cycle))
            cycle = []
            cycle_start = None
    return genome


def two_break_on_genome_graph(graph, i1, i2, i3, i4):
    new_graph = []
    for edge in graph:
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
