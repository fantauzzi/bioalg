from copy import deepcopy


def pretty_print(printme):
    """
    Prints the result of greedy_sorting() in a format accepted by the Stepik challenge.
    :param printme: The list of lists of numbers to be printed.
    """
    fmat = '{:+} ' * len(printme[0])
    fmat = fmat.rstrip(' ')
    for line in printme:
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


def edge_from_block(block):
    assert block != 0
    res = [2 * block - 1, 2 * block] if block > 0 else [2 * -block - 1, 2 * -block]
    if block < 0:
        res[0], res[1] = res[1], res[0]
    return res


def graph_from_permuations(ps, color, adj=None):
    def add_edge(adj, vertex1, vertex2):
        adjs = adj.get(vertex1, [])
        adjs.append((vertex2, color))
        adj[vertex1] = adjs
        adjs = adj.get(vertex2, [])
        adjs.append((vertex1, color))
        adj[vertex2] = adjs

    if adj == None:
        adj = {}
    for p in ps:
        for i in range(0, len(p)):
            block = p[i]
            next_block = p[(i + 1) % len(p)]
            _, v1 = edge_from_block(block)
            v2, _ = edge_from_block(next_block)
            add_edge(adj, v1, v2)
    return adj


def breakpoint_graph(ps, qs, color_ps, color_qs):
    adj = graph_from_permuations(ps, color_ps)
    adj = graph_from_permuations(qs, color_qs, adj=adj)
    return adj


def two_break_dist(ps, qs):
    adj = breakpoint_graph(ps, qs, 'red', 'blue')

    # Visit all vertices along the breakpoint graph, counting cycles

    # Initially, all vertices are (yet) unvisited
    unvisited = set(adj)
    cycles_count = 0
    while unvisited:
        # Begin to follow a cycle, starting from any unvisited vertex; remove it from the set of unvisited vertices.
        starting_vertex = unvisited.pop()
        # Update the cycles counter
        cycles_count += 1
        vertex = None
        prev_color = None
        # Go from vertex to vertex in the graph, until you return to the starting vertex
        while vertex != starting_vertex:
            if vertex is None:
                vertex = starting_vertex
            # Choose any adjacent vertex that hasn't been visited yet, and connected by an edge of alternating color
            for next_vertex, next_color in adj[vertex]:
                if next_color != prev_color and (next_vertex in unvisited or next_vertex == starting_vertex):
                    break
            else:
                assert False
            ''' Remove the chosen adjacent vertex from the visited ones (unless it is the starting vertex,
            closing the cycle, which has been removed from the set already)'''
            if next_vertex != starting_vertex:
                unvisited.remove(next_vertex)
            # Move on to the next vertex, and keep following the cycle from there
            prev_color = next_color
            vertex = next_vertex

    # The distance is the number of synteny blocks minus the number of cycles
    dist = len(adj) / 2 - cycles_count

    return dist
