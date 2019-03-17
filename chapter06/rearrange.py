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


def graph_from_permutations(ps, color, adj=None):
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
    adj = graph_from_permutations(ps, color_ps)
    adj = graph_from_permutations(qs, color_qs, adj=adj)
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


def two_break(adj, v1, v2, v3, v4):
    pass


def permutations_from_breakpoint_graph(adj):
    # Visit all vertices along the breakpoint graph, counting cycles

    # Initially, all vertices are (yet) unvisited
    unvisited = set(adj)
    ps = []
    while unvisited:
        # Initially, all vertices are (yet) unvisited
        # Begin to follow a cycle, starting from any unvisited vertex; remove it from the set of unvisited vertices.
        # starting_vertex = next(iter(unvisited))
        starting_vertex = min(unvisited)
        unvisited.remove(starting_vertex)
        vertex = None
        p = []
        # Go from vertex to vertex in the graph, until you return to the starting vertex
        while vertex != starting_vertex:
            if vertex is None:
                vertex = starting_vertex
            # Choose the adjacent
            next_vertex, _ = adj[vertex][0]
            assert len(adj[vertex]) == 1
            ''' Remove the chosen adjacent vertex from the visited ones (unless it is the starting vertex,
            closing the cycle, which has been removed from the set already)'''
            if next_vertex != starting_vertex:
                p.append(- next_vertex / 2 if next_vertex % 2 == 0 else (next_vertex + 1) / 2)
                if next_vertex not in unvisited:
                    dummy = 42
                unvisited.remove(next_vertex)
            # Move on to the next vertex, and keep following the cycle from there
            vertex = next_vertex + 1 if next_vertex % 2 == 1 else next_vertex - 1
            if vertex != starting_vertex:
                unvisited.remove(vertex)

        p = [p[-1]] + p[:len(p) - 1]
        ps.append(p)
    if len(ps) == 1:
        ps = ps[0]
    return ps


def shortest_rearrangement(ps, qs):
    def find_adj_with_color(adj, vertex, color):
        for vertex2, color2 in adj[vertex]:
            if color2 == color:
                return vertex2
        else:
            assert False

    adj = breakpoint_graph(ps, qs, 'red', 'blue')
    adj_p = graph_from_permutations(ps, 'red')

    # Repeat as long as a non-trivial cycle is found (but at least once)
    found_non_trivial_cycle = None
    while found_non_trivial_cycle is None or found_non_trivial_cycle:
        # Look for a non-trivial cycle, following cycles one at a time starting from unvisited vertices
        unvisited = set(adj)  # At first, all vertices are (yet) unvisited
        found_non_trivial_cycle = False
        while unvisited:
            # Begin to follow a cycle, starting from any unvisited vertex; remove it from the set of unvisited vertices.
            vertex1 = unvisited.pop()
            # Find an adjacent vertex along a red edge
            vertex2 = find_adj_with_color(adj, vertex1, 'red')
            unvisited.remove(vertex2)
            # Find the next one along a blue edge
            vertex3 = find_adj_with_color(adj, vertex2, 'blue')
            # Check if vertex1-vertex2-vertex3 make a trivial cycle
            if vertex3 == vertex1:
                continue
            found_non_trivial_cycle = True
            unvisited.remove(vertex3)
            # Find the next one along a red edge
            vertex4 = find_adj_with_color(adj, vertex3, 'red')
            unvisited.remove(vertex4)
            # Perform the necessary 2-breaks
            adj = two_break(adj, vertex1, vertex2, vertex4, vertex3)
            adj_p = two_break(adj_p, vertex1, vertex2, vertex4, vertex3)
            # Output adj_p
