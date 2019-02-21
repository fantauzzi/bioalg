import numpy as np


def dp_change(money, coins):
    """
    Returns the smallest number of coins that allow to change a given amount of money using available coin denominations.
    :param money: The given amount of money, an integer number.
    :param coins: The available coin denominations, a sequence of integer numbers.
    :return: The smalles number of coins to change the money, an integer.
    """
    min_num_coins = [0]
    for m in range(1, money + 1):
        min_num_coins.append(float('inf'))
        for i in range(0, len(coins)):
            if m >= coins[i]:
                if min_num_coins[m - coins[i]] + 1 < min_num_coins[m]:
                    min_num_coins[m] = min_num_coins[m - coins[i]] + 1
    return min_num_coins[money]


def manhattan_tourist(down, right):
    """
    Returns the longest path for the Manhattan tourist problem.
    :param down: A matrix of the distances to travel downward in the map, a numpy array of integer numbers.
    :param right: A matrix of the distances to travel right in the map, a numpy array of integer numbers.
    :return: The length of the longest path, an integer.
    """
    n, _ = down.shape
    _, m = right.shape
    s = np.zeros((n + 1, m + 1), dtype=np.int)
    for i in range(1, n + 1):
        s[i, 0] = s[i - 1, 0] + down[i - 1, 0]
    for j in range(1, m + 1):
        s[0, j] = s[0, j - 1] + right[0, j - 1]

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s[i, j] = max(s[i - 1, j] + down[i - 1, j], s[i, j - 1] + right[i][j - 1])
    return s[-1, -1]


def topological_ordering(adj):
    topo_order = []

    fan_in = {}
    for vertex, adjs in adj.items():
        if fan_in.get(vertex) is None:
            fan_in[vertex] = 0
        for adj_vertex, _ in adjs:
            current_count = fan_in.get(adj_vertex, 0)
            fan_in[adj_vertex] = current_count + 1

    pending = [vertex for vertex, count in fan_in.items() if count == 0]  # All vertices with no incoming edges
    while pending:
        current = pending.pop()
        topo_order.append(current)
        adjs = adj.get(current)
        if adjs is None:
            continue
        for vertex, _ in adjs:
            vertex_fan_in = fan_in[vertex]
            vertex_fan_in -= 1
            fan_in[vertex] = vertex_fan_in
            if vertex_fan_in == 0:
                pending.append(vertex)

    max_fan_in = max(fan_in.values())
    min_fan_in = min(fan_in.values())
    assert max_fan_in == min_fan_in == 0

    return topo_order


def dag_longest_path(adj, source, sink):
    topo_order = topological_ordering(adj)
    start_idx = topo_order.index(source)
    backtrack = {}
    for current in topo_order[start_idx:]:
        if current == sink:
            break
        adjs = adj.get(current, None)
        if adjs is None:
            continue
        _, current_max_dist = backtrack.get(current, (None, None))
        if current_max_dist is None:
            current_max_dist = 0 if current == source else float('-inf')
        for vertex, distance in adjs:
            _, vertex_max_dist = backtrack.get(vertex, (None, None))
            if vertex_max_dist is None or vertex_max_dist < current_max_dist + distance:
                backtrack[vertex] = (current, current_max_dist + distance)

    longest_path = [sink]
    current = sink
    while current != source:
        current, _ = backtrack.get(current, (None, None))
        if current is None:
            longest_path = None
            break
        longest_path.append(current)
    if longest_path is not None:
        longest_path.reverse()

    _, path_length = backtrack.get(sink, (None, None))

    return path_length, longest_path
