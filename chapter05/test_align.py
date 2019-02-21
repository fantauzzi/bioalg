import align
import numpy as np
from pathlib import Path


def read_weighted_graph(file_name):
    with open(file_name) as input_file:
        lines = input_file.readlines()
    adj = {}
    for line in lines:
        left, right = line.rstrip().split('->')
        vertex = left
        vertex2, weight_s = right.split(':')
        weight = int(weight_s)
        current = adj.get(vertex, [])
        current.append((vertex2, weight))
        adj[vertex] = current
    return adj


def test_dp_change():
    change = align.dp_change(40, [50, 25, 20, 10, 5, 1])
    assert change == 2

    change2 = align.dp_change(19415, [18, 16, 7, 5, 3, 1])
    assert change2 == 1080

    change3 = align.dp_change(16042, [16, 13, 11, 8, 7, 5, 3, 1])
    assert change3 == 1003


def test_manhattan_tourist():
    down = np.array([[1, 0, 2, 4, 3],
                     [4, 6, 5, 2, 1],
                     [4, 4, 5, 2, 1],
                     [5, 6, 8, 5, 3]])

    right = np.array([[3, 2, 4, 0],
                      [3, 2, 4, 2],
                      [0, 7, 3, 3],
                      [3, 3, 0, 2],
                      [1, 3, 2, 2]])

    dist = align.manhattan_tourist(down=down, right=right)
    assert dist == 34

    down2 = np.array([[4, 4, 2, 1, 3, 1, 0, 0, 1],
                      [4, 3, 0, 2, 0, 4, 3, 4, 4],
                      [2, 3, 3, 1, 2, 1, 2, 2, 0],
                      [3, 0, 3, 3, 2, 1, 1, 3, 4],
                      [2, 3, 1, 2, 2, 0, 2, 3, 2],
                      [2, 2, 2, 0, 4, 2, 1, 0, 3],
                      [4, 3, 1, 1, 0, 1, 1, 4, 2],
                      [0, 0, 2, 2, 2, 1, 2, 4, 2],
                      [4, 3, 0, 3, 1, 3, 2, 3, 1],
                      [1, 4, 1, 0, 3, 4, 1, 2, 1],
                      [4, 4, 0, 4, 1, 4, 3, 1, 2],
                      [4, 1, 2, 3, 1, 3, 3, 3, 0],
                      [3, 1, 0, 2, 2, 0, 4, 4, 0],
                      [2, 0, 1, 0, 0, 3, 1, 1, 1],
                      [0, 1, 3, 2, 2, 2, 1, 2, 1],
                      [0, 2, 0, 3, 1, 2, 2, 4, 2],
                      [2, 0, 4, 1, 3, 3, 2, 4, 0],
                      [2, 3, 1, 3, 4, 2, 1, 4, 4]], dtype=np.int)

    right2 = np.array([[3, 3, 1, 1, 3, 4, 4, 4],
                       [4, 0, 3, 1, 0, 3, 4, 4],
                       [2, 2, 2, 3, 3, 1, 1, 4],
                       [1, 3, 1, 4, 4, 2, 0, 1],
                       [0, 2, 0, 3, 3, 3, 1, 0],
                       [3, 2, 0, 4, 1, 4, 4, 3],
                       [3, 0, 1, 1, 0, 3, 3, 0],
                       [3, 1, 1, 0, 2, 3, 4, 0],
                       [2, 4, 2, 1, 1, 3, 1, 2],
                       [1, 0, 4, 3, 0, 3, 3, 0],
                       [2, 3, 2, 4, 4, 3, 3, 0],
                       [3, 1, 2, 0, 3, 4, 3, 2],
                       [0, 0, 4, 1, 4, 0, 3, 4],
                       [3, 2, 3, 2, 0, 1, 2, 1],
                       [4, 3, 3, 2, 0, 1, 1, 2],
                       [0, 0, 4, 1, 2, 4, 0, 3],
                       [3, 4, 0, 1, 2, 3, 0, 1],
                       [4, 0, 2, 4, 2, 2, 4, 0],
                       [4, 3, 4, 2, 2, 3, 2, 3]], dtype=np.int)

    dist2 = align.manhattan_tourist(down=down2, right=right2)
    assert dist2 == 80


def test_dag_longest_path():
    adj = {'0': [('1', 10), ('2', 4)],
           '2': [('3', 2), ('5', 1)],
           '1': [('4', 1)],
           '3': [('4', 3)],
           '5': [('3', 4)]}

    topo_order = align.topological_ordering(adj)
    assert topo_order == ['0', '2', '5', '3', '1', '4']
    distance, path = align.dag_longest_path(adj, '0', '4')
    assert distance == 12
    assert path == ['0', '2', '5', '3', '4']

    adj2 = read_weighted_graph(Path('test/testcase01.txt'))
    topo_order2 = align.topological_ordering(adj2)
    assert topo_order2 == ['0', '2', '3', '1', '4']
    distance2, path2 = align.dag_longest_path(adj2, '0', '4')
    assert distance2 == 9
    assert path2 == ['0', '2', '3', '4']

    adj3 = read_weighted_graph(Path('test/testcase03.txt'))
    distance3, path3 = align.dag_longest_path(adj3, '5', '8')
    topo_order3 = align.topological_ordering(adj3)
    assert topo_order3 == ['4', '3', '2', '1', '5', '7', '6', '8']
    assert distance3 == 5
    assert path3 == ['5', '6', '8']

    adj4 = read_weighted_graph(Path('test/testcase02.txt'))
    topo_order4 = align.topological_ordering(adj4)
    assert topo_order4 == ['9', '10', '3', '19', '2', '6', '4', '0', '7', '13', '8', '1', '11', '5', '14', '15', '23',
                           '26', '21', '17', '27', '12', '25', '28', '20', '24', '16', '18']
    distance4, path4 = align.dag_longest_path(adj4, '5', '20')
    assert path4 == ['5', '14', '15', '20']
    assert distance4 == 36

    adj5 = read_weighted_graph(Path('test/testcase04.txt'))
    topo_order5 = align.topological_ordering(adj5)
    assert topo_order5 == ['2', '9', '5', '8', '3', '6', '19', '7', '0', '4', '11', '1', '12', '18', '13', '21', '22',
                           '26', '10', '14', '15', '16', '17', '20', '23', '24', '27', '25']
    distance5, path5 = align.dag_longest_path(adj5, '12', '27')
    assert distance5 == 94
    assert path5 == ['12', '13', '21', '22', '24', '27']
