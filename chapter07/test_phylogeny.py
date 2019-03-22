import phylogeny
from phylogeny import add_node
from pathlib import Path


def pretty_print_matrix(matrix):
    items = sorted(matrix)
    for row_item in items:
        line = [matrix[row_item][col_item] for col_item in items]
        print(*line, sep=' ')


def parse_stepik_input(file_name):
    tree = {}
    with open(file_name) as input_file:
        input_file.readline().rstrip('\n')
        lines = input_file.readlines()
        for line in lines:
            line = line.rstrip('\n')
            left, right = line.split('->')
            node1 = int(left)
            right_node2, right_weight = right.split(':')
            node2 = int(right_node2)
            weight = int(right_weight)
            adjs = tree.get(node1, [])
            adjs.append((node2, weight))
            tree[node1] = adjs

    return tree


def parse_stepik_result(file_name):
    matrix = {}
    with open(file_name) as input_file:
        lines = input_file.readlines()

    for row_i, line in enumerate(lines):
        matrix[row_i] = {}
        line = line.rstrip('\n')
        for col_i, item in enumerate(line.split(' ')):
            matrix[row_i][col_i] = int(item)

    return matrix


def test_add_node():
    tree = {}
    add_node(tree, 0, 4, 11)
    add_node(tree, 1, 4, 2)
    add_node(tree, 4, 5, 4)
    add_node(tree, 5, 2, 6)
    add_node(tree, 5, 3, 7)
    assert tree == {0: [(4, 11)], 4: [(0, 11), (1, 2), (5, 4)], 1: [(4, 2)], 5: [(4, 4), (2, 6), (3, 7)], 2: [(5, 6)],
                    3: [(5, 7)]}


def test_dist_between_leaves():
    tree = {0: [(4, 11)], 4: [(0, 11), (1, 2), (5, 4)], 1: [(4, 2)], 5: [(4, 4), (2, 6), (3, 7)], 2: [(5, 6)],
            3: [(5, 7)]}
    dist = phylogeny.dist_between_leaves(tree)
    assert dist == {0: {0: 0, 1: 13, 2: 21, 3: 22}, 1: {1: 0, 0: 13, 2: 12, 3: 13}, 2: {2: 0, 0: 21, 1: 12, 3: 13},
                    3: {3: 0, 0: 22, 1: 13, 2: 13}}

    tree = parse_stepik_input(Path('test/testcase00.txt'))
    dist = phylogeny.dist_between_leaves(tree)
    assert dist == {0: {0: 0, 1: 13, 3: 22, 2: 21}, 1: {1: 0, 0: 13, 3: 13, 2: 12}, 3: {3: 0, 0: 22, 1: 13, 2: 13},
                    2: {2: 0, 0: 21, 1: 12, 3: 13}}

    tree = parse_stepik_input(Path('test/testcase01.txt'))
    dist = phylogeny.dist_between_leaves(tree)
    expected = parse_stepik_result(Path('test/testresult01.txt'))
    assert dist == expected

    tree = parse_stepik_input(Path('test/testcase02.txt'))
    dist = phylogeny.dist_between_leaves(tree)
    expected = parse_stepik_result(Path('test/testresult02.txt'))
    assert dist == expected
