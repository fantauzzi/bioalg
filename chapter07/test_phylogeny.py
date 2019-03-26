import phylogeny
from phylogeny import add_node
from pathlib import Path


def pretty_print_matrix(matrix):
    """
    Prints the matrix with distances between nodes in a graph, as returned by dist_between_leaves().
    :param matrix: The matrix.
    """
    items = sorted(matrix)
    for row_item in items:
        line = [matrix[row_item][col_item] for col_item in items]
        print(*line, sep=' ')


def pretty_print_adjacencies(tree):
    for node in sorted(tree):
        for adj_node, adj_dist in tree[node]:
            print(node, '->', adj_node, ':', adj_dist, sep='')


def fetch_stepik_input(file_name):
    """
    Fetches the adjacency lists of a tree from a file, with input from the Stepik challenge "Distances Between Leaves
    Problem"
    :param file_name: The file name, with its relative path.
    :return: The adjacency lists, as expected by dist_between_leaves().
    """
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


def parse_stepik_matrix(input_file):
    """
    Fetches and returns a matrix of distances between nodes in a graph from a given file.
    :param file_name: The file name, with its relative path. The content of the file must be formatted like the
    output of the stepik challenge "Distances Between Leaves Problem".
    :return: The distances matrix, a dictionary of dictionaries, same as returned by dist_between_leaves()
    """
    matrix = {}
    lines = input_file.readlines()

    for row_i, line in enumerate(lines):
        matrix[row_i] = {}
        line = line.rstrip('\n')
        for col_i, item in enumerate(line.split(' ')):
            matrix[row_i][col_i] = int(item)

    return matrix


def fetch_stepik_result(file_name):
    with open(file_name) as input_file:
        matrix = parse_stepik_matrix(input_file)
    return matrix


def fetch_stepik_limb_length_input(file_name):
    with open(file_name) as input_file:
        n = input_file.readline().rstrip('\n')
        n = int(n)
        j = input_file.readline().rstrip('\n')
        j = int(j)
        matrix = parse_stepik_matrix(input_file)
    return n, j, matrix


def fetch_stepik_additive_phylogeny_input(file_name):
    with open(file_name) as input_file:
        n = input_file.readline().rstrip('\n')
        n = int(n)
        matrix = parse_stepik_matrix(input_file)
    return n, matrix


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

    tree = fetch_stepik_input(Path('test/testcase00.txt'))
    dist = phylogeny.dist_between_leaves(tree)
    assert dist == {0: {0: 0, 1: 13, 3: 22, 2: 21}, 1: {1: 0, 0: 13, 3: 13, 2: 12}, 3: {3: 0, 0: 22, 1: 13, 2: 13},
                    2: {2: 0, 0: 21, 1: 12, 3: 13}}

    tree = fetch_stepik_input(Path('test/testcase01.txt'))
    dist = phylogeny.dist_between_leaves(tree)
    expected = fetch_stepik_result(Path('test/testresult01.txt'))
    assert dist == expected

    tree = fetch_stepik_input(Path('test/testcase02.txt'))
    dist = phylogeny.dist_between_leaves(tree)
    expected = fetch_stepik_result(Path('test/testresult02.txt'))
    assert dist == expected


def test_limb_length():
    n, j, d = fetch_stepik_limb_length_input(Path('test/testcase03.txt'))
    assert len(d) == n
    length = phylogeny.limb_length(j, d)
    assert length == 2

    n, j, d = fetch_stepik_limb_length_input(Path('test/testcase04.txt'))
    assert len(d) == n
    length = phylogeny.limb_length(j, d)
    assert length == 534

    n, j, d = fetch_stepik_limb_length_input(Path('test/testcase05.txt'))
    assert len(d) == n
    length = phylogeny.limb_length(j, d)
    assert length == 537


def test_path_in_a_tree():
    tree = {0: [(4, 11)], 4: [(0, 11), (1, 2), (5, 4)], 1: [(4, 2)], 5: [(4, 4), (2, 6), (3, 7)], 2: [(5, 6)],
            3: [(5, 7)]}
    path = phylogeny.path_in_a_tree(tree, 0, 3)
    assert path == [(0, 0), (4, 11), (5, 15), (3, 22)]

    path = phylogeny.path_in_a_tree(tree, 0, 1)
    assert path == [(0, 0), (4, 11), (1, 13)]

    path = phylogeny.path_in_a_tree(tree, 1, 3)
    assert path == [(1, 0), (4, 2), (5, 6), (3, 13)]

    tree = fetch_stepik_input(Path('test/testcase01.txt'))
    path = phylogeny.path_in_a_tree(tree, 0, 61)
    assert path == [(0, 0), (54, 15), (55, 21), (56, 29), (58, 39), (59, 44), (60, 52), (61, 57)]
    path2 = phylogeny.path_in_a_tree(tree, 61, 0)
    assert path2 == [(61, 0), (60, 5), (59, 13), (58, 18), (56, 28), (55, 36), (54, 42), (0, 57)]


def fix_nodes_numbering(tree):
    def corrected_neg_node(node, max_pos):
        assert node < 0
        return abs(node) + max_pos

    max_pos = max(tree)
    nodes = set(tree)
    for node in nodes:
        if node < 0:
            tree[corrected_neg_node(node, max_pos)] = tree[node]
            del tree[node]

    for node in tree:
        new_adj = []
        for adj_node, adj_dist in tree[node]:
            new_adj.append((adj_node, adj_dist) if adj_node >= 0 else (corrected_neg_node(adj_node, max_pos), adj_dist))
        tree[node] = new_adj


def test_additive_phylogeny():
    n, j, d = fetch_stepik_limb_length_input(Path('test/testcase03.txt'))
    assert len(d) == n
    tree = phylogeny.additive_phylogeny(d)
    assert tree == {0: [(-1, 11)], 1: [(-1, 2)], -1: [(0, 11), (1, 2), (-2, 4)], 2: [(-2, 6)],
                    -2: [(-1, 4), (2, 6), (3, 7)], 3: [(-2, 7)]}

    fix_nodes_numbering(tree)
    assert tree == {0: [(4, 11)], 1: [(4, 2)], 2: [(5, 6)], 3: [(5, 7)], 5: [(4, 4), (2, 6), (3, 7)],
                    4: [(0, 11), (1, 2), (5, 4)]}

    n, d = fetch_stepik_additive_phylogeny_input(Path('test/testcase06.txt'))
    assert n == len(d)
    tree = phylogeny.additive_phylogeny(d)
    fix_nodes_numbering(tree)
    assert tree == {0: [(29, 745)], 1: [(36, 156)], 2: [(32, 788)], 3: [(30, 409)], 4: [(31, 280)], 5: [(35, 125)],
                    6: [(33, 492)], 7: [(34, 657)], 8: [(53, 311)], 9: [(43, 820)], 10: [(37, 280)], 11: [(38, 723)],
                    12: [(39, 417)], 13: [(40, 864)], 14: [(41, 236)], 15: [(42, 89)], 16: [(51, 713)], 17: [(44, 445)],
                    18: [(45, 87)], 19: [(46, 441)], 20: [(47, 783)], 21: [(48, 348)], 22: [(49, 922)], 23: [(50, 662)],
                    24: [(51, 375)], 25: [(52, 718)], 26: [(55, 868)], 27: [(54, 841)], 28: [(55, 890)],
                    55: [(53, 965), (26, 868), (28, 890)], 54: [(50, 687), (37, 112), (27, 841)],
                    53: [(44, 416), (8, 311), (55, 965)], 52: [(43, 170), (49, 952), (25, 718)],
                    51: [(48, 381), (16, 713), (24, 375)], 50: [(32, 276), (23, 662), (54, 687)],
                    49: [(48, 683), (22, 922), (52, 952)], 48: [(21, 348), (49, 683), (51, 381)],
                    47: [(31, 110), (34, 866), (20, 783)], 46: [(40, 522), (35, 481), (19, 441)],
                    45: [(29, 355), (33, 809), (18, 87)], 44: [(35, 230), (17, 445), (53, 416)],
                    43: [(42, 982), (9, 820), (52, 170)], 42: [(36, 247), (15, 89), (43, 982)],
                    41: [(30, 698), (31, 656), (14, 236)], 40: [(37, 464), (13, 864), (46, 522)],
                    39: [(29, 323), (30, 64), (12, 417)], 30: [(3, 409), (39, 64), (41, 698)],
                    37: [(10, 280), (40, 464), (54, 112)], 36: [(1, 156), (38, 87), (42, 247)],
                    35: [(5, 125), (44, 230), (46, 481)], 34: [(32, 527), (7, 657), (47, 866)],
                    33: [(6, 492), (38, 884), (45, 809)], 32: [(2, 788), (34, 527), (50, 276)],
                    31: [(4, 280), (41, 656), (47, 110)], 29: [(0, 745), (39, 323), (45, 355)],
                    38: [(33, 884), (36, 87), (11, 723)]}

    n, d = fetch_stepik_additive_phylogeny_input(Path('test/testcase07.txt'))
    assert n == len(d)
    tree = phylogeny.additive_phylogeny(d)
    fix_nodes_numbering(tree)
    assert tree == {0: [(23, 613)], 1: [(21, 174)], 2: [(22, 810)], 3: [(25, 977)], 4: [(39, 683)], 5: [(24, 619)],
                    6: [(34, 801)], 7: [(26, 841)], 8: [(27, 816)], 9: [(28, 221)], 10: [(29, 825)], 11: [(30, 140)],
                    12: [(31, 326)], 13: [(32, 652)], 14: [(33, 639)], 15: [(35, 239)], 16: [(35, 528)], 17: [(36, 82)],
                    18: [(37, 194)], 19: [(38, 421)], 20: [(39, 142)], 39: [(28, 70), (4, 683), (20, 142)],
                    38: [(33, 943), (34, 673), (19, 421)], 37: [(23, 856), (24, 612), (18, 194)],
                    36: [(22, 754), (25, 731), (17, 82)], 35: [(34, 538), (15, 239), (16, 528)],
                    34: [(6, 801), (35, 538), (38, 673)], 33: [(32, 307), (14, 639), (38, 943)],
                    32: [(29, 978), (13, 652), (33, 307)], 31: [(26, 413), (21, 901), (12, 326)],
                    22: [(2, 810), (30, 113), (36, 754)], 29: [(27, 679), (10, 825), (32, 978)],
                    28: [(23, 723), (9, 221), (39, 70)], 27: [(25, 111), (8, 816), (29, 679)],
                    26: [(24, 158), (7, 841), (31, 413)], 25: [(3, 977), (27, 111), (36, 731)],
                    24: [(5, 619), (26, 158), (37, 612)], 23: [(0, 613), (28, 723), (37, 856)],
                    21: [(1, 174), (30, 355), (31, 901)], 30: [(21, 355), (22, 113), (11, 140)]}

    n, d = fetch_stepik_additive_phylogeny_input(Path('test/testcase08.txt'))
    assert n == len(d)
    tree = phylogeny.additive_phylogeny(d)
    fix_nodes_numbering(tree)
    assert tree == {0: [(9, 146)], 1: [(9, 149)], 2: [(10, 149)], 3: [(11, 254)], 4: [(12, 404)], 5: [(12, 414)],
                    6: [(13, 489)], 7: [(15, 9)], 8: [(15, 7)], 10: [(9, 11), (2, 149), (11, 86)],
                    15: [(14, 455), (7, 9), (8, 7)], 14: [(11, 249), (13, 22), (15, 455)],
                    13: [(12, 163), (6, 489), (14, 22)], 12: [(4, 404), (5, 414), (13, 163)],
                    11: [(10, 86), (3, 254), (14, 249)], 9: [(0, 146), (1, 149), (10, 11)]}


def test_upgma():
    n, d = fetch_stepik_additive_phylogeny_input(Path('test/testcase10.txt'))
    assert n == len(d)
    tree = phylogeny.upgma(d)
    assert tree == {0: [(5, 1.5)], 1: [(5, 1.5)], 2: [(4, 1.0)], 3: [(4, 1.0)], 4: [(2, 1.0), (3, 1.0), (6, 1.0)],
                    5: [(0, 1.5), (1, 1.5), (6, 0.5)], 6: [(4, 1.0), (5, 0.5)]}

    n, d = fetch_stepik_additive_phylogeny_input(Path('test/testcase09.txt'))
    assert n == len(d)
    tree = phylogeny.upgma(d)
    assert tree == {0: [(5, 7.0)], 1: [(6, 8.833333333333334)], 2: [(4, 5.0)], 3: [(4, 5.0)],
                    4: [(2, 5.0), (3, 5.0), (5, 2.0)], 5: [(0, 7.0), (4, 2.0), (6, 1.833333333333334)],
                    6: [(1, 8.833333333333334), (5, 1.833333333333334)]}

    n, d = fetch_stepik_additive_phylogeny_input(Path('test/testcase11.txt'))
    assert n == len(d)
    tree = phylogeny.upgma(d)
    assert tree == {0: [(28, 25.5)], 1: [(36, 59.5)], 2: [(36, 59.5)], 3: [(37, 67.75)], 4: [(27, 25.0)],
                    5: [(29, 27.0)], 6: [(33, 45.0)], 7: [(27, 25.0)], 8: [(30, 28.0)], 9: [(39, 117.25)],
                    10: [(44, 168.0)], 11: [(35, 45.75)], 12: [(32, 30.0)], 13: [(31, 28.5)],
                    14: [(41, 131.66666666666666)], 15: [(33, 45.0)], 16: [(44, 168.0)], 17: [(28, 25.5)],
                    18: [(38, 68.5)], 19: [(34, 45.5)], 20: [(32, 30.0)], 21: [(34, 45.5)], 22: [(38, 68.5)],
                    23: [(30, 28.0)], 24: [(29, 27.0)], 25: [(31, 28.5)], 26: [(48, 235.375)],
                    27: [(4, 25.0), (7, 25.0), (35, 20.75)], 28: [(0, 25.5), (17, 25.5), (40, 101.125)],
                    29: [(5, 27.0), (24, 27.0), (37, 40.75)], 30: [(8, 28.0), (23, 28.0), (42, 128.125)],
                    31: [(13, 28.5), (25, 28.5), (45, 201.125)], 32: [(12, 30.0), (20, 30.0), (43, 134.375)],
                    33: [(6, 45.0), (15, 45.0), (39, 72.25)], 34: [(19, 45.5), (21, 45.5), (42, 110.625)],
                    35: [(11, 45.75), (27, 20.75), (41, 85.91666666666666)], 36: [(1, 59.5), (2, 59.5), (43, 104.875)],
                    37: [(3, 67.75), (29, 40.75), (49, 189.15476190476187)], 38: [(18, 68.5), (22, 68.5), (40, 58.125)],
                    39: [(9, 117.25), (33, 72.25), (46, 113.41666666666669)],
                    40: [(28, 101.125), (38, 58.125), (45, 103.0)],
                    41: [(14, 131.66666666666666), (35, 85.91666666666666), (50, 127.06770833333334)],
                    42: [(30, 128.125), (34, 110.625), (48, 79.25)],
                    43: [(32, 134.375), (36, 104.875), (46, 66.29166666666669)],
                    44: [(10, 168.0), (16, 168.0), (47, 67.29166666666666)],
                    45: [(31, 201.125), (40, 103.0), (47, 5.666666666666657)],
                    46: [(39, 113.41666666666669), (43, 66.29166666666669), (49, 26.238095238095184)],
                    47: [(44, 67.29166666666666), (45, 5.666666666666657), (50, 23.442708333333343)],
                    48: [(26, 235.375), (42, 79.25), (52, 59.288636363636385)],
                    49: [(37, 189.15476190476187), (46, 26.238095238095184), (51, 29.87440476190477)],
                    50: [(41, 127.06770833333334), (47, 23.442708333333343), (51, 28.04479166666664)],
                    51: [(49, 29.87440476190477), (50, 28.04479166666664), (52, 7.884469696969745)],
                    52: [(48, 59.288636363636385), (51, 7.884469696969745)]}

    n, d = fetch_stepik_additive_phylogeny_input(Path('test/testcase12.txt'))
    assert n == len(d)
    tree = phylogeny.upgma(d)
    assert tree == {0: [(34, 402.5)], 1: [(29, 393.5)], 2: [(42, 493.5)], 3: [(37, 423.0)], 4: [(33, 401.5)],
                    5: [(41, 481.0)], 6: [(32, 395.0)], 7: [(41, 481.0)], 8: [(54, 640.2592592592592)],
                    9: [(38, 425.5)], 10: [(34, 402.5)], 11: [(32, 395.0)], 12: [(30, 394.0)], 13: [(31, 394.5)],
                    14: [(37, 423.0)], 15: [(42, 493.5)], 16: [(28, 392.0)], 17: [(35, 409.0)], 18: [(29, 393.5)],
                    19: [(30, 394.0)], 20: [(36, 415.0)], 21: [(36, 415.0)], 22: [(35, 409.0)], 23: [(31, 394.5)],
                    24: [(38, 425.5)], 25: [(28, 392.0)], 26: [(33, 401.5)], 27: [(46, 536.0)],
                    28: [(16, 392.0), (25, 392.0), (40, 85.375)], 29: [(1, 393.5), (18, 393.5), (43, 102.75)],
                    30: [(12, 394.0), (19, 394.0), (47, 160.75)], 31: [(13, 394.5), (23, 394.5), (39, 60.375)],
                    32: [(6, 395.0), (11, 395.0), (48, 163.91666666666663)], 33: [(4, 401.5), (26, 401.5), (44, 114.5)],
                    34: [(0, 402.5), (10, 402.5), (39, 52.375)], 35: [(17, 409.0), (22, 409.0), (44, 107.0)],
                    36: [(20, 415.0), (21, 415.0), (43, 81.25)], 37: [(3, 423.0), (14, 423.0), (46, 113.0)],
                    38: [(9, 425.5), (24, 425.5), (40, 51.875)],
                    39: [(31, 60.375), (34, 52.375), (50, 126.10416666666663)],
                    40: [(28, 85.375), (38, 51.875), (45, 48.3125)], 41: [(5, 481.0), (7, 481.0), (45, 44.6875)],
                    42: [(2, 493.5), (15, 493.5), (49, 80.6875)], 43: [(29, 102.75), (36, 81.25), (49, 77.9375)],
                    44: [(33, 114.5), (35, 107.0), (52, 78.72500000000002)],
                    45: [(40, 48.3125), (41, 44.6875), (48, 33.22916666666663)],
                    46: [(27, 536.0), (37, 113.0), (47, 18.75)],
                    47: [(30, 160.75), (46, 18.75), (52, 39.97500000000002)],
                    48: [(32, 163.91666666666663), (45, 33.22916666666663), (51, 30.908333333333417)],
                    49: [(42, 80.6875), (43, 77.9375), (50, 6.791666666666629)],
                    50: [(39, 126.10416666666663), (49, 6.791666666666629), (51, 8.845833333333417)],
                    51: [(48, 30.908333333333417), (50, 8.845833333333417), (53, 11.730555555555611)],
                    52: [(44, 78.72500000000002), (47, 39.97500000000002), (53, 6.830555555555634)],
                    53: [(51, 11.730555555555611), (52, 6.830555555555634), (54, 38.70370370370358)],
                    54: [(8, 640.2592592592592), (53, 38.70370370370358)]}
