# import matplotlib.pyplot as plt
# import networkx.drawing.nx_pylab as nxp
import phylogeny
from phylogeny import add_node
from pathlib import Path
import networkx as nx
from stepik_phylogeny import fetch_small_parsimony_input, \
    fetch_stepik_additive_phylogeny_input, \
    fetch_stepik_input, \
    fetch_stepik_limb_length_input, \
    fetch_stepik_result, \
    fetch_small_parsimony_unrooted, \
    pretty_print_small_parsimony, \
    fetch_nn_input, \
    pretty_print_nn_trees, \
    fetch_small_parsimony_unrooted2


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
    n, j, d = fetch_stepik_limb_length_input(Path('test/testcase34.txt'))
    assert len(d) == n
    length = phylogeny.limb_length(j, d)
    assert length == 6

    n, j, d = fetch_stepik_limb_length_input(Path('test/testcase32.txt'))
    assert len(d) == n
    length = phylogeny.limb_length(j, d)
    assert length == 6

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


def test_neighbor_joining_matrix():
    n, d = fetch_stepik_additive_phylogeny_input(Path('test/testcase33.txt'))
    d_star = phylogeny.neighbor_joining_matrix(d)
    assert d_star == {0: {0: 0, 1: -54, 2: -54, 3: -60}, 1: {0: -54, 1: 0, 2: -60, 3: -54},
                      2: {0: -54, 1: -60, 2: 0, 3: -54}, 3: {0: -60, 1: -54, 2: -54, 3: 0}}

    n, d = fetch_stepik_additive_phylogeny_input(Path('test/testcase13.txt'))
    assert n == len(d)
    d_star = phylogeny.neighbor_joining_matrix(d)
    assert d_star == {0: {0: 0, 1: -68, 2: -60, 3: -60}, 1: {0: -68, 1: 0, 2: -60, 3: -60},
                      2: {0: -60, 1: -60, 2: 0, 3: -68}, 3: {0: -60, 1: -60, 2: -68, 3: 0}}
    tree = phylogeny.neighbor_joining(d)
    assert tree == {0: [(4, 11.0)], 1: [(4, 2.0)], 2: [(5, 6.0)], 3: [(5, 7.0)], 4: [(0, 11.0), (1, 2.0), (5, 4.0)],
                    5: [(2, 6.0), (3, 7.0), (4, 4.0)]}

    n, d = fetch_stepik_additive_phylogeny_input(Path('test/testcase35.txt'))
    tree = phylogeny.neighbor_joining(d)
    assert tree == {0: [(4, 2.25)], 1: [(5, 7.25)], 2: [(4, -0.25)], 3: [(5, 7.75)],
                    4: [(0, 2.25), (2, -0.25), (5, 3.75)], 5: [(1, 7.25), (3, 7.75), (4, 3.75)]}

    n, d = fetch_stepik_additive_phylogeny_input(Path('test/testcase14.txt'))
    assert n == len(d)
    tree = phylogeny.neighbor_joining(d)
    assert tree == {0: [(4, 8.0)], 1: [(5, 13.5)], 2: [(5, 16.5)], 3: [(4, 12.0)], 4: [(0, 8.0), (3, 12.0), (5, 2.0)],
                    5: [(1, 13.5), (2, 16.5), (4, 2.0)]}

    n, d = fetch_stepik_additive_phylogeny_input(Path('test/testcase15.txt'))
    assert n == len(d)
    tree = phylogeny.neighbor_joining(d)
    assert tree == {0: [(37, 535.18)], 1: [(33, 583.051724137931)], 2: [(34, 483.3482142857143)],
                    3: [(39, 559.8097826086956)], 4: [(46, 643.5)], 5: [(35, 565.0648148148148)], 6: [(46, 511.5)],
                    7: [(52, 657.1265625)], 8: [(49, 685.3822115384615)], 9: [(36, 570.7211538461538)],
                    10: [(41, 524.6547619047619)], 11: [(32, 555.9166666666666)], 12: [(32, 601.0833333333334)],
                    13: [(41, 523.3452380952381)], 14: [(44, 456.5208333333333)], 15: [(39, 469.1902173913044)],
                    16: [(42, 500.3875)], 17: [(43, 546.5592105263158)], 18: [(38, 534.9947916666666)],
                    19: [(40, 593.5738636363636)], 20: [(43, 512.4407894736842)], 21: [(42, 565.6125)],
                    22: [(40, 540.4261363636364)], 23: [(35, 553.9351851851852)], 24: [(36, 523.2788461538462)],
                    25: [(53, 590.0833333333334)], 26: [(37, 506.82)], 27: [(33, 520.948275862069)],
                    28: [(53, 563.9166666666666)], 29: [(44, 578.4791666666666)], 30: [(34, 543.6517857142857)],
                    31: [(45, 593.4448529411765)],
                    32: [(11, 555.9166666666666), (12, 601.0833333333334), (38, 66.50520833333334)],
                    33: [(1, 583.051724137931), (27, 520.948275862069), (50, 147.48177083333334)],
                    34: [(2, 483.3482142857143), (30, 543.6517857142857), (51, 165.32102272727272)],
                    35: [(5, 565.0648148148148), (23, 553.9351851851852), (45, 84.55514705882354)],
                    36: [(9, 570.7211538461538), (24, 523.2788461538462), (60, 185.8759765625)],
                    37: [(0, 535.18), (26, 506.82), (48, 90.19642857142857)],
                    38: [(18, 534.9947916666666), (32, 66.50520833333334), (48, 95.05357142857143)],
                    39: [(3, 559.8097826086956), (15, 469.1902173913044), (49, 116.11778846153845)],
                    40: [(19, 593.5738636363636), (22, 540.4261363636364), (47, 99.23333333333333)],
                    41: [(10, 524.6547619047619), (13, 523.3452380952381), (47, 68.26666666666667)],
                    42: [(16, 500.3875), (21, 565.6125), (54, 137.67578125)],
                    43: [(17, 546.5592105263158), (20, 512.4407894736842), (50, 97.01822916666667)],
                    44: [(14, 456.5208333333333), (29, 578.4791666666666), (51, 58.42897727272727)],
                    45: [(31, 593.4448529411765), (35, 84.55514705882354), (59, 156.65494791666666)],
                    46: [(4, 643.5), (6, 511.5), (57, 153.39296875)],
                    47: [(40, 99.23333333333333), (41, 68.26666666666667), (56, 132.25260416666666)],
                    48: [(37, 90.19642857142857), (38, 95.05357142857143), (55, 131.14620535714286)],
                    49: [(8, 685.3822115384615), (39, 116.11778846153845), (52, 49.12343750000002)],
                    50: [(33, 147.48177083333334), (43, 97.01822916666667), (61, 118.2998046875)],
                    51: [(34, 165.32102272727272), (44, 58.42897727272727), (59, 85.72005208333334)],
                    52: [(7, 657.1265625), (49, 49.12343750000002), (55, 40.400669642857146)],
                    53: [(25, 590.0833333333334), (28, 563.9166666666666), (54, 26.32421875)],
                    54: [(42, 137.67578125), (53, 26.32421875), (57, 56.98203125)],
                    55: [(48, 131.14620535714286), (52, 40.400669642857146), (56, 26.677083333333336)],
                    56: [(47, 132.25260416666666), (55, 26.677083333333336), (58, 41.1103515625)],
                    57: [(46, 153.39296875), (54, 56.98203125), (58, 28.9482421875)],
                    58: [(56, 41.1103515625), (57, 28.9482421875), (61, 38.6083984375)],
                    59: [(45, 156.65494791666666), (51, 85.72005208333334), (60, 25.8115234375)],
                    60: [(36, 185.8759765625), (59, 25.8115234375), (61, 14.0595703125)],
                    61: [(50, 118.2998046875), (58, 38.6083984375), (60, 14.0595703125)]}

    n, d = fetch_stepik_additive_phylogeny_input(Path('test/testcase16.txt'))
    assert n == len(d)
    tree = phylogeny.neighbor_joining(d)

    assert tree == {0: [(33, 537.8793103448276)], 1: [(38, 444.2395833333333)], 2: [(39, 546.5434782608696)],
                    3: [(34, 490.35714285714283)], 4: [(42, 527.09375)], 5: [(41, 533.3214285714286)],
                    6: [(35, 531.5092592592592)], 7: [(40, 488.5965909090909)], 8: [(43, 545.6842105263158)],
                    9: [(38, 630.7604166666666)], 10: [(48, 610.6428571428571)], 11: [(35, 506.49074074074076)],
                    12: [(36, 515.3461538461538)], 13: [(44, 525.0763888888889)], 14: [(43, 485.3157894736842)],
                    15: [(33, 523.1206896551724)], 16: [(32, 556.9166666666666)], 17: [(47, 525.225)],
                    18: [(41, 515.6785714285714)], 19: [(32, 482.0833333333333)], 20: [(34, 550.6428571428571)],
                    21: [(37, 600.015)], 22: [(46, 553.75)], 23: [(42, 579.90625)], 24: [(52, 649.103125)],
                    25: [(45, 517.6985294117648)], 26: [(46, 533.25)], 27: [(40, 575.4034090909091)],
                    28: [(45, 539.3014705882352)], 29: [(39, 518.4565217391304)], 30: [(37, 441.985)],
                    31: [(44, 518.9236111111111)],
                    32: [(16, 556.9166666666666), (19, 482.0833333333333), (51, 160.5625)],
                    33: [(0, 537.8793103448276), (15, 523.1206896551724), (36, 44.65384615384616)],
                    34: [(3, 490.35714285714283), (20, 550.6428571428571), (47, 127.775)],
                    35: [(6, 531.5092592592592), (11, 506.49074074074076), (53, 163.63888888888889)],
                    36: [(12, 515.3461538461538), (33, 44.65384615384616), (55, 195.81473214285714)],
                    37: [(21, 600.015), (30, 441.985), (50, 135.40885416666666)],
                    38: [(1, 444.2395833333333), (9, 630.7604166666666), (53, 143.61111111111111)],
                    39: [(2, 546.5434782608696), (29, 518.4565217391304), (48, 95.35714285714283)],
                    40: [(7, 488.5965909090909), (27, 575.4034090909091), (54, 148.1875)],
                    41: [(5, 533.3214285714286), (18, 515.6785714285714), (56, 185.46875)],
                    42: [(4, 527.09375), (23, 579.90625), (59, 196.037109375)],
                    43: [(8, 545.6842105263158), (14, 485.3157894736842), (51, 85.1875)],
                    44: [(13, 525.0763888888889), (31, 518.9236111111111), (49, 90.53846153846153)],
                    45: [(25, 517.6985294117648), (28, 539.3014705882352), (49, 93.46153846153847)],
                    46: [(22, 553.75), (26, 533.25), (54, 110.5625)],
                    47: [(17, 525.225), (34, 127.775), (50, 61.966145833333336)],
                    48: [(10, 610.6428571428571), (39, 95.35714285714283), (52, 51.89687500000002)],
                    49: [(44, 90.53846153846153), (45, 93.46153846153847), (57, 109.5328125)],
                    50: [(37, 135.40885416666666), (47, 61.966145833333336), (58, 142.0791015625)],
                    51: [(32, 160.5625), (43, 85.1875), (57, 116.1546875)],
                    52: [(24, 649.103125), (48, 51.89687500000002), (55, 77.24776785714286)],
                    53: [(35, 163.63888888888889), (38, 143.61111111111111), (56, 56.40625)],
                    54: [(40, 148.1875), (46, 110.5625), (60, 76.83642578125)],
                    55: [(36, 195.81473214285714), (52, 77.24776785714286), (58, 56.4833984375)],
                    56: [(41, 185.46875), (53, 56.40625), (61, 48.02001953125)],
                    57: [(49, 109.5328125), (51, 116.1546875), (59, 21.994140625)],
                    58: [(50, 142.0791015625), (55, 56.4833984375), (61, 25.54248046875)],
                    59: [(42, 196.037109375), (57, 21.994140625), (60, 25.38232421875)],
                    60: [(54, 76.83642578125), (59, 25.38232421875), (61, 9.97607421875)],
                    61: [(56, 48.02001953125), (58, 25.54248046875), (60, 9.97607421875)]}


def test_small_parsimony():
    tree = nx.Graph()
    labels = ['CAAATCCC', 'ATTGCGAC', 'CTGCGCTG', 'ATGGACGA']
    tree.add_nodes_from([(0, {'parent': 4, 'label': labels[0][0]}),
                         (1, {'parent': 4, 'label': labels[1][0]}),
                         (2, {'parent': 5, 'label': labels[2][0]}),
                         (3, {'parent': 5, 'label': labels[3][0]}),
                         (4, {'parent': 6}),
                         (5, {'parent': 6}),
                         (6, {'parent': None})])
    tree.add_edges_from([(4, 0), (4, 1), (5, 2), (5, 3), (6, 4), (6, 5)])

    score = phylogeny.small_parsimony_symbol(tree, 'ACGT')
    # nx.readwrite.write_gpickle(tree, Path('test/testcase15.pickle'))
    expected = nx.readwrite.read_gpickle(Path('test/testcase15.pickle'))
    assert nx.is_isomorphic(tree, expected)
    assert score == 2
    # nxp.draw_planar(result, with_labels=True, font_weight='bold')
    # plt.show()

    tree = nx.Graph()
    # Tree below is the same as in test/testcase17.txt
    labels = ['CAAATCCC', 'ATTGCGAC', 'CTGCGCTG', 'ATGGACGA']
    tree.add_nodes_from([(0, {'parent': 4, 'label': labels[0]}),
                         (1, {'parent': 4, 'label': labels[1]}),
                         (2, {'parent': 5, 'label': labels[2]}),
                         (3, {'parent': 5, 'label': labels[3]}),
                         (4, {'parent': 6}),
                         (5, {'parent': 6}),
                         (6, {'parent': None})])
    tree.add_edges_from([(4, 0), (4, 1), (5, 2), (5, 3), (6, 4), (6, 5)])
    score = phylogeny.small_parsimony(tree, 'ACGT')
    # nx.readwrite.write_gpickle(tree, Path('test/testcase17.pickle'))
    expected = nx.readwrite.read_gpickle(Path('test/testcase17.pickle'))
    assert nx.is_isomorphic(tree, expected)
    assert score == 16

    tree = fetch_small_parsimony_input(Path('test/testcase17.txt'))
    score = phylogeny.small_parsimony(tree, 'ACGT')
    assert score == 16
    assert nx.is_isomorphic(tree, expected)

    labels = ['ACGTAGGCCT', 'ATGTAAGACT', 'TCGAGAGCAC', 'TCGAAAGCAT']
    tree.add_nodes_from([(0, {'parent': 4, 'label': labels[0]}),
                         (1, {'parent': 4, 'label': labels[1]}),
                         (2, {'parent': 5, 'label': labels[2]}),
                         (3, {'parent': 5, 'label': labels[3]}),
                         (4, {'parent': 6}),
                         (5, {'parent': 6}),
                         (6, {'parent': None})])
    tree.add_edges_from([(4, 0), (4, 1), (5, 2), (5, 3), (6, 4), (6, 5)])
    score = phylogeny.small_parsimony(tree, 'ACGT')
    assert score == 8
    # nx.readwrite.write_gpickle(tree, Path('test/testcase20.pickle'))
    expected = nx.readwrite.read_gpickle(Path('test/testcase20.pickle'))
    assert nx.is_isomorphic(tree, expected)

    tree = fetch_small_parsimony_input(Path('test/testcase18.txt'))
    score = phylogeny.small_parsimony(tree, 'ACGT')
    assert score == 11342
    # nx.readwrite.write_gpickle(tree, Path('test/testcase18.pickle'))
    expected = nx.readwrite.read_gpickle(Path('test/testcase18.pickle'))
    assert nx.is_isomorphic(tree, expected)

    tree = fetch_small_parsimony_input(Path('test/testcase19.txt'))
    score = phylogeny.small_parsimony(tree, 'ACGT')
    assert score == 11008
    # nx.readwrite.write_gpickle(tree, Path('test/testcase19.pickle'))
    expected = nx.readwrite.read_gpickle(Path('test/testcase19.pickle'))
    assert nx.is_isomorphic(tree, expected)


def test_root_unrooted():
    tree = nx.Graph()
    # Tree below is the same as in test/testcase17.txt, but unrooted
    labels = ['CAAATCCC', 'ATTGCGAC', 'CTGCGCTG', 'ATGGACGA']
    tree.add_nodes_from([(0, {'label': labels[0]}),
                         (1, {'label': labels[1]}),
                         (2, {'label': labels[2]}),
                         (3, {'label': labels[3]}),
                         (4, {}),
                         (5, {})])
    tree.add_edges_from([(4, 0), (4, 1), (5, 2), (5, 3), (4, 5)])
    phylogeny.root_unrooted(tree)
    # nx.readwrite.write_gpickle(tree, Path('test/testcase21.pickle'))
    expected = nx.readwrite.read_gpickle(Path('test/testcase21.pickle'))
    assert nx.is_isomorphic(tree, expected)


def test_small_parsimony_unrooted():
    tree = nx.Graph()
    # Tree below is the same as in test/testcase17.txt
    labels = ['CAAATCCC', 'ATTGCGAC', 'CTGCGCTG', 'ATGGACGA']
    tree.add_nodes_from([(0, {'label': labels[0]}),
                         (1, {'label': labels[1]}),
                         (2, {'label': labels[2]}),
                         (3, {'label': labels[3]}),
                         (4, {}),
                         (5, {})])
    tree.add_edges_from([(4, 0), (4, 1), (5, 2), (5, 3), (4, 5)])
    score = phylogeny.small_parsimony_unrooted(tree, 'ACGT')
    assert score == 16
    # nx.readwrite.write_gpickle(tree, Path('test/testcase25.pickle'))
    expected = nx.readwrite.read_gpickle(Path('test/testcase25.pickle'))
    assert nx.is_isomorphic(tree, expected)

    tree = fetch_small_parsimony_unrooted(Path('test/testcase22.txt'))
    score = phylogeny.small_parsimony_unrooted(tree, 'ACGT')
    assert score == 17
    # nx.readwrite.write_gpickle(tree, Path('test/testcase22.pickle'))
    expected = nx.readwrite.read_gpickle(Path('test/testcase22.pickle'))
    assert nx.is_isomorphic(tree, expected)

    tree = fetch_small_parsimony_unrooted(Path('test/testcase23.txt'))
    score = phylogeny.small_parsimony_unrooted(tree, 'ACGT')
    assert score == 570
    # nx.readwrite.write_gpickle(tree, Path('test/testcase23.pickle'))
    expected = nx.readwrite.read_gpickle(Path('test/testcase23.pickle'))
    assert nx.is_isomorphic(tree, expected)

    tree = fetch_small_parsimony_unrooted(Path('test/testcase24.txt'))
    score = phylogeny.small_parsimony_unrooted(tree, 'ACGT')
    assert score == 441
    # nx.readwrite.write_gpickle(tree, Path('test/testcase24.pickle'))
    expected = nx.readwrite.read_gpickle(Path('test/testcase24.pickle'))
    assert nx.is_isomorphic(tree, expected)


def test_nearest_neighbor_trees():
    a, b, tree = fetch_nn_input(Path('test/testcase26.txt'))
    tree1, tree2 = phylogeny.nearest_neighbor_trees(tree, a, b)
    # nx.readwrite.write_gpickle(tree1, Path('test/testcase26-1.pickle'))
    # nx.readwrite.write_gpickle(tree1, Path('test/testcase26-2.pickle'))
    expected1 = nx.readwrite.read_gpickle(Path('test/testcase26-1.pickle'))
    expected2 = nx.readwrite.read_gpickle(Path('test/testcase26-2.pickle'))
    assert nx.is_isomorphic(tree1, expected1)
    assert nx.is_isomorphic(tree2, expected2)

    a, b, tree = fetch_nn_input(Path('test/testcase27.txt'))
    tree1, tree2 = phylogeny.nearest_neighbor_trees(tree, a, b)
    # nx.readwrite.write_gpickle(tree1, Path('test/testcase27-1.pickle'))
    expected1 = nx.readwrite.read_gpickle(Path('test/testcase27-1.pickle'))
    assert nx.is_isomorphic(tree1, expected1)
    # nx.readwrite.write_gpickle(tree1, Path('test/testcase27-2.pickle'))
    expected2 = nx.readwrite.read_gpickle(Path('test/testcase27-2.pickle'))
    assert nx.is_isomorphic(tree2, expected2)

    a, b, tree = fetch_nn_input(Path('test/testcase28.txt'))
    tree1, tree2 = phylogeny.nearest_neighbor_trees(tree, a, b)
    # nx.readwrite.write_gpickle(tree1, Path('test/testcase28-1.pickle'))
    expected1 = nx.readwrite.read_gpickle(Path('test/testcase28-1.pickle'))
    assert nx.is_isomorphic(tree1, expected1)
    # nx.readwrite.write_gpickle(tree1, Path('test/testcase28-2.pickle'))
    expected2 = nx.readwrite.read_gpickle(Path('test/testcase28-2.pickle'))
    assert nx.is_isomorphic(tree2, expected2)


def test_nearest_neighbor_interchange():
    tree = fetch_small_parsimony_unrooted2(Path('test/testcase29.txt'))
    lp_tree, score = phylogeny.nearest_neighbor_interchange(tree, 'ACGT')
    assert score == 21
    # nx.readwrite.write_gpickle(lp_tree, Path('test/testcase29.pickle'))
    expected = nx.readwrite.read_gpickle(Path('test/testcase29.pickle'))
    assert nx.is_isomorphic(lp_tree, expected)

    tree = fetch_small_parsimony_unrooted2(Path('test/testcase30.txt'))
    lp_tree, score = phylogeny.nearest_neighbor_interchange(tree, 'ACGT')
    assert score == 156
    # nx.readwrite.write_gpickle(lp_tree, Path('test/testcase30.pickle'))
    expected = nx.readwrite.read_gpickle(Path('test/testcase30.pickle'))
    assert nx.is_isomorphic(lp_tree, expected)

    tree = fetch_small_parsimony_unrooted2(Path('test/testcase31.txt'))
    lp_tree, score = phylogeny.nearest_neighbor_interchange(tree, 'ACGT')
    assert score == 163
    # nx.readwrite.write_gpickle(lp_tree, Path('test/testcase31.pickle'))
    expected = nx.readwrite.read_gpickle(Path('test/testcase31.pickle'))
    assert nx.is_isomorphic(lp_tree, expected)


def test_discrepancy():
    d1 = [[0, 13, 16, 10],
          [13, 0, 21, 15],
          [16, 21, 0, 18],
          [10, 15, 18, 0]]

    d2 = [[0, 14, 17, 11],
          [14, 0, 21, 15],
          [17, 21, 0, 18],
          [11, 15, 18, 0]]

    res = phylogeny.discrepancy(d1, d2)
    assert res == 3


def test_cluster_dist():
    d = [[0, 20, 9, 11],
         [20, 0, 17, 11],
         [9, 17, 0, 8],
         [11, 11, 8, 0]]

    dist = phylogeny.clusters_distance([0, 1], [2, 3], d)
    assert dist == 12
