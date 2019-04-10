from pathlib import Path
import pickle
import alignment
from stepik_alignment import pretty_print_trie_adj_lists, serialise_suffix_tree, NodeInfo, ChildInfo


def test_trie_from_strings():
    with open(Path('test/testcase01.txt')) as input_file:
        strings = input_file.readlines()
        strings = [string.rstrip('\n') for string in strings]

    trie = alignment.trie_from_strings(strings)
    assert trie == {0: {'A': 1, 'G': 7}, 1: {'T': 2}, 2: {'A': 3, 'C': 6}, 3: {'G': 4}, 4: {'A': 5}, 5: {}, 6: {},
                    7: {'A': 8}, 8: {'T': 9}, 9: {}}

    with open(Path('test/testcase02.txt')) as input_file:
        strings = input_file.readlines()
        strings = [string.rstrip('\n') for string in strings]
    trie = alignment.trie_from_strings(strings)
    # with open(Path('test/testcase02.pickle'), 'wb') as f:
    #    pickle.dump(trie, f, pickle.HIGHEST_PROTOCOL)
    with open(Path('test/testcase02.pickle'), 'rb') as pickled:
        expected = pickle.load(pickled)
    assert trie == expected

    with open(Path('test/testcase03.txt')) as input_file:
        strings = input_file.readlines()
        strings = [string.rstrip('\n') for string in strings]
    trie = alignment.trie_from_strings(strings)
    # with open(Path('test/testcase03.pickle'), 'wb') as f:
    #     pickle.dump(trie, f, pickle.HIGHEST_PROTOCOL)
    with open(Path('test/testcase03.pickle'), 'rb') as pickled:
        expected = pickle.load(pickled)
    assert trie == expected


def test_prefix_trie_matching():
    trie = {0: {'A': 1, 'G': 7},
            1: {'T': 2},
            2: {'A': 3, 'C': 6},
            3: {'G': 4},
            4: {'A': 5},
            5: {},
            6: {},
            7: {'A': 8},
            8: {'T': 9},
            9: {}}
    match = alignment.prefix_trie_matching('ATAGA', trie)
    assert match == 'ATAGA'

    match = alignment.prefix_trie_matching('ATC', trie)
    assert match == 'ATC'

    match = alignment.prefix_trie_matching('GAT', trie)
    assert match == 'GAT'

    match = alignment.prefix_trie_matching('GACCCATAGAC', trie)
    assert match == None

    match = alignment.prefix_trie_matching('AT', trie)
    assert match == None


def test_trie_matching():
    trie = {0: {'A': 1, 'G': 7},
            1: {'T': 2},
            2: {'A': 3, 'C': 6},
            3: {'G': 4},
            4: {'A': 5},
            5: {},
            6: {},
            7: {'A': 8},
            8: {'T': 9},
            9: {}}

    matches = alignment.trie_matching('GAT', trie)
    assert matches == [('GAT', 0)]

    matches = alignment.trie_matching('GATTA', trie)
    assert matches == [('GAT', 0)]

    matches = alignment.trie_matching('AATC', trie)
    assert matches == [('ATC', 1)]

    matches = alignment.trie_matching('AATCCGAT', trie)
    assert matches == [('ATC', 1), ('GAT', 5)]

    matches = alignment.trie_matching('ATATGAAT', trie)
    assert matches == []

    trie = alignment.trie_from_strings(['ATCG', 'GGGT'])
    matches = alignment.trie_matching('AATCGGGTTCAATCGGGGT', trie)
    assert matches == [('ATCG', 1), ('GGGT', 4), ('ATCG', 11), ('GGGT', 15)]

    with open(Path('test/testcase04.txt')) as input_file:
        text = input_file.readline().rstrip('\n')
        strings = input_file.readlines()
        strings = [s.rstrip('\n') for s in strings]
    trie = alignment.trie_from_strings(strings)
    matches = alignment.trie_matching(text, trie)
    assert matches == [('TCCTTACTC', 239), ('TCCTTACTC', 246), ('ACTGCTGAC', 383), ('ACTGCTGAC', 390),
                       ('ACTGCTGAC', 669), ('ACTGCTGAC', 676), ('GCTGCGCGC', 963), ('GCGCCACGC', 1757),
                       ('ACTGCTGAC', 1935), ('GCTGCGCGC', 1993), ('GGTCGGGGG', 2065), ('GCTGCGCGC', 2465),
                       ('GTAGTAAGT', 2599), ('GTAGTAAGT', 2606), ('GTAGTAAGT', 2914), ('ACTGCTGAC', 3125),
                       ('GCTGCGCGC', 3626), ('GCTGCGCGC', 3633), ('GCGCCACGC', 3873), ('GGTCGGGGG', 3937),
                       ('GGTCGGGGG', 3944), ('GGTCGGGGG', 4088), ('GGTCGGGGG', 4095), ('TCCTTACTC', 4503),
                       ('TCCTTACTC', 4576), ('GTAGTAAGT', 4978), ('TCCTTACTC', 5104), ('TCCTTACTC', 5111),
                       ('GCGCCACGC', 5419), ('ACTGCTGAC', 5985), ('GTAGTAAGT', 6035), ('GCGCCACGC', 6053),
                       ('GCGCCACGC', 6060), ('GCTGCGCGC', 6079), ('GCGCCACGC', 6659), ('CTCTTCCCT', 7092),
                       ('CTCTTCCCT', 7111), ('CTCTTCCCT', 7118), ('GCTGCGCGC', 7528), ('GCTGCGCGC', 7535),
                       ('CTCTTCCCT', 7613), ('CTCTTCCCT', 7620), ('CTCTTCCCT', 8079)]

    with open(Path('test/testcase05.txt')) as input_file:
        text = input_file.readline().rstrip('\n')
        strings = input_file.readlines()
        strings = [s.rstrip('\n') for s in strings]
    trie = alignment.trie_from_strings(strings)
    matches = alignment.trie_matching(text, trie)
    assert matches == [('TGTCCGGTG', 870), ('TGGAGGATG', 998), ('TGGAGGATG', 1005), ('TGGAGGATG', 1093),
                       ('CAATCCCCA', 1892), ('ATACACGAT', 2218), ('TGTCCGGTG', 2657), ('TGTCCGGTG', 2664),
                       ('ACATCGGAC', 2844), ('ACATCGGAC', 2851), ('TAAAGCATA', 2948), ('TGGAGGATG', 3385),
                       ('AATAATGAA', 3915), ('CAATCCCCA', 4269), ('CAATCCCCA', 4276), ('TGTCCGGTG', 4356),
                       ('TGTCCGGTG', 4363), ('AATAATGAA', 4444), ('AATAATGAA', 4708), ('CAATCCCCA', 5056),
                       ('TAAAGCATA', 5239), ('TAAAGCATA', 5246), ('TGGAGGATG', 5545), ('TGGAGGATG', 5552),
                       ('CAATCCCCA', 5602), ('TGTCCGGTG', 6038), ('ATACACGAT', 6623), ('AATAATGAA', 7504),
                       ('ACATCGGAC', 7677), ('ACATCGGAC', 7684), ('TAAAGCATA', 8019), ('TAAAGCATA', 8026),
                       ('ACATCGGAC', 8274), ('ACATCGGAC', 8281), ('ATACACGAT', 8767), ('AATAATGAA', 8971),
                       ('TGTCCGGTG', 9110), ('TGTCCGGTG', 9117), ('CAATCCCCA', 9325), ('CAATCCCCA', 9332)]


def test_suffix_trie_from_text():
    text = 'panamabananas'
    root = alignment.suffix_trie_from_text(text)
    serialised = serialise_suffix_tree(root)
    # with open(Path('test/testcase06.pickle'), 'wb') as f:
    #    pickle.dump(serialised, f, pickle.HIGHEST_PROTOCOL)
    with open(Path('test/testcase06.pickle'), 'rb') as pickled:
        expected = pickle.load(pickled)
    assert serialised == expected


def test_suffix_tree_from_text():
    text = 'panamabananas'
    root = alignment.suffix_tree_from_text(text)
    serialised = serialise_suffix_tree(root)
    assert serialised == [NodeInfo(data=0, parent_data=None, label=None,
                                   children=[ChildInfo(data=14, symbol='p', weight=None, position=0, length=14),
                                             ChildInfo(data=15, symbol='a', weight=None, position=1, length=1),
                                             ChildInfo(data=29, symbol='n', weight=None, position=2, length=2),
                                             ChildInfo(data=59, symbol='m', weight=None, position=4, length=10),
                                             ChildInfo(data=75, symbol='b', weight=None, position=6, length=8),
                                             ChildInfo(data=91, symbol='s', weight=None, position=12, length=2),
                                             ChildInfo(data=92, symbol='$', weight=None, position=13, length=1)]),
                          NodeInfo(data=14, parent_data=0, label=0, children=[]),
                          NodeInfo(data=15, parent_data=0, label=None,
                                   children=[ChildInfo(data=17, symbol='n', weight=None, position=2, length=2),
                                             ChildInfo(data=49, symbol='m', weight=None, position=4, length=10),
                                             ChildInfo(data=67, symbol='b', weight=None, position=6, length=8),
                                             ChildInfo(data=89, symbol='s', weight=None, position=12, length=2)]),
                          NodeInfo(data=17, parent_data=15, label=None,
                                   children=[ChildInfo(data=27, symbol='m', weight=None, position=4, length=10),
                                             ChildInfo(data=79, symbol='n', weight=None, position=10, length=4),
                                             ChildInfo(data=85, symbol='s', weight=None, position=12, length=2)]),
                          NodeInfo(data=27, parent_data=17, label=1, children=[]),
                          NodeInfo(data=29, parent_data=0, label=None,
                                   children=[ChildInfo(data=39, symbol='m', weight=None, position=4, length=10),
                                             ChildInfo(data=83, symbol='n', weight=None, position=10, length=4),
                                             ChildInfo(data=87, symbol='s', weight=None, position=12, length=2)]),
                          NodeInfo(data=39, parent_data=29, label=2, children=[]),
                          NodeInfo(data=49, parent_data=15, label=3, children=[]),
                          NodeInfo(data=59, parent_data=0, label=4, children=[]),
                          NodeInfo(data=67, parent_data=15, label=5, children=[]),
                          NodeInfo(data=75, parent_data=0, label=6, children=[]),
                          NodeInfo(data=79, parent_data=17, label=7, children=[]),
                          NodeInfo(data=83, parent_data=29, label=8, children=[]),
                          NodeInfo(data=85, parent_data=17, label=9, children=[]),
                          NodeInfo(data=87, parent_data=29, label=10, children=[]),
                          NodeInfo(data=89, parent_data=15, label=11, children=[]),
                          NodeInfo(data=91, parent_data=0, label=12, children=[]),
                          NodeInfo(data=92, parent_data=0, label=13, children=[])]

    text = 'ATAAATG'
    root = alignment.suffix_tree_from_text(text)
    serialised = serialise_suffix_tree(root)
    assert serialised == [NodeInfo(data=0, parent_data=None, label=None,
                                   children=[ChildInfo(data=1, symbol='A', weight=None, position=0, length=1),
                                             ChildInfo(data=9, symbol='T', weight=None, position=1, length=1),
                                             ChildInfo(data=29, symbol='G', weight=None, position=6, length=2),
                                             ChildInfo(data=30, symbol='$', weight=None, position=7, length=1)]),
                          NodeInfo(data=1, parent_data=0, label=None,
                                   children=[ChildInfo(data=2, symbol='T', weight=None, position=1, length=1),
                                             ChildInfo(data=16, symbol='A', weight=None, position=3, length=1)]),
                          NodeInfo(data=2, parent_data=1, label=None,
                                   children=[ChildInfo(data=8, symbol='A', weight=None, position=2, length=6),
                                             ChildInfo(data=25, symbol='G', weight=None, position=6, length=2)]),
                          NodeInfo(data=8, parent_data=2, label=0, children=[]),
                          NodeInfo(data=9, parent_data=0, label=None,
                                   children=[ChildInfo(data=15, symbol='A', weight=None, position=2, length=6),
                                             ChildInfo(data=27, symbol='G', weight=None, position=6, length=2)]),
                          NodeInfo(data=15, parent_data=9, label=1, children=[]),
                          NodeInfo(data=16, parent_data=1, label=None,
                                   children=[ChildInfo(data=20, symbol='A', weight=None, position=4, length=4),
                                             ChildInfo(data=23, symbol='T', weight=None, position=5, length=3)]),
                          NodeInfo(data=20, parent_data=16, label=2, children=[]),
                          NodeInfo(data=23, parent_data=16, label=3, children=[]),
                          NodeInfo(data=25, parent_data=2, label=4, children=[]),
                          NodeInfo(data=27, parent_data=9, label=5, children=[]),
                          NodeInfo(data=29, parent_data=0, label=6, children=[]),
                          NodeInfo(data=30, parent_data=0, label=7, children=[])]


def test_longest_repeat():
    repeat = alignment.longest_repeat('panamabananas')
    assert repeat == 'ana'

    repeat = alignment.longest_repeat('ATATCGTTTTATCGTT')
    assert repeat == 'TATCGTT'

    """
    text = 'AATTTCCGACTTGCATGACGAGTCAGCGTTCCATCTGATCGAGTCTCCGAAGAACAAATACCCCTACTCAGTTGTGAGCCCCTTTACCGTGAGGACAGGGTCCTTGATGTCGTCTCCTAATTTGCGTTGCGGCTCAACATGTTGTACATAGTGGGGCCAGCCCCAGGGATTTTGTAATTTCTACACTCCATATACGGGACAAGGGTGAGCATTTCCGGGCTTGGATAGGGGCTGCAAGAAAATATCTGGACGTAAGAACTTAATGCCATTCCTACATCCTCGATACCTCGTCTGTCAGAGCAATGAGCTGGTTAGAGGACAGTATTGGTCGGTCATCCTCAGATTGGGGACACATCCGTCTCTATGTGCGTTCCGTTGCCTTGTGCTGACCTTGTCGAACGTACCCCATCTTCGAGCCGCACGCTCGACCAGCTAGGTCCCAGCAGTGGCCTGATAGAAAAATTACCTACGGGCCTCCCAATCGTCCTCCCAGGGTGTCGAACTCTCAAAATTCCCGCATGGTCGTGCTTCCGTACGAATTATGCAAACTCCAGAACCCGGATCTATTCCACGCTCAACGAGTCCTTCACGCTTGGTAGAATTTCATGCTCGTCTTTTGTATCCGTGTAAGTAGGAGGCCGCTGTACGGGTATCCCAGCCTTCGCGCTCTGCTGCAGGGACGTTAACACTCCGAACTTTCCATATACGGGACAAGGGTGAGCATTTCCGGGCTTGGATAGGGGCTGCAAGAAAATATCTGGACGTAAGAAGCTCTGAGGGATCCTCACGGAGTTAGATTTATTTTCCATATACGGGACAAGGGTGAGCATTTCCGGGCTTGGATAGGGGCTGCAAGAAAATATCTGGACGTAAGAAGAGTGATGTTTGGAATGCCAACTTCCATGCACGCCAATTGAGCAATCAGGAGAATCGAGTGCTGTTGACCTAGACCTTGTCAGAAGTATGAATTAACCGCGCGTGTAGGTTTGTCGCTCGACCTGCAAGGGTGCACAATCTGGACTGTCGTCGGCGAACGCTTTCATACGCCTACAAACCGCGTTGCTGGTCGAATCGATCTCACCACCGGCCTTGCAGGATTCTAATTATTCTCTCTCGGTGAGACTGCCGGCGGTCCATGGGTCTGTGTTTCGCTTCAAGCAGTGATATACTGGCGTTTTGTGACACATGGCCACGCACGCCTCTCGTTACTCCCAAT'
    repeat = alignment.longest_repeat(text)
    assert repeat == 'TTTCCATATACGGGACAAGGGTGAGCATTTCCGGGCTTGGATAGGGGCTGCAAGAAAATATCTGGACGTAAGAAG'
    """

    text = 'ATGGTAGGTGCCCGGCAGTAATGTGAAGCTCAAGCATCATTTGGACGCCCGGGGTACCTAAGAATTTATAGGCAAGCTCACTAAGAATG'
    repeat = alignment.longest_repeat(text)
    assert repeat == 'CTAAGAAT'


def test_color_suffix_tree():
    root = alignment.suffix_tree_from_text('panama#banana')
    counts = alignment.color_suffix_tree(root, no_of_blue_leaves=7)
    pass