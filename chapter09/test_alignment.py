from pathlib import Path
import pickle
import alignment
from alignment import Edge


# from stepik_alignment import pretty_print_trie


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
    trie, leaf_labels = alignment.suffix_trie_from_text(text)
    assert trie == {0: {'p': Edge(node=1, weight=0), 'a': Edge(node=15, weight=1), 'n': Edge(node=28, weight=2),
                        'm': Edge(node=50, weight=4), 'b': Edge(node=68, weight=6), 's': Edge(node=90, weight=12),
                        '$': Edge(node=92, weight=13)}, 1: {'a': Edge(node=2, weight=1)},
                    2: {'n': Edge(node=3, weight=2)}, 3: {'a': Edge(node=4, weight=3)},
                    4: {'m': Edge(node=5, weight=4)}, 5: {'a': Edge(node=6, weight=5)},
                    6: {'b': Edge(node=7, weight=6)}, 7: {'a': Edge(node=8, weight=7)},
                    8: {'n': Edge(node=9, weight=8)}, 9: {'a': Edge(node=10, weight=9)},
                    10: {'n': Edge(node=11, weight=10)}, 11: {'a': Edge(node=12, weight=11)},
                    12: {'s': Edge(node=13, weight=12)}, 13: {'$': Edge(node=14, weight=13)}, 14: {},
                    15: {'n': Edge(node=16, weight=2), 'm': Edge(node=40, weight=4), 'b': Edge(node=60, weight=6),
                         's': Edge(node=88, weight=12)}, 16: {'a': Edge(node=17, weight=3)},
                    17: {'m': Edge(node=18, weight=4), 'n': Edge(node=76, weight=10), 's': Edge(node=84, weight=12)},
                    18: {'a': Edge(node=19, weight=5)}, 19: {'b': Edge(node=20, weight=6)},
                    20: {'a': Edge(node=21, weight=7)}, 21: {'n': Edge(node=22, weight=8)},
                    22: {'a': Edge(node=23, weight=9)}, 23: {'n': Edge(node=24, weight=10)},
                    24: {'a': Edge(node=25, weight=11)}, 25: {'s': Edge(node=26, weight=12)},
                    26: {'$': Edge(node=27, weight=13)}, 27: {}, 28: {'a': Edge(node=29, weight=3)},
                    29: {'m': Edge(node=30, weight=4), 'n': Edge(node=80, weight=10), 's': Edge(node=86, weight=12)},
                    30: {'a': Edge(node=31, weight=5)}, 31: {'b': Edge(node=32, weight=6)},
                    32: {'a': Edge(node=33, weight=7)}, 33: {'n': Edge(node=34, weight=8)},
                    34: {'a': Edge(node=35, weight=9)}, 35: {'n': Edge(node=36, weight=10)},
                    36: {'a': Edge(node=37, weight=11)}, 37: {'s': Edge(node=38, weight=12)},
                    38: {'$': Edge(node=39, weight=13)}, 39: {}, 40: {'a': Edge(node=41, weight=5)},
                    41: {'b': Edge(node=42, weight=6)}, 42: {'a': Edge(node=43, weight=7)},
                    43: {'n': Edge(node=44, weight=8)}, 44: {'a': Edge(node=45, weight=9)},
                    45: {'n': Edge(node=46, weight=10)}, 46: {'a': Edge(node=47, weight=11)},
                    47: {'s': Edge(node=48, weight=12)}, 48: {'$': Edge(node=49, weight=13)}, 49: {},
                    50: {'a': Edge(node=51, weight=5)}, 51: {'b': Edge(node=52, weight=6)},
                    52: {'a': Edge(node=53, weight=7)}, 53: {'n': Edge(node=54, weight=8)},
                    54: {'a': Edge(node=55, weight=9)}, 55: {'n': Edge(node=56, weight=10)},
                    56: {'a': Edge(node=57, weight=11)}, 57: {'s': Edge(node=58, weight=12)},
                    58: {'$': Edge(node=59, weight=13)}, 59: {}, 60: {'a': Edge(node=61, weight=7)},
                    61: {'n': Edge(node=62, weight=8)}, 62: {'a': Edge(node=63, weight=9)},
                    63: {'n': Edge(node=64, weight=10)}, 64: {'a': Edge(node=65, weight=11)},
                    65: {'s': Edge(node=66, weight=12)}, 66: {'$': Edge(node=67, weight=13)}, 67: {},
                    68: {'a': Edge(node=69, weight=7)}, 69: {'n': Edge(node=70, weight=8)},
                    70: {'a': Edge(node=71, weight=9)}, 71: {'n': Edge(node=72, weight=10)},
                    72: {'a': Edge(node=73, weight=11)}, 73: {'s': Edge(node=74, weight=12)},
                    74: {'$': Edge(node=75, weight=13)}, 75: {}, 76: {'a': Edge(node=77, weight=11)},
                    77: {'s': Edge(node=78, weight=12)}, 78: {'$': Edge(node=79, weight=13)}, 79: {},
                    80: {'a': Edge(node=81, weight=11)}, 81: {'s': Edge(node=82, weight=12)},
                    82: {'$': Edge(node=83, weight=13)}, 83: {}, 84: {'$': Edge(node=85, weight=13)}, 85: {},
                    86: {'$': Edge(node=87, weight=13)}, 87: {}, 88: {'$': Edge(node=89, weight=13)}, 89: {},
                    90: {'$': Edge(node=91, weight=13)}, 91: {}, 92: {}}

    assert leaf_labels == {14: 0, 27: 1, 39: 2, 49: 3, 59: 4, 67: 5, 75: 6, 79: 7, 83: 8, 85: 9, 87: 10, 89: 11, 91: 12,
                           92: 13}
