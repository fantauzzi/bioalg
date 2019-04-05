from pathlib import Path
import pickle
import alignment


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
