from pathlib import Path
import pickle
import alignment
from stepik_alignment import serialise_suffix_tree, NodeInfo, ChildInfo, fetch_sequence_of_int, fetch_string, \
    fetch_BW_matching_input, fetch_find_all_input, fetch_approx_match_input


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
    text = 'panamabananas$'
    root = alignment.suffix_trie_from_text(text)
    serialised = serialise_suffix_tree(root)
    # with open(Path('test/testcase06.pickle'), 'wb') as f:
    #    pickle.dump(serialised, f, pickle.HIGHEST_PROTOCOL)
    with open(Path('test/testcase06.pickle'), 'rb') as pickled:
        expected = pickle.load(pickled)
    assert serialised == expected


def test_suffix_tree_from_text():
    text = 'panamabananas$'
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

    text = 'ATAAATG$'
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
    root = alignment.suffix_tree_from_text('panama#bananas$')
    counts = alignment.color_suffix_tree(root, no_of_blue_leaves=7)
    assert counts == {'red': 8, 'blue': 7, 'purple': 4}
    serialised = serialise_suffix_tree(root)
    assert serialised == [NodeInfo(data=0, parent_data=None, label=None,
                                   children=[ChildInfo(data=15, symbol='p', weight=None, position=0, length=15),
                                             ChildInfo(data=16, symbol='a', weight=None, position=1, length=1),
                                             ChildInfo(data=31, symbol='n', weight=None, position=2, length=2),
                                             ChildInfo(data=64, symbol='m', weight=None, position=4, length=11),
                                             ChildInfo(data=82, symbol='#', weight=None, position=6, length=9),
                                             ChildInfo(data=90, symbol='b', weight=None, position=7, length=8),
                                             ChildInfo(data=106, symbol='s', weight=None, position=13, length=2),
                                             ChildInfo(data=107, symbol='$', weight=None, position=14, length=1)]),
                          NodeInfo(data=15, parent_data=0, label=0, children=[]),
                          NodeInfo(data=16, parent_data=0, label=None,
                                   children=[ChildInfo(data=18, symbol='n', weight=None, position=2, length=2),
                                             ChildInfo(data=53, symbol='m', weight=None, position=4, length=11),
                                             ChildInfo(data=73, symbol='#', weight=None, position=6, length=9),
                                             ChildInfo(data=104, symbol='s', weight=None, position=13, length=2)]),
                          NodeInfo(data=18, parent_data=16, label=None,
                                   children=[ChildInfo(data=29, symbol='m', weight=None, position=4, length=11),
                                             ChildInfo(data=94, symbol='n', weight=None, position=11, length=4),
                                             ChildInfo(data=100, symbol='s', weight=None, position=13, length=2)]),
                          NodeInfo(data=29, parent_data=18, label=1, children=[]),
                          NodeInfo(data=31, parent_data=0, label=None,
                                   children=[ChildInfo(data=42, symbol='m', weight=None, position=4, length=11),
                                             ChildInfo(data=98, symbol='n', weight=None, position=11, length=4),
                                             ChildInfo(data=102, symbol='s', weight=None, position=13, length=2)]),
                          NodeInfo(data=42, parent_data=31, label=2, children=[]),
                          NodeInfo(data=53, parent_data=16, label=3, children=[]),
                          NodeInfo(data=64, parent_data=0, label=4, children=[]),
                          NodeInfo(data=73, parent_data=16, label=5, children=[]),
                          NodeInfo(data=82, parent_data=0, label=6, children=[]),
                          NodeInfo(data=90, parent_data=0, label=7, children=[]),
                          NodeInfo(data=94, parent_data=18, label=8, children=[]),
                          NodeInfo(data=98, parent_data=31, label=9, children=[]),
                          NodeInfo(data=100, parent_data=18, label=10, children=[]),
                          NodeInfo(data=102, parent_data=31, label=11, children=[]),
                          NodeInfo(data=104, parent_data=16, label=12, children=[]),
                          NodeInfo(data=106, parent_data=0, label=13, children=[]),
                          NodeInfo(data=107, parent_data=0, label=14, children=[])]


def test_longest_shared_substring():
    s1 = 'TCGGTAGATTGCGCCCACTC'
    s2 = 'AGGGGCTCGCAGTGTAAGAA'
    substring = alignment.longest_shared_substring(s1, s2)
    assert substring == 'TCG'

    with open(Path('test/testcase07.txt')) as input_file:
        s1 = input_file.readline().rstrip('\n')
        s2 = input_file.readline().rstrip('\n')
    substring = alignment.longest_shared_substring(s1, s2)
    assert substring == 'CATCTGGT'

    with open(Path('test/testcase08.txt')) as input_file:
        s1 = input_file.readline().rstrip('\n')
        s2 = input_file.readline().rstrip('\n')
    substring = alignment.longest_shared_substring(s1, s2)
    assert substring == 'GTTCAGACG'


def test_shortes_substring_not_appearing():
    s1 = 'CCAAGCTGCTAGAGG'
    s2 = 'CATGCTGGGCTGGCT'
    res = alignment.shortest_not_shared(s1, s2)
    assert res == 'CC'

    with open(Path('test/testcase09.txt')) as input_file:
        s1 = input_file.readline().rstrip('\n')
        s2 = input_file.readline().rstrip('\n')
    res = alignment.shortest_not_shared(s1, s2)
    assert res == 'GCGAT'

    with open(Path('test/testcase10.txt')) as input_file:
        s1 = input_file.readline().rstrip('\n')
        s2 = input_file.readline().rstrip('\n')
    res = alignment.shortest_not_shared(s1, s2)
    assert res == 'GAAAT'


def test_suffix_array_for_text():
    suffixes = alignment.suffix_array_for_text('papaya$')
    print()
    print(*suffixes, sep=' ')

    suffixes = alignment.suffix_array_for_text('AACGATAGCGGTAGA$')
    assert suffixes == [15, 14, 0, 1, 12, 6, 4, 2, 8, 13, 3, 7, 9, 10, 11, 5]

    text = fetch_string(Path('test/testcase11.txt'))
    suffixes = alignment.suffix_array_for_text(text)
    expecetd = fetch_sequence_of_int(Path('test/expected11.txt'))
    assert suffixes == expecetd

    text = fetch_string(Path('test/testcase12.txt'))
    suffixes = alignment.suffix_array_for_text(text)
    expecetd = fetch_sequence_of_int(Path('test/expected12.txt'))
    assert suffixes == expecetd


def test_burrows_wheeler_transform():
    trans = alignment.burrows_wheeler_transform('TCAGGGCTTG$')
    assert trans == 'GCTGTGGA$TC'

    trans = alignment.burrows_wheeler_transform('panamabananas$')
    assert trans == 'smnpbnnaaaaa$a'

    trans = alignment.burrows_wheeler_transform('GCGTGCCTGGTCA$')
    assert trans == 'ACTGGCT$TGCGGC'

    trans = alignment.burrows_wheeler_transform('abracadabra$')
    assert trans == 'ard$rcaaaabb'

    text = fetch_string(Path('test/testcase13.txt'))
    expected = fetch_string(Path('test/expected13.txt'))
    trans = alignment.burrows_wheeler_transform(text)
    assert trans == expected

    text = fetch_string(Path('test/testcase14.txt'))
    expected = fetch_string(Path('test/expected14.txt'))
    trans = alignment.burrows_wheeler_transform(text)
    assert trans == expected


def test_inverted_burrow_wheeler():
    trans = alignment.inverted_burrow_wheeler('TTCCATTGGA$')
    assert trans == 'TGTACCATGT$'

    trans = alignment.inverted_burrow_wheeler('ard$rcaaaabb')
    assert trans == 'abracadabra$'

    trans = alignment.inverted_burrow_wheeler('smnpbnnaaaaa$a')
    assert trans == 'panamabananas$'

    text = fetch_string(Path('test/expected13.txt'))
    expected = fetch_string(Path('test/testcase13.txt'))
    trans = alignment.inverted_burrow_wheeler(text)
    assert trans == expected

    text = fetch_string(Path('test/expected14.txt'))
    expected = fetch_string(Path('test/testcase14.txt'))
    trans = alignment.inverted_burrow_wheeler(text)
    assert trans == expected

    text = fetch_string(Path('test/testcase15.txt'))
    expected = fetch_string(Path('test/expected15.txt'))
    trans = alignment.inverted_burrow_wheeler(text)
    assert trans == expected


def test_match_counts():
    text = 'TCCTCTATGAGATCCTATTCTATGAAACCTTCA$GACCAAAATTCTCCGGC'
    pattern = ['CCT', 'CAC', 'GAG', 'CAG', 'ATC']
    counts = alignment.count_matches(text, pattern)
    assert counts == [2, 1, 1, 0, 1]

    text = 'TCCTCTATGAGATCCTATTCTATGAAACCTTCA$GACCAAAATTCTCCGGC'
    pattern = 'CCT'
    counts = alignment.count_matches(text, pattern)
    assert counts == 2

    text, pattern = fetch_BW_matching_input(Path('test/testcase16.txt'))
    counts = alignment.count_matches(text, pattern)
    expected = fetch_sequence_of_int(Path('test/expected16.txt'))
    assert counts == expected

    text, pattern = fetch_BW_matching_input(Path('test/testcase17.txt'))
    counts = alignment.count_matches(text, pattern)
    expected = fetch_sequence_of_int(Path('test/expected17.txt'))
    assert counts == expected


def test_better_count_matches():
    text = 'GGCGCCGC$TAGTCACACACGCCGTA'
    patterns = ['ACC', 'CCG', 'CAG']
    res = alignment.better_matches_count(text, patterns)
    assert res == [1, 2, 1]

    text, pattern = fetch_BW_matching_input(Path('test/testcase18.txt'))
    counts = alignment.better_matches_count(text, pattern)
    expected = fetch_sequence_of_int(Path('test/expected18.txt'))
    assert counts == expected

    text, pattern = fetch_BW_matching_input(Path('test/testcase19.txt'))
    counts = alignment.better_matches_count(text, pattern)
    expected = fetch_sequence_of_int(Path('test/expected19.txt'))
    assert counts == expected


def test_find_all():
    text = 'AATCGGGTTCAATCGGGGT$'
    patterns = ['ATCG', 'GGGT']
    positions = alignment.find_all(text, patterns)
    assert positions == [1, 4, 11, 15]

    text = 'AAAAAA$'
    patterns = ['AA']
    positions = alignment.find_all(text, patterns)
    assert positions == [0, 1, 2, 3, 4]

    text, patterns = fetch_find_all_input(Path('test/testcase20.txt'))
    text = text + '$'
    positions = alignment.find_all(text, patterns)
    expected = fetch_sequence_of_int(Path('test/expected20.txt'))
    assert positions == sorted(expected)

    text, patterns = fetch_find_all_input(Path('test/testcase21.txt'))
    text = text + '$'
    positions = alignment.find_all(text, patterns)
    expected = fetch_sequence_of_int(Path('test/expected21.txt'))
    assert positions == sorted(expected)


def test_find_all_approx():
    text = 'GTATTCTATA$'
    patterns = ['TATT']
    pos = alignment.find_all(text, patterns, d=1)
    assert pos == [1, 6]

    text = 'GTATTCTATA$'
    patterns = ['TATG']
    pos = alignment.find_all(text, patterns, d=1)
    assert pos == [1, 6]

    text = 'ACATGCTACTTT$'
    patterns = ['ATT', 'GCC', 'GCTA', 'TATT']
    pos = alignment.find_all(text, patterns, d=1)
    assert pos == [2, 4, 4, 6, 7, 8, 9]

    text, patterns, d = fetch_approx_match_input(Path('test/testcase22.txt'))
    text = text + '$'
    pos = alignment.find_all(text, patterns, d)
    expected = fetch_sequence_of_int(Path('test/expected22.txt'))
    assert pos == sorted(expected)

    text, patterns, d = fetch_approx_match_input(Path('test/testcase23.txt'))
    text = text + '$'
    pos = alignment.find_all(text, patterns, d)
    expected = fetch_sequence_of_int(Path('test/expected23.txt'))
    assert pos == sorted(expected)
