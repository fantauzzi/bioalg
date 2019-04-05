from pathlib import Path
import pickle
import alignment
from stepik_alignment import pretty_print_trie


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
