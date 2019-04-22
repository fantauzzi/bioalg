# import matplotlib.pyplot as plt
from pathlib import Path
import pickle
import networkx as nx
import networkx.algorithms.isomorphism as iso
# import networkx.drawing.nx_pylab as nxp
from stepik_hmm import fetch_hmm, fetch_alignment, ugly_print_matrices, fetch_profile_alignment, pretty_alignment_print
import hmm


def test_make_graph():
    emissions, model = fetch_hmm(Path('test/testcase01.txt'))
    graph = hmm.make_graph(emissions, model)
    # nx.readwrite.write_gpickle(graph, Path('test/testcase01.pickle'))
    expected = nx.readwrite.read_gpickle(Path('test/testcase01.pickle'))
    em = iso.numerical_edge_match('weight', 1)
    # nxp.draw_planar(graph, with_labels=True, font_weight='bold')
    # plt.show()
    assert nx.is_isomorphic(graph, expected, edge_match=em)

    emissions, model = fetch_hmm(Path('test/testcase03.txt'))
    graph = hmm.make_graph(emissions, model)
    # nx.readwrite.write_gpickle(graph, Path('test/testcase03.pickle'))
    expected = nx.readwrite.read_gpickle(Path('test/testcase03.pickle'))
    em = iso.numerical_edge_match('weight', 1)
    assert nx.is_isomorphic(graph, expected, edge_match=em)


def test_viterbi():
    emissions, model = fetch_hmm(Path('test/testcase01.txt'))
    path = hmm.viterbi(emissions=emissions, model=model)
    assert path == 'AAABBAAAAA'

    emissions, model = fetch_hmm(Path('test/testcase02.txt'))
    path = hmm.viterbi(emissions=emissions, model=model)
    assert path == 'AAAAAAAAAAAAAABBBBBBBBBBBAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABBBBBBBBBBBAAAAAAAAAAAAAAAAAAAAABBBBBBBBBBAAA'

    emissions, model = fetch_hmm(Path('test/testcase03.txt'))
    path = hmm.viterbi(emissions=emissions, model=model)
    assert path == 'ABACCBABBABABACCABCCBABAABBBAABABCCBABBABABACCCCCCCCCCBBBBBABACCBABBACCCCCCCCCCCCCCCCBABABACBABAACCC'

    emissions, model = fetch_hmm(Path('test/testcase04.txt'))
    path = hmm.viterbi(emissions=emissions, model=model)
    assert path == 'CCCCCAAAAAAAAABABCAAAAAAABCCCAABAAAAAAAAAAABABAAABAAAAAAAAAAAAABABAAABAAAABAAABCABAAAABCAAABAAABCCCC'


def test_make_profile_HMM():
    theta, sigma, alphabet, alignment = fetch_alignment(Path('test/testcase07.txt'))
    the_HMM = hmm.make_profile_HMM(theta=theta, sigma=sigma, alphabet=alphabet, alignment=alignment)
    # with open(Path('test/testcase07.pickle'), 'wb') as f:
    #    pickle.dump(the_HMM, f, pickle.HIGHEST_PROTOCOL)
    with open(Path('test/testcase07.pickle'), 'rb') as f:
        expected = pickle.load(f)
    assert the_HMM == expected

    theta, sigma, alphabet, alignment = fetch_alignment(Path('test/testcase08.txt'))
    the_HMM = hmm.make_profile_HMM(theta=theta, sigma=sigma, alphabet=alphabet, alignment=alignment)
    # with open(Path('test/testcase08.pickle'), 'wb') as f:
    #    pickle.dump(the_HMM, f, pickle.HIGHEST_PROTOCOL)
    with open(Path('test/testcase08.pickle'), 'rb') as f:
        expected = pickle.load(f)
    assert the_HMM == expected

    theta, sigma, alphabet, alignment = fetch_alignment(Path('test/testcase05.txt'))
    the_HMM = hmm.make_profile_HMM(theta=theta, sigma=sigma, alphabet=alphabet, alignment=alignment)
    # with open(Path('test/testcase05.pickle'), 'wb') as f:
    #     pickle.dump(the_HMM, f, pickle.HIGHEST_PROTOCOL)
    with open(Path('test/testcase05.pickle'), 'rb') as f:
        expected = pickle.load(f)
    assert the_HMM == expected

    theta, sigma, alphabet, alignment = fetch_alignment(Path('test/testcase10.txt'))
    the_HMM = hmm.make_profile_HMM(theta=theta, sigma=sigma, alphabet=alphabet, alignment=alignment)
    # with open(Path('test/testcase10.pickle'), 'wb') as f:
    #    pickle.dump(the_HMM, f, pickle.HIGHEST_PROTOCOL)
    with open(Path('test/testcase10.pickle'), 'rb') as f:
        expected = pickle.load(f)
    assert the_HMM == expected

    # TODO for Rosalind
    # Save in file and use tab separated
    # Replace theta wiht 1-theta
    theta, sigma, alphabet, alignment = fetch_alignment(Path('test/testcase09.txt'))
    the_HMM = hmm.make_profile_HMM(theta=theta, sigma=sigma, alphabet=alphabet, alignment=alignment)
    # print()
    # ugly_print_matrices(the_HMM.transition, the_HMM.emission, the_HMM.transition.keys(), the_HMM.alphabet)


def test_profile_alignment():
    text, theta, sigma, alphabet, alignment = fetch_profile_alignment(Path('test/testcase11.txt'))
    path, score = hmm.align(emissions=text, theta=theta, sigma=sigma, alphabet=alphabet, alignment=alignment)
    assert path == [('M', 1), ('D', 2), ('D', 3), ('M', 4), ('M', 5), ('I', 5), ('M', 6), ('M', 7), ('M', 8)]
    assert score == -8.893769303174695

    text, theta, sigma, alphabet, alignment = fetch_profile_alignment(Path('test/testcase12.txt'))
    path, score = hmm.align(emissions=text, theta=theta, sigma=sigma, alphabet=alphabet, alignment=alignment)
    assert path == [('M', 1), ('M', 2), ('M', 3), ('M', 4), ('M', 5), ('M', 6), ('M', 7), ('M', 8), ('M', 9), ('D', 10),
                    ('M', 11), ('M', 12), ('I', 12), ('I', 12), ('M', 13), ('M', 14), ('M', 15), ('M', 16), ('M', 17),
                    ('M', 18), ('D', 19), ('M', 20), ('M', 21), ('M', 22), ('M', 23), ('I', 23), ('M', 24), ('M', 25),
                    ('M', 26), ('I', 26), ('I', 26), ('M', 27), ('D', 28), ('M', 29), ('M', 30), ('M', 31), ('I', 31),
                    ('I', 31), ('D', 32), ('M', 33), ('M', 34), ('I', 34), ('M', 35), ('M', 36), ('M', 37), ('M', 38),
                    ('I', 38), ('M', 39), ('M', 40), ('D', 41), ('M', 42), ('M', 43), ('D', 44), ('I', 44), ('I', 44)]
    assert score == -90.94050212069119

    text, theta, sigma, alphabet, alignment = fetch_profile_alignment(Path('test/testcase13.txt'))
    path, score = hmm.align(emissions=text, theta=theta, sigma=sigma, alphabet=alphabet, alignment=alignment)
    assert path == [('M', 1), ('M', 2), ('M', 3), ('I', 3), ('M', 4), ('M', 5), ('M', 6), ('M', 7), ('M', 8), ('M', 9),
                    ('D', 10), ('M', 11), ('M', 12), ('M', 13), ('M', 14), ('M', 15), ('I', 15), ('I', 15), ('M', 16),
                    ('M', 17), ('I', 17), ('M', 18), ('M', 19), ('M', 20), ('M', 21), ('M', 22), ('M', 23), ('M', 24),
                    ('M', 25), ('M', 26), ('I', 26), ('M', 27), ('I', 27), ('M', 28), ('M', 29), ('M', 30), ('I', 30),
                    ('M', 31), ('M', 32), ('M', 33), ('M', 34), ('I', 34), ('M', 35), ('M', 36), ('M', 37), ('M', 38),
                    ('M', 39), ('M', 40), ('M', 41), ('M', 42)]
    assert score == -63.01881096946894

    text, theta, sigma, alphabet, alignment = fetch_profile_alignment(Path('test/testcase14.txt'))
    path, score = hmm.align(emissions=text, theta=theta, sigma=sigma, alphabet=alphabet, alignment=alignment)
    assert path == [('M', 1), ('M', 2), ('M', 3), ('M', 4), ('M', 5), ('M', 6), ('I', 6), ('M', 7), ('M', 8), ('I', 8),
                    ('I', 8), ('I', 8), ('I', 8), ('M', 9), ('M', 10), ('I', 10), ('M', 11), ('I', 11), ('M', 12),
                    ('M', 13), ('M', 14), ('D', 15), ('I', 15), ('I', 15), ('M', 16), ('I', 16), ('M', 17), ('M', 18),
                    ('M', 19), ('I', 19), ('M', 20), ('M', 21), ('I', 21), ('M', 22), ('I', 22), ('I', 22), ('I', 22),
                    ('I', 22), ('I', 22), ('M', 23), ('M', 24), ('I', 24), ('I', 24), ('I', 24), ('I', 24), ('I', 24),
                    ('M', 25), ('M', 26), ('D', 27), ('I', 27), ('M', 28), ('D', 29)]
    assert score == -94.51230211004193
