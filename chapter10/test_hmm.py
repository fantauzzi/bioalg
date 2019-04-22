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
    assert path == [('M', 1, 1), ('D', 2, 1), ('D', 3, 1), ('M', 4, 2), ('M', 5, 3), ('I', 5, 4), ('M', 6, 5),
                    ('M', 7, 6), ('M', 8, 7)]
    assert score == -8.893769303174695

    text, theta, sigma, alphabet, alignment = fetch_profile_alignment(Path('test/testcase12.txt'))
    path, score = hmm.align(emissions=text, theta=theta, sigma=sigma, alphabet=alphabet, alignment=alignment)
    assert path == [('M', 1, 1), ('M', 2, 2), ('M', 3, 3), ('M', 4, 4), ('M', 5, 5), ('M', 6, 6), ('M', 7, 7),
                    ('M', 8, 8), ('M', 9, 9), ('D', 10, 9), ('M', 11, 10), ('M', 12, 11), ('I', 12, 12), ('I', 12, 13),
                    ('M', 13, 14), ('M', 14, 15), ('M', 15, 16), ('M', 16, 17), ('M', 17, 18), ('M', 18, 19),
                    ('D', 19, 19), ('M', 20, 20), ('M', 21, 21), ('M', 22, 22), ('M', 23, 23), ('I', 23, 24),
                    ('M', 24, 25), ('M', 25, 26), ('M', 26, 27), ('I', 26, 28), ('I', 26, 29), ('M', 27, 30),
                    ('D', 28, 30), ('M', 29, 31), ('M', 30, 32), ('M', 31, 33), ('I', 31, 34), ('I', 31, 35),
                    ('D', 32, 35), ('M', 33, 36), ('M', 34, 37), ('I', 34, 38), ('M', 35, 39), ('M', 36, 40),
                    ('M', 37, 41), ('M', 38, 42), ('I', 38, 43), ('M', 39, 44), ('M', 40, 45), ('D', 41, 45),
                    ('M', 42, 46), ('M', 43, 47), ('D', 44, 47), ('I', 44, 48), ('I', 44, 49)]
    assert score == -90.94050212069119

    text, theta, sigma, alphabet, alignment = fetch_profile_alignment(Path('test/testcase13.txt'))
    path, score = hmm.align(emissions=text, theta=theta, sigma=sigma, alphabet=alphabet, alignment=alignment)
    assert path == [('M', 1, 1), ('M', 2, 2), ('M', 3, 3), ('I', 3, 4), ('M', 4, 5), ('M', 5, 6), ('M', 6, 7),
                    ('M', 7, 8), ('M', 8, 9), ('M', 9, 10), ('D', 10, 10), ('M', 11, 11), ('M', 12, 12), ('M', 13, 13),
                    ('M', 14, 14), ('M', 15, 15), ('I', 15, 16), ('I', 15, 17), ('M', 16, 18), ('M', 17, 19),
                    ('I', 17, 20), ('M', 18, 21), ('M', 19, 22), ('M', 20, 23), ('M', 21, 24), ('M', 22, 25),
                    ('M', 23, 26), ('M', 24, 27), ('M', 25, 28), ('M', 26, 29), ('I', 26, 30), ('M', 27, 31),
                    ('I', 27, 32), ('M', 28, 33), ('M', 29, 34), ('M', 30, 35), ('I', 30, 36), ('M', 31, 37),
                    ('M', 32, 38), ('M', 33, 39), ('M', 34, 40), ('I', 34, 41), ('M', 35, 42), ('M', 36, 43),
                    ('M', 37, 44), ('M', 38, 45), ('M', 39, 46), ('M', 40, 47), ('M', 41, 48), ('M', 42, 49)]
    assert score == -63.01881096946894

    text, theta, sigma, alphabet, alignment = fetch_profile_alignment(Path('test/testcase14.txt'))
    path, score = hmm.align(emissions=text, theta=theta, sigma=sigma, alphabet=alphabet, alignment=alignment)
    assert path == [('M', 1, 1), ('M', 2, 2), ('M', 3, 3), ('M', 4, 4), ('M', 5, 5), ('M', 6, 6), ('I', 6, 7),
                    ('M', 7, 8), ('M', 8, 9), ('I', 8, 10), ('I', 8, 11), ('I', 8, 12), ('I', 8, 13), ('M', 9, 14),
                    ('M', 10, 15), ('I', 10, 16), ('M', 11, 17), ('I', 11, 18), ('M', 12, 19), ('M', 13, 20),
                    ('M', 14, 21), ('D', 15, 21), ('I', 15, 22), ('I', 15, 23), ('M', 16, 24), ('I', 16, 25),
                    ('M', 17, 26), ('M', 18, 27), ('M', 19, 28), ('I', 19, 29), ('M', 20, 30), ('M', 21, 31),
                    ('I', 21, 32), ('M', 22, 33), ('I', 22, 34), ('I', 22, 35), ('I', 22, 36), ('I', 22, 37),
                    ('I', 22, 38), ('M', 23, 39), ('M', 24, 40), ('I', 24, 41), ('I', 24, 42), ('I', 24, 43),
                    ('I', 24, 44), ('I', 24, 45), ('M', 25, 46), ('M', 26, 47), ('D', 27, 47), ('I', 27, 48),
                    ('M', 28, 49), ('D', 29, 49)]
    assert score == -94.51230211004193
