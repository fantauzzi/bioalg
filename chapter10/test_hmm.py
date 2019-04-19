# import matplotlib.pyplot as plt
from pathlib import Path
import pickle
import networkx as nx
import networkx.algorithms.isomorphism as iso
# import networkx.drawing.nx_pylab as nxp
from stepik_hmm import fetch_hmm, fetch_alignment, ugly_print_matrices
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
