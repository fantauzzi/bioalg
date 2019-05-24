# import matplotlib.pyplot as plt
from pathlib import Path
import pickle
import networkx as nx
import networkx.algorithms.isomorphism as iso
from numpy import isclose
# import networkx.drawing.nx_pylab as nxp
from stepik_hmm import fetch_hmm, fetch_alignment, ugly_print_matrices, fetch_profile_alignment, pretty_alignment_print, \
    fetch_hmm_path, fetch_hmm_outcome, pretty_print_path, fetch_parameter_estimation, ugly_print_matrix
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

    emissions, model = fetch_hmm(Path('test/testcase21.txt'))
    path = hmm.viterbi(emissions=emissions, model=model)
    assert path == 'AAABBAAAAA'

    emissions, model = fetch_hmm(Path('test/testcase23.txt'))
    path = hmm.viterbi(emissions=emissions, model=model)
    assert path == 'CCCDABBBBBBBBBBBBBBBBBBBBBBCDACDACCCDABBBBBDACDACDABBBBBBBBBBBBBBBBBBBBBBBBBBBBBDADACCDADACCDADADADA'


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

    theta, sigma, alphabet, alignment = fetch_alignment(Path('test/testcase09.txt'))
    the_HMM = hmm.make_profile_HMM(theta=theta, sigma=sigma, alphabet=alphabet, alignment=alignment)
    # with open(Path('test/testcase09.pickle'), 'wb') as f:
    #    pickle.dump(the_HMM, f, pickle.HIGHEST_PROTOCOL)
    with open(Path('test/testcase09.pickle'), 'rb') as f:
        expected = pickle.load(f)
    assert the_HMM == expected

    theta, sigma, alphabet, alignment = fetch_alignment(Path('test/testcase28.txt'))
    the_HMM = hmm.make_profile_HMM(theta=theta, sigma=sigma, alphabet=alphabet, alignment=alignment)
    # with open(Path('test/testcase28.pickle'), 'wb') as f:
    #    pickle.dump(the_HMM, f, pickle.HIGHEST_PROTOCOL)
    with open(Path('test/testcase28.pickle'), 'rb') as f:
        expected = pickle.load(f)
    assert the_HMM == expected

    theta, sigma, alphabet, alignment = fetch_alignment(Path('test/testcase29.txt'))
    the_HMM = hmm.make_profile_HMM(theta=theta, sigma=sigma, alphabet=alphabet, alignment=alignment)
    print()
    # ugly_print_matrices(the_HMM.transition, the_HMM.emission, the_HMM.transition.keys(), the_HMM.alphabet)
    # with open(Path('test/testcase29.pickle'), 'wb') as f:
    #     pickle.dump(the_HMM, f, pickle.HIGHEST_PROTOCOL)
    with open(Path('test/testcase29.pickle'), 'rb') as f:
        expected = pickle.load(f)
    assert the_HMM == expected


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

    text, theta, sigma, alphabet, alignment = fetch_profile_alignment(Path('test/testcase30.txt'))
    path, score = hmm.align(emissions=text, theta=theta, sigma=sigma, alphabet=alphabet, alignment=alignment)
    assert path == [('M', 1), ('I', 1), ('M', 2), ('M', 3), ('M', 4), ('M', 5), ('M', 6), ('M', 7), ('M', 8), ('I', 8),
                    ('M', 9), ('M', 10), ('I', 10), ('M', 11), ('M', 12), ('M', 13), ('M', 14), ('I', 14), ('I', 14),
                    ('M', 15), ('I', 15), ('M', 16), ('M', 17), ('I', 17), ('I', 17), ('M', 18), ('M', 19), ('M', 20),
                    ('M', 21), ('M', 22), ('I', 22), ('M', 23), ('M', 24), ('M', 25), ('M', 26), ('M', 27), ('D', 28),
                    ('M', 29), ('I', 29), ('I', 29), ('D', 30), ('M', 31), ('I', 31), ('M', 32), ('M', 33), ('I', 33),
                    ('M', 34), ('M', 35), ('M', 36), ('M', 37), ('M', 38)]
    assert score == -74.42156226598739


def test_hidden_path_prob():
    path, _, transitions = fetch_hmm_path(Path('test/testcase15.txt'))
    prob = hmm.hidden_path_prob(path, transitions)
    assert isclose(prob, 0.000384928691755)

    path, _, transitions = fetch_hmm_path(Path('test/testcase16.txt'))
    prob = hmm.hidden_path_prob(path, transitions)
    assert isclose(prob, 3.26233331904e-21)

    path, _, transitions = fetch_hmm_path(Path('test/testcase17.txt'))
    prob = hmm.hidden_path_prob(path, transitions)
    assert isclose(prob, 1.5860533198043927e-19)


def test_outcome_prob():
    emissions, path, emission_matrix = fetch_hmm_outcome(Path('test/testcase18.txt'))
    prob = hmm.outcome_prob(emissions, path, emission_matrix)
    assert isclose(prob, 3.59748954746e-06)

    emissions, path, emission_matrix = fetch_hmm_outcome(Path('test/testcase19.txt'))
    prob = hmm.outcome_prob(emissions, path, emission_matrix)
    assert isclose(prob, 3.42316482177e-35)

    emissions, path, emission_matrix = fetch_hmm_outcome(Path('test/testcase20.txt'))
    prob = hmm.outcome_prob(emissions, path, emission_matrix)
    assert isclose(prob, 3.660029947725436e-27)


def test_outcome_likelyhood():
    emissions, model = fetch_hmm(Path('test/testcase24.txt'))
    prob = hmm.outcome_likelyhood(emissions=emissions, model=model)
    assert isclose(prob, 1.1005510319694847e-06)

    emissions, model = fetch_hmm(Path('test/testcase25.txt'))
    prob = hmm.outcome_likelyhood(emissions=emissions, model=model)
    assert isclose(prob, 4.08210708381e-55)

    emissions, model = fetch_hmm(Path('test/testcase26.txt'))
    prob = hmm.outcome_likelyhood(emissions=emissions, model=model)
    assert isclose(prob, 4.577544141262532e-49)

    emissions, model = fetch_hmm(Path('test/testcase27.txt'))
    prob = hmm.outcome_likelyhood(emissions=emissions, model=model)
    assert isclose(prob, 9.461855360076486e-51)


def test_hmm_parameter_estimation():
    emissions, alphabet, path, states = fetch_parameter_estimation(Path('test/testcase31.txt'))
    trans_p, em_p = hmm.hmm_parameter_estimation(emissions, alphabet, path, states)
    # with open(Path('test/testcase31.pickle'), 'wb') as f:
    #    pickle.dump((trans_p, em_p), f, pickle.HIGHEST_PROTOCOL)
    with open(Path('test/testcase31.pickle'), 'rb') as f:
        expected = pickle.load(f)
    assert trans_p, em_p == expected

    emissions, alphabet, path, states = fetch_parameter_estimation(Path('test/testcase32.txt'))
    trans_p, em_p = hmm.hmm_parameter_estimation(emissions, alphabet, path, states)
    # with open(Path('test/testcase32.pickle'), 'wb') as f:
    #     pickle.dump((trans_p, em_p), f, pickle.HIGHEST_PROTOCOL)
    with open(Path('test/testcase32.pickle'), 'rb') as f:
        expected = pickle.load(f)
    assert trans_p, em_p == expected

    emissions, alphabet, path, states = fetch_parameter_estimation(Path('test/testcase33.txt'))
    trans_p, em_p = hmm.hmm_parameter_estimation(emissions, alphabet, path, states)
    # with open(Path('test/testcase33.pickle'), 'wb') as f:
    #     pickle.dump((trans_p, em_p), f, pickle.HIGHEST_PROTOCOL)
    with open(Path('test/testcase33.pickle'), 'rb') as f:
        expected = pickle.load(f)
    assert trans_p, em_p == expected

    emissions, alphabet, path, states = fetch_parameter_estimation(Path('test/testcase34.txt'))
    trans_p, em_p = hmm.hmm_parameter_estimation(emissions, alphabet, path, states)
    # print()
    # ugly_print_matrix(trans_p, row_labels=states, col_labels=states)
    # print('--------')
    # ugly_print_matrix(em_p, row_labels=states, col_labels=alphabet)
    with open(Path('test/testcase34.pickle'), 'wb') as f:
        pickle.dump((trans_p, em_p), f, pickle.HIGHEST_PROTOCOL)
    # with open(Path('test/testcase34.pickle'), 'rb') as f:
    #    expected = pickle.load(f)
    assert trans_p, em_p == expected
