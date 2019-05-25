# import matplotlib.pyplot as plt
from pathlib import Path
import pickle
import networkx as nx
import networkx.algorithms.isomorphism as iso
from numpy import isclose
# import networkx.drawing.nx_pylab as nxp
from stepik_hmm import fetch_hmm, fetch_alignment, ugly_print_matrices, fetch_profile_alignment, pretty_alignment_print, \
    fetch_hmm_path, fetch_hmm_outcome, pretty_print_path, fetch_parameter_estimation, ugly_print_matrix, \
    fetch_viterbi_learning, pretty_print_cond_prob
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


def test_viterbi_learning():
    n_iterations, emissions, alphabet, states, transision_matrix, emission_matrix = fetch_viterbi_learning(
        Path('test/testcase35.txt'))
    trans_m, em_m = hmm.viterbi_learning(n_iterations, emissions, transision_matrix, emission_matrix)
    # with open(Path('test/testcase35.pickle'), 'wb') as f:
    #     pickle.dump((trans_m, em_m), f, pickle.HIGHEST_PROTOCOL)
    with open(Path('test/testcase35.pickle'), 'rb') as f:
        expected = pickle.load(f)
    assert trans_m, em_m == expected

    n_iterations, emissions, alphabet, states, transision_matrix, emission_matrix = fetch_viterbi_learning(
        Path('test/testcase36.txt'))
    trans_m, em_m = hmm.viterbi_learning(n_iterations, emissions, transision_matrix, emission_matrix)
    # with open(Path('test/testcase36.pickle'), 'wb') as f:
    #    pickle.dump((trans_m, em_m), f, pickle.HIGHEST_PROTOCOL)
    with open(Path('test/testcase36.pickle'), 'rb') as f:
        expected = pickle.load(f)
    assert trans_m, em_m == expected

    n_iterations, emissions, alphabet, states, transision_matrix, emission_matrix = fetch_viterbi_learning(
        Path('test/testcase37.txt'))
    trans_m, em_m = hmm.viterbi_learning(n_iterations, emissions, transision_matrix, emission_matrix)
    # with open(Path('test/testcase37.pickle'), 'wb') as f:
    #     pickle.dump((trans_m, em_m), f, pickle.HIGHEST_PROTOCOL)
    with open(Path('test/testcase37.pickle'), 'rb') as f:
        expected = pickle.load(f)
    assert trans_m, em_m == expected
    # print()
    # ugly_print_matrix(trans_m, states, states)
    # print('--------')
    # ugly_print_matrix(em_m, states, alphabet)


def test_soft_decoding():
    emissions, model = fetch_hmm(Path('test/testcase38.txt'))
    prob = hmm.soft_decoding(emissions, model)
    assert prob == {10: {'A': 0.8167149234719551, 'B': 0.18328507652804493},
                    9: {'A': 0.8737316291602154, 'B': 0.12626837083978454},
                    8: {'A': 0.9640271511628047, 'B': 0.03597284883719514},
                    7: {'A': 0.9153809707336851, 'B': 0.08461902926631495},
                    6: {'A': 0.9891304185958109, 'B': 0.01086958140418904},
                    5: {'A': 0.9956659237728355, 'B': 0.004334076227164441},
                    4: {'A': 0.9936237458533239, 'B': 0.006376254146675927},
                    3: {'A': 0.964681556054572, 'B': 0.035318443945427914},
                    2: {'A': 0.6491713244181421, 'B': 0.35082867558185776},
                    1: {'A': 0.543808916049139, 'B': 0.456191083950861}}

    emissions, model = fetch_hmm(Path('test/testcase39.txt'))
    prob = hmm.soft_decoding(emissions, model)
    assert prob == {
        10: {'A': 0.223289369266927, 'B': 0.13556464995539308, 'C': 0.16567729338903842, 'D': 0.47552591055508003},
        9: {'A': 0.3695792530135062, 'B': 0.05784449613333948, 'C': 0.19779995810755133, 'D': 0.37483351591204145},
        8: {'A': 0.5087855851743058, 'B': 0.2689760725735352, 'C': 0.1975573577850791, 'D': 0.024738207633518347},
        7: {'A': 0.27886945530939805, 'B': 0.15359335093539414, 'C': 0.21386244022381676, 'D': 0.3537319766978296},
        6: {'A': 0.3031460771869307, 'B': 0.2213051307048867, 'C': 0.23388690769439863, 'D': 0.2417191075802224},
        5: {'A': 0.4413839668701327, 'B': 0.2628282870020594, 'C': 0.26734516932400904, 'D': 0.028499799970237175},
        4: {'A': 0.12970777817517634, 'B': 0.5359490137873554, 'C': 0.15417972574071118, 'D': 0.18022070546319552},
        3: {'A': 0.15108839833338034, 'B': 0.12510807820668837, 'C': 0.15534609034541605, 'D': 0.5685146562809537},
        2: {'A': 0.3647848535970413, 'B': 0.05300301759771121, 'C': 0.1909222870924323, 'D': 0.3913470648792535},
        1: {'A': 0.5002916546988827, 'B': 0.2113798200881886, 'C': 0.26619795818566405, 'D': 0.022187790193703014}}

    emissions, model = fetch_hmm(Path('test/testcase40.txt'))
    prob = hmm.soft_decoding(emissions, model)
    assert prob == {10: {'A': 0.27251473363412476, 'B': 0.06877055178886632, 'C': 0.6587147145770089},
                    9: {'A': 0.24644905928986852, 'B': 0.06775208504530535, 'C': 0.6857988556648262},
                    8: {'A': 0.8564177157596157, 'B': 0.0485827763793196, 'C': 0.09499950786106474},
                    7: {'A': 0.49947992359126636, 'B': 0.04092795555338618, 'C': 0.4595921208553475},
                    6: {'A': 0.24692594463437936, 'B': 0.0630143174963245, 'C': 0.6900597378692964},
                    5: {'A': 0.5194205792262436, 'B': 0.09938843946662584, 'C': 0.38119098130713064},
                    4: {'A': 0.23676809114626354, 'B': 0.07827750667115782, 'C': 0.6849544021825789},
                    3: {'A': 0.48129185400062685, 'B': 0.11777609952517239, 'C': 0.400932046474201},
                    2: {'A': 0.3950922233250534, 'B': 0.1268233812850449, 'C': 0.4780843953899019},
                    1: {'A': 0.7006610552681448, 'B': 0.2189828497060014, 'C': 0.08035609502585392}}

    emissions, model = fetch_hmm(Path('test/testcase41.txt'))
    prob = hmm.soft_decoding(emissions, model)
    assert prob == {10: {'A': 0.60230736750807, 'B': 0.39769263249192993},
                    9: {'A': 0.5256510452843234, 'B': 0.47434895471567656},
                    8: {'A': 0.514915220951736, 'B': 0.48508477904826386},
                    7: {'A': 0.5672049370701465, 'B': 0.43279506292985337},
                    6: {'A': 0.4982841969762023, 'B': 0.5017158030237975},
                    5: {'A': 0.5181739029263772, 'B': 0.4818260970736225},
                    4: {'A': 0.40346161694074234, 'B': 0.5965383830592574},
                    3: {'A': 0.3332677918449609, 'B': 0.666732208155039},
                    2: {'A': 0.315520101298843, 'B': 0.6844798987011569},
                    1: {'A': 0.3438907578307824, 'B': 0.6561092421692175}}
