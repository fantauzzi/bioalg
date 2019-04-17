import matplotlib.pyplot as plt
from pathlib import Path
import networkx as nx
import networkx.algorithms.isomorphism as iso
import networkx.drawing.nx_pylab as nxp
from stepik_hmm import fetch_hmm
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
