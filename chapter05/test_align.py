import align
import numpy as np
from pathlib import Path


def read_weighted_graph(file_name):
    with open(file_name) as input_file:
        lines = input_file.readlines()
    adj = {}
    for line in lines:
        left, right = line.rstrip().split('->')
        vertex = left
        vertex2, weight_s = right.split(':')
        weight = int(weight_s)
        current = adj.get(vertex, [])
        current.append((vertex2, weight))
        adj[vertex] = current
    return adj


def test_dp_change():
    change = align.dp_change(40, [50, 25, 20, 10, 5, 1])
    assert change == 2

    change2 = align.dp_change(19415, [18, 16, 7, 5, 3, 1])
    assert change2 == 1080

    change3 = align.dp_change(16042, [16, 13, 11, 8, 7, 5, 3, 1])
    assert change3 == 1003


def test_manhattan_tourist():
    down = np.array([[1, 0, 2, 4, 3],
                     [4, 6, 5, 2, 1],
                     [4, 4, 5, 2, 1],
                     [5, 6, 8, 5, 3]])

    right = np.array([[3, 2, 4, 0],
                      [3, 2, 4, 2],
                      [0, 7, 3, 3],
                      [3, 3, 0, 2],
                      [1, 3, 2, 2]])

    dist = align.manhattan_tourist(down=down, right=right)
    assert dist == 34

    down2 = np.array([[4, 4, 2, 1, 3, 1, 0, 0, 1],
                      [4, 3, 0, 2, 0, 4, 3, 4, 4],
                      [2, 3, 3, 1, 2, 1, 2, 2, 0],
                      [3, 0, 3, 3, 2, 1, 1, 3, 4],
                      [2, 3, 1, 2, 2, 0, 2, 3, 2],
                      [2, 2, 2, 0, 4, 2, 1, 0, 3],
                      [4, 3, 1, 1, 0, 1, 1, 4, 2],
                      [0, 0, 2, 2, 2, 1, 2, 4, 2],
                      [4, 3, 0, 3, 1, 3, 2, 3, 1],
                      [1, 4, 1, 0, 3, 4, 1, 2, 1],
                      [4, 4, 0, 4, 1, 4, 3, 1, 2],
                      [4, 1, 2, 3, 1, 3, 3, 3, 0],
                      [3, 1, 0, 2, 2, 0, 4, 4, 0],
                      [2, 0, 1, 0, 0, 3, 1, 1, 1],
                      [0, 1, 3, 2, 2, 2, 1, 2, 1],
                      [0, 2, 0, 3, 1, 2, 2, 4, 2],
                      [2, 0, 4, 1, 3, 3, 2, 4, 0],
                      [2, 3, 1, 3, 4, 2, 1, 4, 4]], dtype=np.int)

    right2 = np.array([[3, 3, 1, 1, 3, 4, 4, 4],
                       [4, 0, 3, 1, 0, 3, 4, 4],
                       [2, 2, 2, 3, 3, 1, 1, 4],
                       [1, 3, 1, 4, 4, 2, 0, 1],
                       [0, 2, 0, 3, 3, 3, 1, 0],
                       [3, 2, 0, 4, 1, 4, 4, 3],
                       [3, 0, 1, 1, 0, 3, 3, 0],
                       [3, 1, 1, 0, 2, 3, 4, 0],
                       [2, 4, 2, 1, 1, 3, 1, 2],
                       [1, 0, 4, 3, 0, 3, 3, 0],
                       [2, 3, 2, 4, 4, 3, 3, 0],
                       [3, 1, 2, 0, 3, 4, 3, 2],
                       [0, 0, 4, 1, 4, 0, 3, 4],
                       [3, 2, 3, 2, 0, 1, 2, 1],
                       [4, 3, 3, 2, 0, 1, 1, 2],
                       [0, 0, 4, 1, 2, 4, 0, 3],
                       [3, 4, 0, 1, 2, 3, 0, 1],
                       [4, 0, 2, 4, 2, 2, 4, 0],
                       [4, 3, 4, 2, 2, 3, 2, 3]], dtype=np.int)

    dist2 = align.manhattan_tourist(down=down2, right=right2)
    assert dist2 == 80


def test_dag_longest_path():
    adj = {'0': [('1', 10), ('2', 4)],
           '2': [('3', 2), ('5', 1)],
           '1': [('4', 1)],
           '3': [('4', 3)],
           '5': [('3', 4)]}

    topo_order = align.topological_ordering(adj)
    assert topo_order == ['0', '2', '5', '3', '1', '4']
    path, distances = align.dag_longest_path(adj, '0', '4')
    assert distances[-1] == 12
    assert path == ['0', '2', '5', '3', '4']

    adj2 = read_weighted_graph(Path('test/testcase01.txt'))
    topo_order2 = align.topological_ordering(adj2)
    assert topo_order2 == ['0', '2', '3', '1', '4']
    path2, distance2 = align.dag_longest_path(adj2, '0', '4')
    assert distance2[-1] == 9
    assert path2 == ['0', '2', '3', '4']

    adj3 = read_weighted_graph(Path('test/testcase03.txt'))
    path3, distances3 = align.dag_longest_path(adj3, '5', '8')
    topo_order3 = align.topological_ordering(adj3)
    assert topo_order3 == ['4', '3', '2', '1', '5', '7', '6', '8']
    assert distances3[-1] == 5
    assert path3 == ['5', '6', '8']

    adj4 = read_weighted_graph(Path('test/testcase02.txt'))
    topo_order4 = align.topological_ordering(adj4)
    assert topo_order4 == ['9', '10', '3', '19', '2', '6', '4', '0', '7', '13', '8', '1', '11', '5', '14', '15', '23',
                           '26', '21', '17', '27', '12', '25', '28', '20', '24', '16', '18']
    path4, distances4 = align.dag_longest_path(adj4, '5', '20')
    assert path4 == ['5', '14', '15', '20']
    assert distances4[-1] == 36

    adj5 = read_weighted_graph(Path('test/testcase04.txt'))
    topo_order5 = align.topological_ordering(adj5)
    assert topo_order5 == ['2', '9', '5', '8', '3', '6', '19', '7', '0', '4', '11', '1', '12', '18', '13', '21', '22',
                           '26', '10', '14', '15', '16', '17', '20', '23', '24', '27', '25']
    path5, distances5 = align.dag_longest_path(adj5, '12', '27')
    assert distances5[-1] == 94
    assert path5 == ['12', '13', '21', '22', '24', '27']


def test_longest_common_string():
    string1 = 'ATGTTATA'
    string2 = 'ATCGTCC'
    res = align.longest_common_string(string1, string2)
    assert res == 'ATGT'

    string1 = 'AACCTTGG'
    string2 = 'ACACTGTGA'
    res = align.longest_common_string(string1, string2)
    assert res == 'ACCTTG'

    string1 = 'AGCAGTTCCCTGATTGTTTAGTATTTGACTCCGTAGTTGAGCCTATATCGTAATTCTGCCAAGGAA'
    string2 = 'ATTATAATCCCCCGGACAAGAACCTAGTGGGCGCGTGGGACGGGACAAGAGCGACGTTCCGTAGGTGTGAGAGCCGTACTCAATTTTGGTTATTACCAT'
    res = align.longest_common_string(string1, string2)
    assert res == 'AATTCCCGATTGTAGAGACTCCGTAGTTGAGCCTACTAATTTGAAA'


def test_best_protein_alignment():
    alphabet, scoring_matrix = align.get_blosum62()

    score, alignment = align.best_alignment('PLEASANTLY', 'MEANLY', scoring_matrix, alphabet, sigma=5, local=False)

    assert score == 8
    assert alignment == ('PLEASANTLY', '-ME--AN-LY')

    score, alignment = align.best_alignment(
        'KHLGRRPTYGFPFWYMVWDFQCQDDKEQKFFCKPRHVPCTWLGCEVTDEMWMDLHVEVQPQFCLVRQEFWHIFPPFSSIYWMYFDPSDVNRIMHDD',
        'KPTYGFPFWYMDWDFQCQDEWKKEIRFCKEQKFFCKPRHVPCWWLGCEVMDLHTQRYFWH',
        scoring_matrix,
        alphabet,
        sigma=5,
        local=False)

    assert score == 43
    assert alignment[
               0] == 'KHLGRRPTYGFPFWYMVWDFQCQD----D----KEQKFFCKPRHVPCTWLGCEVTDEMWMDLHVEVQPQFCLVRQEFWHIFPPFSSIYWMYFDPSDVNRIMHDD'

    assert alignment[
               1] == 'K-----PTYGFPFWYMDWDFQCQDEWKKEIRFCKEQKFFCKPRHVPCWWLGCEV-----MDLH--T--Q----R------Y--F----W------------H--'


def test_local_protein_alignment():
    alphabet, scoring_matrix = align.get_pam250()

    score, alignment = align.best_alignment('MEANLY',
                                            'PENALTY',
                                            scoring_matrix=scoring_matrix,
                                            alphabet=alphabet,
                                            sigma=5,
                                            local=True)

    assert score == 15
    assert alignment == ('EANL-Y', 'ENALTY')

    s1 = 'AMTAFRYRQGNPRYVKHFAYEIRLSHIWLLTQMPWEFVMGIKMPEDVFQHWRVYSVCTAEPMRSDETYPCELFTVFDDIFTAEPVVCSCFYDDPM'
    s2 = 'WQEKAVDGTVPSRHQYREKEDRQGNEIGKEFRRGPQVCEYSCNSHSCGWMPIFCIVCMSYVAFYCGLEYPMSRKTAKSQFIEWCDWFCFNHEFIPWVLRRYVVYDKIRYNYSYRNSASMEFV'

    score, alignment = align.best_alignment(s1,
                                            s2,
                                            scoring_matrix=scoring_matrix,
                                            alphabet=alphabet,
                                            sigma=5,
                                            local=True)

    assert score == 56
    assert alignment[0] == 'KHFAYEIRLSHIWLLTQMPWEFVMGIKMPE-DVFQH---W-RVYSVCTAEPMRSDETYPCEL-FTVFDDIFTAEPVV-CS--CFYDD'
    assert alignment[1] == 'RH-QY--REKEDRQGNEIGKEFRRGPQVCEYSCNSHSCGWMPIF--CIV-CMSYVAFY-CGLEYPMSRKTAKSQFIEWCDWFCFNHE'

    s1 = 'LLQWKRYAVWQFNHLQFPVHAAAAAAAAAADVAGTQCA'
    s2 = 'ANQCQYLSDGRCIRINHPFYVNSLVAAAAAAAAAAAKGAAAAAAAAAAGPPWVYIPAQWTGYVETRIGPT'

    score, alignment = align.best_alignment(s1,
                                            s2,
                                            scoring_matrix=scoring_matrix,
                                            alphabet=alphabet,
                                            sigma=5,
                                            local=True)

    assert score == 32
    assert alignment[0] == 'RYAVWQFNHLQFPVHAAAAAAAAAADVAG'
    assert alignment[1] == 'RC-IR-INH-PFYVNSLVAAAAAAAAAAA'


def test_edit_distance():
    dist = align.edit_distance('PLEASANTLY', 'MEANLY')
    assert dist == 5

    s1 = 'GGACRNQMSEVNMWGCWWASVWVSWCEYIMPSGWRRMKDRHMWHWSVHQQSSPCAKSICFHETKNQWNQDACGPKVTQHECMRRRLVIAVKEE'
    s2 = 'GMWGFVQVSTQSRFRHMWHWSVHQQSSECAKSICHHEWKNQWNQDACGPKVTQHECMANMPMHKCNNWFWRLVIAVKEEKVRETKMLDLIHRHWLVLNQGRMNEHNVTLRKSPCVKRIMHKWKSRTTFHR'
    dist = align.edit_distance(s1, s2)
    assert dist == 97

    dist = align.edit_distance('AC', 'AC')
    assert dist == 0

    dist = align.edit_distance('AT', 'G')
    assert dist == 2

    dist = align.edit_distance('CAGACCGAGTTAG', 'CGG')
    assert dist == 10

    dist = align.edit_distance('CGT', 'CAGACGGTGACG')
    assert dist == 9


def test_fit_aligh():
    long = 'CCAT'
    short = 'AT'
    score, (aligned1, aligned2) = align.fit_align(long, short)
    assert score == 2
    assert aligned1 == aligned2 == 'AT'

    long = 'GTAGGCTTAAGGTTA'
    short = 'TAGATA'
    score, (aligned1, aligned2) = align.fit_align(long, short)
    assert score == 2
    assert aligned1 == 'TAAGGTTA'
    assert aligned2 == 'T-A-GATA'

    long = 'CAATCACCCCAATCCCTCAATCCTGGCCCCACGCATAGGCTAATGCCAATCGCGGCCAGGGTATAACCGCCATAACTGTGGGTCAGAAGGGATAAGTTCCACAATCCTATTTTCCTCGAGGCGCTTCGATGCGTTAACGCGTACACTCTGTCGGCCAACCGTGTGGGAGCCGAATTGGCTGGGCTGTTGAACATTCTATCAGTAGATAAACGAAGGTACATCCGAGGTTGTCGATCGACCGCGGGGTCGTAGCGCGTGCATGTTCCTTTCAGGCCCACATACTCCGGAACGGTTCATATCACGACTATTCTTGCACAATCGGACAACGGTGTACCATGGTGGACACCGTAGGAGACCAATACTGCGTAAATCATAAGCATTGGAGAGTGGACTGCTAGCGAGGCTCACCATGGAGTCTCGGTCGGCATCTCCTGACTGCTGTTCCATCGCGTTTTTCTTTTACTCACGCAATAAATCAATACCCCCTAACACAGGCCTGCTCCAGCCTTATTAAGGCCATAGTAGCTCTACATGTAGACCGAACGGAAGCACAGTTTGGTAGAAATTCTTAATCGACTATGGTCCGTGCAGGCCAAAAAAGGAATAATCTTCGAATTCTCACGCCTTCATTAGGGCGCACATGGTGGGGTAAATCACTGCACTCTGTTCGCAGTTAAGCGTTGCAATCAATATCGGCAGAACTCGGAGTCCGTATAAAGCCGCCTCAGCGTGCACACGCCCGTGCGGCACGTCATTAGACGAGGATTCCGGGGGACTGGCCTGTTCGTAATCCACTAAAACAATGGTCCTACCATCTAAAACGCACCGTGTTCCCCTCTACGGGAACCCCCTAGAT'
    short = 'AGAGCGCAGAGAAGTCATTAGAACATGTAGCACATCGCTTATTAAGGGTCAATACCTAAAGGGCCTAACTATACGCCACACGGAACAGCTC'
    score, (aligned1, aligned2) = align.fit_align(long, short)
    assert score == 22
    assert aligned1 == 'AGGGCGCACATG--GTGGGGTA-AATCA-CT-GCAC-TCTG-TTCGCAGTTAAGCGTTGCAATCAATATCGGC-AGAACTCGGAGTCCGT-A-TAAAGCCGCCTCAGCGTGCACACGC-C'
    assert aligned2 == 'AGAGCGCAGA-GAAGT-CATTAGAA-CATGTAGCACATC-GCTT---A-TTAAG-G--G---TCAATA-C--CTA-AA---GG-G-CC-TAACTATA--CGCCACA-CG-GAACA-GCTC'
