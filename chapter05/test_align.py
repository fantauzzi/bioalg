import numpy as np
from pathlib import Path
from pandas import DataFrame
import align


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

    adj = {'a': [('b', 3), ('c', 6), ('d', 5)],
           'b': [('c', 2), ('f', 4)],
           'c': [('e', 4), ('f', 3), ('g', 7)],
           'd': [('e', 4), ('f', 5)],
           'e': [('g', 2)],
           'f': [('g', 1)]}
    topo_order = align.topological_ordering(adj)
    assert topo_order == ['a', 'd', 'b', 'c', 'f', 'e', 'g']
    path, distances = align.dag_longest_path(adj, 'a', 'g')
    assert path == ['a', 'c', 'g']
    assert distances[-1] == 13

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
    string1 = 'ACGATACGT'
    string2 = 'CCCATTAAGT'
    res1 = align.longest_common_string(string1, string2)
    assert res1 == 'CATAGT'

    string1 = 'ACGATACGT'
    string2 = 'GACTATAGAA'
    res2 = align.longest_common_string(string1, string2)
    assert res2 == 'ACATAG'

    string1 = 'CCCATTAAGT'
    string2 = 'GACTATAGAA'
    res3 = align.longest_common_string(string1, string2)
    assert res3 == 'ATAAG'

    string1 = 'CTCGAT'
    string2 = 'TACGTC'
    res4 = align.longest_common_string(string1, string2)
    assert res4 == 'TCGT'

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


    score, alignment = align.best_alignment('AMMY',
                                            'PMMTN',
                                            scoring_matrix=scoring_matrix,
                                            alphabet=alphabet,
                                            sigma=5,
                                            local=True)
    assert score == 13
    assert alignment[0] == 'AMYM'
    assert alignment[1] == 'PM-M'


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


def test_fit_align():
    long = 'GTTGGATTACGAATCGATATCTGTTTG'
    short = 'ACGTCG'
    score, (aligned1, aligned2) = align.fit_align(long, short)
    assert score == 4
    assert (aligned1, aligned2) == ('ACGAATCG', 'ACG--TCG')

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

    long = 'GACTTCTACTGACTTACTAGCAACGCTACACGTACGATCTATGCATACGTAGCTATCCTACGCACGTACGCGC'
    short = 'GGCACGCTAACGTACCCGATCTATAA'
    score, (aligned1, aligned2) = align.fit_align(long, short)
    assert score == 15
    assert aligned1 == 'AGCAACGCTACACGTA--CGATCTATGCA'
    assert aligned2 == 'GGC-ACGCTA-ACGTACCCGATCTAT-AA'


def test_overlap_alignment():
    s1 = 'PAWHEAE'
    s2 = 'HEAGAWGHEE'
    score, (aligned1, aligned2) = align.overlap_align(s1, s2)
    assert score == 1
    assert (aligned1, aligned2) == ('HEAE', 'HEA-')

    s1 = 'GTACGTTCAGCTAATCCTAGCATGGTATAGCTATAGCCGACGCCGCCTCCCGACATAGCTAA'
    s2 = 'GACTAGCTACTAGCCGACGCGCCTCCCGACATAGCTAATCATAATCAGACCATCGACGGACTAGCCATGC'
    score, (aligned1, aligned2) = align.overlap_align(s1, s2)
    assert score == 28
    assert aligned1 == 'GTA-TAGCTA-TAGCCGACGCCGCCTCCCGACATAGCTAA'
    assert aligned2 == 'G-ACTAGCTACTAGCCGACG-CGCCTCCCGACATAGCTAA'


def test_align_with_gap_penalties():
    alphabet, scoring_matrix = align.get_blosum62()

    ammino1 = 'PRTEINS'
    ammino2 = 'PRTWPSEIN'
    score, (aligned1, aligned2) = align.align_with_gap_penalties(ammino1,
                                                                 ammino2,
                                                                 alphabet,
                                                                 scoring_matrix,
                                                                 gap_open_penalty=11,
                                                                 gap_ext_penalty=1)
    assert score == 8
    assert aligned1 == 'PRT---EINS'
    assert aligned2 == 'PRTWPSEIN-'

    ammino1 = 'WC'
    ammino2 = 'WHC'
    score, (aligned1, aligned2) = align.align_with_gap_penalties(ammino1,
                                                                 ammino2,
                                                                 alphabet,
                                                                 scoring_matrix,
                                                                 gap_open_penalty=11,
                                                                 gap_ext_penalty=1)
    assert score == 9
    assert aligned1 == 'W-C'
    assert aligned2 == 'WHC'

    ammino1 = 'YHFDVPDCWAHRYWVENPQAIAQMEQICFNWFPSMMMKQPHVFKVDHHMSCRWLPIRGKKCSSCCTRMRVRTVWE'
    ammino2 = 'YHEDVAHEDAIAQMVNTFGFVWQICLNQFPSMMMKIYWIAVLSAHVADRKTWSKHMSCRWLPIISATCARMRVRTVWE'
    score, (aligned1, aligned2) = align.align_with_gap_penalties(ammino1,
                                                                 ammino2,
                                                                 alphabet,
                                                                 scoring_matrix,
                                                                 gap_open_penalty=11,
                                                                 gap_ext_penalty=1)
    assert score == 144
    assert aligned1 == 'YHFDVPDCWAHRYWVENPQAIAQME-------QICFNWFPSMMMK-------QPHVF---KVDHHMSCRWLPIRGKKCSSCCTRMRVRTVWE'
    assert aligned2 == 'YHEDV----AHE------DAIAQMVNTFGFVWQICLNQFPSMMMKIYWIAVLSAHVADRKTWSKHMSCRWLPI----ISATCARMRVRTVWE'

    ammino1 = 'GLWFLNEPSSPQYEFPWDYIHRQNKIEKIPVNMKVFNSNPPMPQTFTEQIIIACMADCAFCYGAGQRWCSYVP'
    ammino2 = 'GLWFLNEPSSPQYEAPSDYGTDPHLQNKIEKIPIYALAPNMWTVDPEKVFNSNPPMPQNIMFQEQIIIACMVWCSYVP'
    score, (aligned1, aligned2) = align.align_with_gap_penalties(ammino1,
                                                                 ammino2,
                                                                 alphabet,
                                                                 scoring_matrix,
                                                                 gap_open_penalty=11,
                                                                 gap_ext_penalty=1)
    assert score == 226
    assert aligned1 == 'GLWFLNEPSSPQYEFPWDY---IHRQNKIEKIPV-----NM------KVFNSNPPMPQT--FTEQIIIACMADCAFCYGAGQRWCSYVP'
    assert aligned2 == 'GLWFLNEPSSPQYEAPSDYGTDPHLQNKIEKIPIYALAPNMWTVDPEKVFNSNPPMPQNIMFQEQIIIACMV-----------WCSYVP'


def test_scoring_matrix_as_dict():
    alphabet, scoring_matrix = align.get_blosum62()
    as_dict = align.scoring_matrix_as_dict(alphabet, scoring_matrix)
    scoring_matrix = DataFrame(scoring_matrix, columns=alphabet, index=alphabet)
    for item1 in alphabet:
        for item2 in alphabet:
            assert as_dict[item1][item2] == scoring_matrix[item1][item2]


def test_middle_edge():
    alphabet, scoring_matrix = align.get_blosum62()

    res = align.middle_edge('FP', 'PF', alphabet, scoring_matrix, sigma=5)
    assert res == (-3, ((2, 1), (2, 2)))

    res = align.middle_edge('FPPF', 'FFPF', alphabet, scoring_matrix, sigma=5)
    assert res == (15, ((2, 2), (3, 3)))

    s1 = 'FPPPFPFFF'
    s2 = 'PPFPPPFFF'
    res = align.middle_edge(s1, s2, alphabet, scoring_matrix, sigma=5)
    assert res == (32, ((3, 4), (4, 5)))

    s1 = 'PLEASANTLY'
    s2 = 'MEASNLY'
    res = align.middle_edge(s1, s2, alphabet, scoring_matrix, sigma=5)
    assert res == (17, ((4, 3), (5, 4)))

    s1 = 'ACTTAATT'
    s2 = 'GAGCAATT'
    res = align.middle_edge(s1, s2, alphabet, scoring_matrix, sigma=5)
    assert res == (15, ((4, 4), (5, 5)))

    s1 = 'TWLNSACYGVNFRRLNPMNKTKWDCWTWVPMVMAAQYLCRIFIPVMDHWEFFGDWGLETWRLGIHDHVKIPNFRWSCELHI'
    s2 = 'LWFKFLQCIFQYFKDQQETNCIWTFSPFSEHICQRVCQVYWNWNTPSSRTSDPRELFANSTIHNNRCGEWR'
    res = align.middle_edge(s1, s2, alphabet, scoring_matrix, sigma=5)
    assert res == (14, ((40, 35), (41, 36)))

    s1 = 'TWLNSACYGVNFRRLNPMNKTKWDCWTWVPMVMAAQYLCRIFIPVMDHWEFFGDWGLETWRLGIHDHVKIPNFRWSCELHIREHGHHFKTRFLKHNQFTQCYGLMPDPQFHRSYDVACQWEVTMSQGLMRFHRQNQIEKQRDRTSTYCMMTIGPGFTSNGYDPFVTITITPVQEPVENWFTPGGSMGFMIISRYMQMFFYLTRFSDMTYLVGVHCENYVCWNNVAKFLNGNLQGIFDQGERAYHQFVTWHSYSQYSRCSVGRYACEQAMSRVNSKMTWHWPIRDQGHEHFSEQYLSEKRNPPCNPRIGNAGQHFYEIHRIAHRVAMCNWAPQGQHPGGPTPHDVETCLWLWSLCLKGSDRGYVDRPWMFLADQLGEANLTLITMFHGCTRGCLMWFMDWEECVCSYSVVNPRCHGSEQWSVQNLGWRTCDTLISLWEPECDKHNTPPCLHWEFEDHPSQLRPVMMCDKYVQSIPTDAKWAWTYSKDFVISHWLIWTPIKLEECVFPQINRLWGTACNQGSQKIVIQNVWLRPSSFFQERSKCSDSSCILNVGGSNVNITGKETRTHVPILHMHEIDLISTASSGMRHNLILPHGMLMLHMNWHHSTRAMNPYSSLKLIPWTFQVCETDDRDQNVATHVADPCHKGEDQEIRCCKGGVDHQWKGDRMWMMCMPDMNYVKQDQAPSGTCEGACENYPADKDKCYMIFTIVFDYRRCTKKVCIWISGFPVDAFNLISIANAGFFCCWLEPTELKWRRTFYLGKGTQGWMCTFPHRNIIPVIICAGFGRWVQGEVPFRPVAQISAHSSDRRQGHHPPGTNMCHDYGDQYPIKRVGMQVEEDDGASYCDCAADWKLADMYEADHLSIGVIDFTDWIYPKNGGIWSEIIKSHFHWYHWETPQNTVGAFNTIVGINGSDMCIYHGNTQWEFGWCWKWLNHGHMRNQGPCHLGILEGRISKFAQVTSWWWQTKHDKDWSIEPYGRHWGEAGRPYTYNYCWMRWAIVYNHGNVISVELVPFMDEYPGKCNKEDVQFELFSPMQA'
    s2 = 'LWFKFLQCIFQYFKDQQETNCIWTFSPFSEHICQRVCQVYWNWNTPSSRTSDPRELFANSTIHNNRCGEWRYMFYHTRTLVQTAPLMKETLHSDGKHSMYCEQRHFFRSSYLIKVNYDVSHYLELYTFSEIPWKLTTHGWDGFSWFLLVNSCCTFDIDGKCGILSQCGMSRAFRTRQEDAYHFQTSLMHLHLHLHVQEGKHEKADLFAQFYNMLPMHGGTCGRNTEPSDLFDSATMNKYMAEHPASCKACPNVSKECFVYWWSHDFTKKHKLIEFSCGRDTGQTTQRTWNVDENEGGKWIWRFHYFMRAKALQIDPKFKPYWNEPRAIMRPGHVTAAPCICAQHSQNETAVCNRDQMHIHAIEFQQYHSRAFGEVQTWCDIGKENENDFIYEQHWWLVGGTEGMAGVIWKFVCARCRTQDCDFWKTCLTYSAQPMMKVYDTIFYVNSINPWEFEDHPSQCDKCVQSIPTDAKYAICGKFVISHWLYWTPQKFEECVHNNVRCAPMGNRLWGTACMVIQNVWLRPSMGSHFSCILNVGGSNINIQGKETWTHVPILHMHEIDLISTASSGMETCKPCFLSGPTIHMGFSYEIRAQPYSRDYFCMDWMQEADEVDHNRCETVQPTLPLLQQFEWKTSCMGQRWITIFCDHCQIVCFSTFFCVMPTFLPNTSILDKFYCIYLSISWTHYCNVHALGFIMRLHYSYMGWKEHKRMHAWDIGLDELWAQEGIQRAQLWCGDEFEVAKYPEWITEARTAIATRPWFHNCYIKPWWIREKHLWFGKESKLDHGHRGAMFTPVANDNTEWMHHWYMFCWAGSKNRLKRQIKEKLIFIIKFMITEFGLFLMIDYTQCYIAWMWAYTGIACYIDWEKCLKHDLTTTDLGCCVYRLFKWYEVRHRAPPQVNTRLPWSQIPMVAIQCNIVDECKEQWHFSYKASFVVEYLCPGCCTNGNRWQWYQVKETPFMYAFAASIFGFHHENLVVFITGSVTIPNGLFGCIAWTSPKPVQKTPASANTIIAYDKCILMG'
    res = align.middle_edge(s1, s2, alphabet, scoring_matrix, sigma=5)
    assert res[1] == ((512, 510), (513, 511))


def test_linear_space_alignment():
    alphabet, scoring_matrix = align.get_blosum62()

    s1 = 'PS'
    s2 = 'IP'
    score, path = align.linear_space_alignment(s1, s2, alphabet, scoring_matrix, sigma=5)
    aligned = align.aligned_strings_from_path(s1, s2, path)
    assert aligned == ('-PS', 'IP-')

    s1 = 'PLE'
    s2 = 'MEA'
    score, path = align.linear_space_alignment(s1, s2, alphabet, scoring_matrix, sigma=5)
    aligned = align.aligned_strings_from_path(s1, s2, path)
    assert score == -3
    assert aligned == ('PLE-', '-MEA')

    s1 = 'PLEASANTLY'
    s2 = 'MEANLY'
    score, path = align.linear_space_alignment(s1, s2, alphabet, scoring_matrix, sigma=5)
    aligned = align.aligned_strings_from_path(s1, s2, path)
    assert score == 8
    assert aligned == ('PLEASANTLY', '-MEA--N-LY')

    s1 = 'PTGQQVPFPTVDIVCCTGIKCPFNYHMASIMDSYVFLQVPFPTTVDDVICCTGIKEPMNVGYDQQMTCFCKNYHMSVKDAAYDGDKEMDGMTKWCVMPNCMWENEAQDQMQAWDS'
    s2 = 'QVPFPTVVDVIVCCTGIKCEPMNVGYDQQMKDCFICTREYDIRRLHTIVCGSEWACRLWIEADWEDCEKSFRDFDAPINIVQYAVWRANVE'
    score, path = align.linear_space_alignment(s1, s2, alphabet, scoring_matrix, sigma=5)
    aligned = align.aligned_strings_from_path(s1, s2, path)
    assert score == 13
    assert aligned == (
        'PTGQQVPFPT-VD-IVCCTGIKC-PFNY-HMASIMDSYVFLQVPFP-TTVDDVICCTGIKEPMNVGYDQQMTCFCKNYHMSVKDAAYDGDKEMDGMTKWCVMPNCMWENEAQDQMQAWDS',
        '---Q-VPFPTVVDVIVCCTGIKCEPMNVGYDQQMKDCFICTR-EYDIRRLHTIVC--GSEWACRLWIEADWED-CEK---SFRD--FDAP--IN-IVQYAV-----WR--A-N-VE----')


def test_three_way_alignment():
    s1 = 'ACGATACGT'
    s2 = 'CCCATTAAGT'
    s3 = 'GACTATAGAA'
    score, alignment = align.three_way_alignment(s1, s2, s3)

    s1 = 'ATATCCG'
    s2 = 'TCCGA'
    s3 = 'ATGTACTG'
    score, alignment = align.three_way_alignment(s1, s2, s3)
    assert score == 3
    assert alignment == ('AT-ATCCG-', '-T---CCGA', 'ATGTACTG-')

    s1 = 'TGTTTAAAAATGTCCGCAACCATTTC'
    s2 = 'GATATAAAACAGGGATAACTGCAATGG'
    s3 = 'CCTGCTACTTTATGCCGTCTCCATATGCG'
    score, alignment = align.three_way_alignment(s1, s2, s3)
    assert score == 11
    assert alignment == (
        '--TGTTTA--AAAATG---T--CCGCAACCATTTC', '---G-ATATAAAACAGGGATAACTGC-A--AT-GG',
        'CCTG-CTACTTTA-TGCCGT--C-TCCA-TATGCG')

    s1 = 'A'
    s2 = 'AT'
    s3 = 'A'
    score, alignment = align.three_way_alignment(s1, s2, s3)
    assert score == 1
    assert alignment == ('A-', 'AT', 'A-')

    s1 = 'AAAAT'
    s2 = 'CCCCT'
    s3 = 'T'
    score, alignment = align.three_way_alignment(s1, s2, s3)
    assert score == 1
    assert alignment == ('AAAAT', 'CCCCT', '----T')

    s1 = 'AT'
    s2 = 'ACCT'
    s3 = 'AGGGGT'
    score, alignment = align.three_way_alignment(s1, s2, s3)
    assert score == 2
    assert alignment == ('A----T', 'A--CCT', 'AGGGGT')

    s1 = 'GGAG'
    s2 = 'TT'
    s3 = 'CCCC'
    score, alignment = align.three_way_alignment(s1, s2, s3)
    assert score == 0
    assert alignment == ('GGAG', '--TT', 'CCCC')

    s1 = 'T'
    s2 = 'T'
    s3 = 'T'
    score, alignment = align.three_way_alignment(s1, s2, s3)
    assert score == 1
    assert alignment == ('T', 'T', 'T')

    s1 = 'CGCGGGCATG'
    s2 = 'TATACGCGCG'
    s3 = 'CAGTACGT'
    score, alignment = align.three_way_alignment(s1, s2, s3)
    assert score == 4
    assert alignment == ('----C-G--CGGGCATG', 'TATAC-G--C--G--CG', '----CAGTAC--G---T')
