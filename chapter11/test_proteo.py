from pathlib import Path
import networkx as nx
from numpy import isclose
import proteo
from stepik_proteo import fetch_ints, pretty_print_adj, fetch_ints_and_string, fetch_psm_search_input, \
    fetch_spectral_dict_input, fetch_spectral_alignment_input


def test_graph_from_spectrum():
    spectrum = [57, 71, 154, 185, 301, 332, 415, 429, 486]
    graph = proteo.graph_from_spectrum(spectrum)
    # nx.readwrite.write_gpickle(graph, Path('test/testcase01.pickle'))
    expected = nx.readwrite.read_gpickle(Path('test/testcase01.pickle'))
    assert nx.is_isomorphic(graph, expected)

    spectrum = [113, 137, 184, 268, 297, 381, 425, 496, 581, 624, 695, 752, 832, 889, 988, 1045, 1102, 1160, 1215, 1217,
                1314, 1373, 1401, 1460, 1557, 1559, 1614, 1672, 1729, 1786, 1885, 1942, 2022, 2079, 2150, 2193, 2278,
                2349, 2393, 2477, 2506, 2590, 2637, 2661, 2774]
    graph = proteo.graph_from_spectrum(spectrum)
    # nx.readwrite.write_gpickle(graph, Path('test/testcase02.pickle'))
    expected = nx.readwrite.read_gpickle(Path('test/testcase02.pickle'))
    assert nx.is_isomorphic(graph, expected)

    spectrum = [115, 147, 216, 260, 287, 358, 416, 521, 563, 620, 652, 781, 806, 937, 953, 1008, 1081, 1095, 1152, 1194,
                1281, 1323, 1394, 1444, 1525, 1575, 1646, 1688, 1775, 1817, 1874, 1888, 1961, 2016, 2032, 2163, 2188,
                2317, 2349, 2406, 2448, 2553, 2611, 2682, 2709, 2753, 2822, 2854, 2969]
    graph = proteo.graph_from_spectrum(spectrum)
    # nx.readwrite.write_gpickle(graph, Path('test/testcase03.pickle'))
    expected = nx.readwrite.read_gpickle(Path('test/testcase03.pickle'))
    assert nx.is_isomorphic(graph, expected)


def test_ideal_spectrum():
    spectrum = proteo.ideal_spectrum('GPG')
    assert spectrum == [0, 57, 57, 154, 154, 211]

    spectrum = proteo.ideal_spectrum('GPFNA')
    assert spectrum == [0, 57, 71, 154, 185, 301, 332, 415, 429, 486]

    spectrum = proteo.ideal_spectrum('CRQCSLAMQRASQHYVYVWPQETFGFVCRM')
    assert spectrum == [0, 103, 131, 259, 287, 387, 390, 489, 490, 577, 636, 690, 693, 761, 840, 892, 941, 1020, 1070,
                        1176, 1198, 1247, 1295, 1334, 1462, 1481, 1580, 1599, 1743, 1762, 1842, 1861, 2005, 2024, 2123,
                        2142, 2270, 2309, 2357, 2406, 2428, 2534, 2584, 2663, 2712, 2764, 2843, 2911, 2914, 2968, 3027,
                        3114, 3115, 3214, 3217, 3317, 3345, 3473, 3501, 3604]


def test_peptide_from_ideal_spectrum():
    spectrum = [57, 71, 154, 185, 301, 332, 415, 429, 486]
    peptide = proteo.peptide_from_ideal_spectrum(spectrum)
    assert peptide == 'GPFNA'

    spectrum = [56, 71, 154, 185, 301, 332, 415, 429, 486]
    peptide = proteo.peptide_from_ideal_spectrum(spectrum)
    assert peptide == None

    spectrum = [103, 131, 259, 287, 387, 390, 489, 490, 577, 636, 690, 693, 761, 840, 892, 941, 1020, 1070,
                1176, 1198, 1247, 1295, 1334, 1462, 1481, 1580, 1599, 1743, 1762, 1842, 1861, 2005, 2024, 2123,
                2142, 2270, 2309, 2357, 2406, 2428, 2534, 2584, 2663, 2712, 2764, 2843, 2911, 2914, 2968, 3027,
                3114, 3115, 3214, 3217, 3317, 3345, 3473, 3501, 3604]
    peptide = proteo.peptide_from_ideal_spectrum(spectrum)
    assert peptide == 'CRQCSLAMQRASQHYVYVWPQETFGFVCRM'

    spectrum = [57, 185, 186, 300, 348, 429, 504, 585, 651, 722, 814, 835, 917, 938, 1009, 1014, 1096, 1101, 1209,
                1216, 1287, 1372, 1402, 1459, 1473, 1530, 1560, 1658, 1674, 1771, 1787, 1885, 1915, 1972, 1986, 2043,
                2073, 2158, 2229, 2236, 2344, 2349, 2431, 2436, 2507, 2528, 2610, 2631, 2723, 2794, 2860, 2941, 3016,
                3097, 3145, 3259, 3260, 3388, 3445]
    peptide = proteo.peptide_from_ideal_spectrum(spectrum)
    assert peptide == 'GQYRFYCPSDADGAQLNSTYLSACLHRENW'


def test_peptide_vector_from_peptide():
    pept_vect = proteo.vector_from_peptide('TFPRGPHSPRVVDLRCCKQMNDHKSIDWKYSLYFM')
    expected = fetch_ints(Path('test/testcase04.txt'))
    assert pept_vect == expected

    pept_vect = proteo.vector_from_peptide('YRWAIPQIVYRFWDELWLVWNGQVDFMDP')
    expected = fetch_ints(Path('test/testcase05.txt'))
    assert pept_vect == expected


def test_peptide_from_vector():
    vector = fetch_ints(Path('test/testcase06.txt'))
    peptide = proteo.peptide_from_vector(vector)
    assert peptide == 'FRVLRCLEDQPLGLMAQPHRTMPDHPQFTGDAHHMYATC'

    vector = fetch_ints(Path('test/testcase07.txt'))
    peptide = proteo.peptide_from_vector(vector)
    assert peptide == 'CNWSTDLSQEFTCAQFTYYYMLMQQLTW'

    vector = [0] * 55 + [1]
    peptide = proteo.peptide_from_vector(vector)
    assert peptide is None


def test_peptide_from_spectral_vector():
    spectrum = fetch_ints(Path('test/testcase08.txt'))
    peptide = proteo.peptide_from_spectral_vector(spectrum)
    assert peptide == 'GGPGGPGGAGG'

    spectrum = fetch_ints(Path('test/testcase09.txt'))
    peptide = proteo.peptide_from_spectral_vector(spectrum)
    assert peptide == 'DAGGGGGGTGGV'

    spectrum = fetch_ints(Path('test/testcase10.txt'))
    peptide = proteo.peptide_from_spectral_vector(spectrum)
    assert peptide == 'CGEAGLQG'


def test_identify_peptide_from_proteome():
    spectrum, proteome = fetch_ints_and_string(Path('test/testcase11.txt'))
    peptide, score = proteo.identify_peptide_from_proteome(spectrum, proteome)
    assert peptide == 'KLEAARSCFSTRNE'
    assert score == 274

    spectrum, proteome = fetch_ints_and_string(Path('test/testcase12.txt'))
    peptide, score = proteo.identify_peptide_from_proteome(spectrum, proteome)
    assert peptide == 'SQSVIKFTESATGGN'
    assert score == 234

    spectrum, proteome = fetch_ints_and_string(Path('test/testcase13.txt'))
    peptide, score = proteo.identify_peptide_from_proteome(spectrum, proteome)
    assert peptide == 'LQKTIIAFHSHVHT'
    assert score == 142


def test_psm_search():
    spectra, proteome, threshold = fetch_psm_search_input(Path('test/testcase14.txt'))
    peptides = proteo.psm_search(proteome, spectra, threshold)
    assert sorted(peptides) == sorted(['QQCGVHEYFWVSKK',
                                       'HTNGPDCSQYQLLK',
                                       'VIAAGAHPADGQGVRGP',
                                       'NGMPFCCMCWDVVM',
                                       'AAPVCLQQMQPKAVL',
                                       'SIAQIMVEYTVHGH',
                                       'KMARKRHIHKFLSP',
                                       'NRAEQFDMTKYCV',
                                       'ADMCRPCQACTGKAFG',
                                       'CKFADFDSKTMGVITQ',
                                       'DETTVPHLVCPWHD',
                                       'IFWVHEMMYHCE',
                                       'GWKRGTYEIIFCPP',
                                       'DGQGVRGPHQIILMVR',
                                       'TCFAAGAHVMRKGCH',
                                       'DCQNYMLMHMVETG',
                                       'CYCMFHTNTARGERK'])

    spectra, proteome, threshold = fetch_psm_search_input(Path('test/testcase15.txt'))
    peptides = proteo.psm_search(proteome, spectra, threshold)
    assert sorted(peptides) == sorted(['MIVALRDMFFFPR',
                                       'AGFQCLEGVDHAMKK',
                                       'PYTCKGAHADCPPCAG',
                                       'CQLTTTKHCWSWDP',
                                       'DPTCVNMSLAISFVYCG',
                                       'HSDVQTENSNNPAVPM',
                                       'MWKYGDFVTCIDP',
                                       'HHRANTPLVAEGAIIC',
                                       'LKQKDWCGISRADD',
                                       'LPSDQIKILERVM',
                                       'MAVWTCWHCGHNAT',
                                       'LNFEVVVTMWALWL',
                                       'IGRIPVEHQQFMAC',
                                       'KHEGCYRPECTVW',
                                       'AWLYPPYRFESFC',
                                       'NCGQFARGWCGCTTA',
                                       'TSQPPIYIRSHVNKT',
                                       'AKYEVHIVSTDRGYN',
                                       'MAMEWLSFHMQR',
                                       'VWVQGNAAIAKRHRF'])


def test_size_of_spectral_dict():
    spectrum, threshold, max_score = fetch_spectral_dict_input(Path('test/testcase16.txt'))
    size = proteo.size_of_spectral_dict(spectrum, threshold, max_score)
    assert size == 330

    spectrum, threshold, max_score = fetch_spectral_dict_input(Path('test/testcase17.txt'))
    size = proteo.size_of_spectral_dict(spectrum, threshold, max_score)
    assert size == 336

    spectrum, threshold, max_score = fetch_spectral_dict_input(Path('test/testcase18.txt'))
    size = proteo.size_of_spectral_dict(spectrum, threshold, max_score)
    assert size == 293


def test_prob_of_spectral_dict():
    spectrum, threshold, max_score = fetch_spectral_dict_input(Path('test/testcase19.txt'))
    prob = proteo.prob_of_spectral_dict(spectrum, threshold, max_score)
    assert isclose(prob, 0.00132187890625)

    spectrum, threshold, max_score = fetch_spectral_dict_input(Path('test/testcase20.txt'))
    prob = proteo.prob_of_spectral_dict(spectrum, threshold, max_score)
    assert isclose(prob, .00012864921875)

    spectrum, threshold, max_score = fetch_spectral_dict_input(Path('test/testcase21.txt'))
    prob = proteo.prob_of_spectral_dict(spectrum, threshold, max_score)
    assert isclose(prob, 8.381796875e-05)


def test_spectral_alignment():
    peptide, spectrum, k = fetch_spectral_alignment_input(Path('test/testcase23.txt'))
    alignment, score = proteo.spectral_alignment(peptide, spectrum, k)
    assert alignment == 'L(-61)VW(-9)STE(+69)'
    assert score == 67

    peptide, spectrum, k = fetch_spectral_alignment_input(Path('test/testcase24.txt'))
    alignment, score = proteo.spectral_alignment(peptide, spectrum, k)
    assert (alignment, score) == ('Y(-69)CRC(+69)N', 45)

    peptide, spectrum, k = fetch_spectral_alignment_input(Path('test/testcase25.txt'))
    alignment, score = proteo.spectral_alignment(peptide, spectrum, k)
    assert (alignment, score) == ('A(-20)EE(-9)IN(+29)', 62)
