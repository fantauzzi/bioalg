from pathlib import Path
import networkx as nx
import proteo
from stepik_proteo import fetch_ints, pretty_print_adj, fetch_ints_and_string


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
    peptide = proteo.identify_peptide_from_proteome(spectrum, proteome)
    assert peptide == 'KLEAARSCFSTRNE'

    spectrum, proteome = fetch_ints_and_string(Path('test/testcase12.txt'))
    peptide = proteo.identify_peptide_from_proteome(spectrum, proteome)
    assert peptide == 'SQSVIKFTESATGGN'

    spectrum, proteome = fetch_ints_and_string(Path('test/testcase13.txt'))
    peptide = proteo.identify_peptide_from_proteome(spectrum, proteome)
    assert peptide == 'LQKTIIAFHSHVHT'
