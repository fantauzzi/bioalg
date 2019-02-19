import sequence
from pathlib import Path


def fetch_line_from_file(file_name):
    with open(file_name) as input_file:
        line = input_file.readline().rstrip().split(' ')
    int_list = [int(item) for item in line]
    return int_list


def test_translate_RNA():
    rna = 'AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA'
    expected = 'MAMAPRTEINSTRING'
    res = sequence.translate_RNA(rna)
    assert res == expected

    assert sequence.translate_RNA('GAA') == 'E'
    assert sequence.translate_RNA('UAG') == ''


def test_peptide_encoding():
    dna = 'ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA'
    ammino = 'MA'
    res = sequence.peptide_encoding(dna, ammino)
    assert sorted(res) == sorted(['ATGGCC', 'GGCCAT', 'ATGGCC'])


def test_cyclic_spectrum():
    peptide = 'LEQN'
    res = sequence.peptide_spectrum(peptide)
    assert res == [0, 113, 114, 128, 129, 227, 242, 242, 257, 355, 356, 370, 371, 484]

    peptide2 = 'YYTMDKDVWRAP'
    res2 = sequence.peptide_spectrum(peptide2)
    assert res2 == [0, 71, 97, 99, 101, 115, 115, 128, 131, 156, 163, 163, 168, 186, 214, 227, 232, 243, 243, 246, 260,
                    264, 285, 324, 326, 331, 342, 342, 347, 358, 374, 395, 400, 413, 423, 427, 441, 457, 475, 487, 489,
                    494, 510, 510, 512, 524, 528, 556, 558, 588, 590, 595, 609, 627, 638, 643, 650, 655, 673, 673, 684,
                    689, 724, 726, 751, 753, 755, 770, 772, 774, 799, 801, 836, 841, 852, 852, 870, 875, 882, 887, 898,
                    916, 930, 935, 937, 967, 969, 997, 1001, 1013, 1015, 1015, 1031, 1036, 1038, 1050, 1068, 1084, 1098,
                    1102, 1112, 1125, 1130, 1151, 1167, 1178, 1183, 1183, 1194, 1199, 1201, 1240, 1261, 1265, 1279,
                    1282, 1282, 1293, 1298, 1311, 1339, 1357, 1362, 1362, 1369, 1394, 1397, 1410, 1410, 1424, 1426,
                    1428, 1454, 1525]


def test_linear_spectrum():
    peptide = 'NQEL'
    res = sequence.peptide_spectrum(peptide, cyclic=False)
    assert res == [0, 113, 114, 128, 129, 242, 242, 257, 370, 371, 484]

    res2 = sequence.peptide_spectrum('VNLMFIITTVLEWGFTKAFHDLNCMNKVGRFPSVYGGCRKTSHESG', cyclic=False)
    expected = fetch_line_from_file(Path('test/testcase01.txt'))

    assert res2 == expected


def print_cyclopeptide_sequences(seqs):
    for item in seqs:
        print(*item, sep='-', end='')
        print(' ', sep='', end='')
    print()


def test_sequence_cyclopeptide():
    spectrum = [0, 97, 97, 99, 101, 103, 196, 198, 198, 200, 202, 295, 297, 299, 299, 301, 394, 396, 398, 400, 400, 497]
    res = sequence.sequence_cyclopeptide(spectrum)
    assert sorted(res) == sorted(
        [[97, 99, 103, 97, 101], [97, 101, 97, 99, 103], [97, 101, 97, 103, 99], [97, 103, 99, 97, 101],
         [99, 97, 101, 97, 103], [99, 103, 97, 101, 97], [101, 97, 99, 103, 97], [101, 97, 103, 99, 97],
         [103, 97, 101, 97, 99], [103, 99, 97, 101, 97]])

    spectrum2 = [0, 113, 128, 186, 241, 299, 314, 427]
    res2 = sequence.sequence_cyclopeptide(spectrum2)
    assert sorted(res2) == sorted(
        [[113, 128, 186], [113, 186, 128], [128, 113, 186], [128, 186, 113], [186, 113, 128], [186, 128, 113]])

    spectrum3 = [0, 87, 99, 113, 115, 115, 137, 137, 147, 147, 202, 228, 234, 236, 246, 252, 252, 260, 284, 339, 347,
                 349, 351, 365, 375, 383, 383, 399, 438, 462, 462, 486, 498, 498, 512, 512, 520, 577, 585, 585, 599,
                 599, 611, 635, 635, 659, 698, 714, 714, 722, 732, 746, 748, 750, 758, 813, 837, 845, 845, 851, 861,
                 863, 869, 895, 950, 950, 960, 960, 982, 982, 984, 998, 1010, 1097]
    res3 = sequence.sequence_cyclopeptide(spectrum3)
    assert sorted(res3) == sorted(
        [[87, 115, 137, 99, 147, 137, 115, 113, 147], [87, 147, 113, 115, 137, 147, 99, 137, 115],
         [99, 137, 115, 87, 147, 113, 115, 137, 147], [99, 147, 137, 115, 113, 147, 87, 115, 137],
         [113, 115, 137, 147, 99, 137, 115, 87, 147], [113, 147, 87, 115, 137, 99, 147, 137, 115],
         [115, 87, 147, 113, 115, 137, 147, 99, 137], [115, 113, 147, 87, 115, 137, 99, 147, 137],
         [115, 137, 99, 147, 137, 115, 113, 147, 87], [115, 137, 147, 99, 137, 115, 87, 147, 113],
         [137, 99, 147, 137, 115, 113, 147, 87, 115], [137, 115, 87, 147, 113, 115, 137, 147, 99],
         [137, 115, 113, 147, 87, 115, 137, 99, 147], [137, 147, 99, 137, 115, 87, 147, 113, 115],
         [147, 87, 115, 137, 99, 147, 137, 115, 113], [147, 99, 137, 115, 87, 147, 113, 115, 137],
         [147, 113, 115, 137, 147, 99, 137, 115, 87], [147, 137, 115, 113, 147, 87, 115, 137, 99]])


def test_score_peptide():
    pept = 'NQEL'
    spectrum = [0, 99, 113, 114, 128, 227, 257, 299, 355, 356, 370, 371, 484]
    score = sequence.score_peptide(pept, spectrum)
    assert score == 11

    pept2 = 'EERSDSDAQFLWAGYRCYPWPQKRTPLRGDAAQERGSWVAMKCSHYRS'
    spectrum2 = fetch_line_from_file(Path('test/testcase02.txt'))
    score2 = sequence.score_peptide(pept2, spectrum2)
    assert score2 == 679

    score3 = sequence.score_peptide('NQEL', spectrum, cyclic=False)
    assert score3 == 8

    spectrum3 = fetch_line_from_file(Path('test/testcase03.txt'))
    score4 = sequence.score_peptide('ADEQMQSPHEIDPMRYVSAQLADRTPWPVRLLRGHSD', spectrum3, cyclic=False)
    assert score4 == 201


def test_leaderboard_peptide_sequence():
    spectrum = [0, 71, 113, 129, 147, 200, 218, 260, 313, 331, 347, 389, 460]
    res = sequence.leaderboard_peptide_sequence(spectrum, 10)
    assert res == [113, 147, 71, 129]

