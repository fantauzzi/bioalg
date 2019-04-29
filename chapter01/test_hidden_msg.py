from random import choices, seed, randint
from pathlib import Path
import hidden_msg as hm
from stepik_hidden_msg import fetch_clumps_finding_input


def test_skew():
    genome = 'CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAACTAGTTCACACGGCTTGATGGCAAATGGTTTTTCCGGCGACCGTAATCGTCCACCGAG'
    running_skew = [0] + hm.compute_skew(genome)
    mins = hm.argsmin(running_skew)
    assert mins == [53, 97]

    with open('skew_dataset.txt') as genome_file:
        genome = genome_file.read()
    running_skew = [0] + hm.compute_skew(genome)
    mins = hm.argsmin(running_skew)
    assert sorted(mins) == sorted([89969, 89970, 89971, 90345, 90346])

    genome = 'TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT'
    running_skew = [0] + hm.compute_skew(genome)
    mins = hm.argsmin(running_skew)
    assert mins == [11, 24]

    genome = 'CATTCCAGTACTTCGATGATGGCGTGAAGA'
    running_skew = [0] + hm.compute_skew(genome)
    mins = hm.argsmin(running_skew)
    assert mins == [14]


def test_DNA_Complement():
    res = hm.DNA_rev_complement('GATTACA')
    assert res == 'TGTAATC'

    res = hm.DNA_rev_complement('TTGTGTC')
    assert res == 'GACACAA'

    seed(42)
    for i in range(0, 10):
        dna = ''.join(choices(['A', 'T', 'C', 'G'], k=randint(50, 100)))
        compl = hm.DNA_rev_complement(dna)
        assert dna == hm.DNA_rev_complement(compl)

    res = hm.DNA_rev_complement('AAAACCCGGT')
    assert res == 'ACCGGGTTTT'


def test_pattern_matching():
    text = 'GCGCG'
    pattern = 'GCG'
    count = hm.approximate_pattern_count(text, pattern, 0)
    assert count == 2

    text = 'CGCGATACGTTACATACATGATAGACCGCGCGCGATCATATCGCGATTATC'
    count = hm.approximate_pattern_count(text, 'CGCG', 0)
    count == 3

    text = 'ACCACAAAACAATAAACAATTGTAAAAACAATAAACAATAAACAATAAACAATAAACAATAAACAATACCCGTGAAACAATTAAACAATATTGAAACAATGATTGATTAAAACAATGCCTCTTATCTAGTCGAAACAATAAACAATAAACAATCAAACAATACCCGAAACAATCTGACCCAAAACAATAAACAATTGCAAACAATACCAAACAATAAACAATACCAAACAATGAAACAATAAACAATAGTTCAAACAATGTCATACTGGAAACAATATTGTAAACAATCCTTTTGAACAAAACAATAGCCGATGGGGAGAAAAACAATCAAACAATAAACAATAAACAATTACTAAACAATTGTCGCAAACAATTAAACAATAAAGAAAAACAATAAGGTGAAACAATAAACAATCCCAATACAAACAATGAGAAACAATAAAACAATGCTTAAACAATATAAACAATGGGAAAGAAACAATAAACAATCTGTACCAAACAATTAAACAATCTCATCAAAACAATTAAACAATAAACAATAAAACAATGAAACAATAAACAATACAAACAATACAAACAATTAAACAATAAACAATGAAATAAACAATCAAACAATGTCGCTAGACCCAAACAATGGGGAAACAATAGAAACAATGTTCTAAACAATAAAAAACAATAAACAATCTGTCTATAAACAATAAACAATAAACAATAAACAATAAACAATCGCAAACAATGCAAAACAATTTCGAAACAATCAAACAATAAACAATGTAAACAATAAACAATTTGGAAACAATAAAACAATTAAACAATGAGAAACAATAAAACAATCAAACAATAAAACAATAAACAATCAAGGTAAACAATGTAAACAATAGGACAGTAAACAATTATAAACAATAAACAATGAAACAATCGGTTACACACTAAACAATTAAAAACAATAGAAACAATATCAAACAATAGAAACAATTCCCGAAGGAAACAATAAACAATATAAAACAATAAACAATAAACAATTCGTAAACAAT'
    pattern = 'AAACAATAA'
    count = hm.approximate_pattern_count(text, pattern, 0)
    assert count == 38


def test_frequent_words_with_mismatches():
    text = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
    kmers = hm.frequent_words_with_mismatches(text, k=4, d=0)
    assert sorted(kmers) == ['CATG', 'GCAT']

    text = 'GAACTTCCAAGCGCGAACTTCGAACTTCCTACATGACTGCGTCAAGCTACATGACCTACATGACCTACATGACTGCGTCAAGTGCGTCAAGCTACATGACTGCGTCAAGCTACATGACCAAGCGCGAACTTCGAACTTCCTACATGACAACATGCCAAGCGCGAACTTCGAACTTCTGCGTCAAGGAACTTCCAAGCGCCTACATGACGAACTTCGAACTTCCTACATGACGAACTTCCTACATGACCAAGCGCAACATGCAACATGCCTACATGACCAAGCGCAACATGCCTACATGACGAACTTCTGCGTCAAGCAAGCGCGAACTTCCAAGCGCAACATGCTGCGTCAAGCAAGCGCCTACATGACCAAGCGCGAACTTCTGCGTCAAGCAAGCGCCTACATGACTGCGTCAAGCTACATGACCTACATGACGAACTTCCAAGCGCGAACTTCCTACATGACTGCGTCAAGCAAGCGCGAACTTCCTACATGACTGCGTCAAGCAAGCGCTGCGTCAAGAACATGCCAAGCGCCTACATGACCAAGCGCCAAGCGCAACATGCTGCGTCAAGCTACATGACTGCGTCAAGCTACATGACAACATGCCAAGCGCCAAGCGCGAACTTCAACATGCGAACTTCAACATGCTGCGTCAAGGAACTTCTGCGTCAAGCTACATGACCAAGCGCCTACATGACAACATGCCTACATGACCTACATGACGAACTTCTGCGTCAAGAACATGCCTACATGACTGCGTCAAGTGCGTCAAGCAAGCGCAACATGCGAACTTCTGCGTCAAGTGCGTCAAGAACATGC'
    kmers = hm.frequent_words_with_mismatches(text, k=11, d=0)
    assert kmers == ['CTGCGTCAAGC']

    text = 'CGGAAGCGAGATTCGCGTGGCGTGATTCCGGCGGGCGTGGAGAAGCGAGATTCATTCAAGCCGGGAGGCGTGGCGTGGCGTGGCGTGCGGATTCAAGCCGGCGGGCGTGATTCGAGCGGCGGATTCGAGATTCCGGGCGTGCGGGCGTGAAGCGCGTGGAGGAGGCGTGGCGTGCGGGAGGAGAAGCGAGAAGCCGGATTCAAGCAAGCATTCCGGCGGGAGATTCGCGTGGAGGCGTGGAGGCGTGGAGGCGTGCGGCGGGAGATTCAAGCCGGATTCGCGTGGAGAAGCGAGAAGCGCGTGCGGAAGCGAGGAGGAGAAGCATTCGCGTGATTCCGGGAGATTCAAGCATTCGCGTGCGGCGGGAGATTCAAGCGAGGAGGCGTGAAGCAAGCAAGCAAGCGCGTGGCGTGCGGCGGGAGAAGCAAGCGCGTGATTCGAGCGGGCGTGCGGAAGCGAGCGG'
    kmers = hm.frequent_words_with_mismatches(text, k=12, d=0)
    assert sorted(kmers) == ['CGGCGGGAGATT', 'CGGGAGATTCAA', 'CGTGCGGCGGGA', 'CGTGGAGGCGTG', 'CGTGGCGTGCGG',
                             'GCGTGCGGCGGG', 'GCGTGGAGGCGT', 'GCGTGGCGTGCG', 'GGAGAAGCGAGA', 'GGAGATTCAAGC',
                             'GGCGGGAGATTC', 'GGGAGATTCAAG', 'GTGCGGCGGGAG', 'TGCGGCGGGAGA']

    text = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
    k = 4
    d = 1
    kmers = hm.frequent_words_with_mismatches(text, k, d)
    assert sorted(kmers) == sorted(['GATG', 'ATGC', 'ATGT'])

    text = 'TCACGACGACTGATTCACGACGAATGATTGATTCACTGATGACGCCGCCTCACTCACGACGACGCCGACTCACGAAGCCGACGAAGACGAAGAATGATTGATTCACGACGCCGCCGACGAATCACGAAGCCTCACGCCGCCG'
    k = 6
    d = 3
    kmers = hm.frequent_words_with_mismatches(text, k, d)
    assert kmers == ['CCCCCC']


def test_approximate_pattern_matching():
    res = hm.approximate_pattern_matching(text='GATATATGCATATACTT', pattern='ATAT', d=0)
    assert res == [1, 3, 9]

    res = hm.approximate_pattern_matching(text='ATGACTTCGCTGTTACGCGC', pattern='CGC', d=0)
    assert res == [7, 15, 17]

    text = 'GCCTCCGATTTCCGATTTCCGATTTCCGATTATCCGATTCTCCGATTAATACTCCGATTGTTCCGATTTCCGATTTCCGATTTCCGATTTTGACATCCGATTTCCGATTCTGGTGAGCATCCGATTTCCGATTTCCGATTTCCGATTGGTTCCGATTTACTTCCGATTTCCGATTTGGGACAGTCCGATTTCCGATTTGGTCCGATTTCCGATTATCCGATTTTCCGATTGTTTTCGGTCCGATTAGTCCGATTTCCGATTTTCCGATTTCCGATTATCCGATTCGTTCCGATTTTCCGATTTCCGATTGTCTCCGATTGTGCCCAGTTTCCGATTCCATTCCGATTATCAAAGCTCCGATTTCCGATTTCTCCGATTGGGGAGGTGCTCCGATTCTCCGATTCAGTCCGATTATGTTTTCCGATTACGGTCTCCGATTTCCGATTTCCGATTTCCGATTGATGCCTATCCGATTAGAGTCCGATTGTCATAGCTTTATCCGATTTCCGATTGGCCGGAACTCCGATTCGATCCGATTTTGGATCCGATTGTACGTCCGATTATCCGATTTCCGATTATTCCGATTTTCCGATTTCCGATTTGATCCGATTCCCGCGCTCATCCGATTTCCGATTTCCGATTGTCCGATTTCCGATTTTAGATTCCGATTTCCGATTAGCCTCTCCGATTAATTGGTCAATTCCGATTTCCGATTTCCGATTATGATCCGATTCTGCTCCGATTCTCCGATTTCCGATTAGGGGTACTCCATCCGATTTCCGATTCTCCGATTTTCCGATTGTTCCGATTGGTCCGATTTCCGATTTTCTCCGATTTCCGATTTGCAATTCCGATTGGTCCGATTATCCGATTCTGTCCGATTTCCGATTTCCGATTATCCGATTCTCCGATTGACATACCTGTCCGATTAGTCCGATTGCCTCCGATTTATCCGATTGTCCGATTTCCGATTATCCGATTAGTACTCCGATTAGTTCCGATTATCTCAGTCCGATTCGTCCGATTGTCCGATTGGTCCGATTTTCCGATTTCCGATTATTTTCCGATTTCCGATTTCCGATTACCCTTTCCGATTATTCCGATTTTCCGATTTCCGATTATCCGATTAGTCCGATTAGACATCCGATTGCAGGCTTATTCCGATTTGAGGTCCGATTCTCCGATTATGTCCGATTCTCCGATTTCTCCGATTGTCCGATTTATCCGATTATTAGTCCGATTTCCGATTATCCGATTTGTTCCGATTTCCGATTCATCCGATTCTCCGATTTCCGATTAAATTTTCCGATTTCCGATTGTCCGATTGTTCCGATTCTTCCGATTCTCTCATAGCCAATCCGATTCGCTCCGATTACTCCGATTGCACAATATGTTCTCCGATTACGATCCGATTACTCCGATTTCCGATTCATCCGATTTCCGATTGTCCGATTTTCCGATTGTGCCAATCCGATTTCCGATTGCTCCGATTTTCCGATTTCCGATTTCCGATTTCTGACACATCCGATTGTCCGATTTGAGTCCGATTTCCGATTATCCGATTGTGTGCACTTTCCGATTCCAGGATCCGATTGTAATCCGATTGGCTGTCCGATTTCCGATTGGTACGGCAAGAGTAGTCCGAACTTACTCCGGTAAGTCCGATTTCCGATTTTCCGATTGCATCCGATTCTCCGATTAGTCCGATTAGTTCCGATTCATCCGATTCATCCGATTTCCGATTCTCCGATTTTCCGTCGCGGTGATGCTTCCGATTTTCCGATTTCCGATTCTTCATCCGATTTGCTCCGATTATCCTCGGATTCCGATTAGAATAGTTCCGATTGTTCCGATTTCCGATTCTTCCGATTCCGCCTCCGATTTCCGATTTCCGATTGTCCGATTAGGATCCGATTCTTTTTGTAATCCGATTTATTTCTATCCGATTGCATCCGATTCCTTCCGATTTTCCGATTTCCGATTAGTGAGTCCGATTAATCCGATTGGATCGGTCCGATTTGTTCCGATTAAATCTTCCGATTTAATCCGATTTCCGATTCACGTCCGATTGGAGCAGAGTCCGATTGTCGGTATAACTCCGATTTCGGAATTCCGATTCCGTCCGATTTCCGATTTATAACCAGGCAGCTCCGATTTTATGATCCGATTCTTTCCGATTTTTTCCGATTGATACTCTCCGATTTCCGATTTTCCGATTCCTCCGATTTCTCCGATTTCCGATTGACTCTTCCGATTGGGGGATTCCGATTTGTCCGATTTCCGATTAAGCTGTCCGATTTCCGATTTCCGATTTCCGATTCGTTCCGATTAATTTCCGATTGCTCCGATTATCCGATTTTCCGATTGCTAACTCCGATTCGTTATTTTCCGATTGCCGTCCGATTTTCCGATTAGACGCGCAAATCCGATTTTCCGATTCTCCGATTTAGTCCGATTTACATCCGATTTCCGATTTTGGTTCCGATTTCCGATTTATCCGATTTTCCGATTATTTCCGATTTCCGATTATTCCGATTAAACCATCCGATTAATCCGATTTCCGATTATCCGATTACCCCGCCTCCGATTGACTCCGATTTGAACATCCGATTCGGGCTCCGATTAAGACGCCGTCCGATTCTTTTCCGATTGCATTTTGTGTCCGATTATCCGATTCCGTCCGATTCAATAGCTCCGATTTCCGATTTTTCCGATTTCTCCGATTTTCCGATTTTAAGTCCGATTTCCGATTAGTCCGATTCGGGCGAGATTGCACTCCGATTCATTCCGATTTTCCGATTTTGTCAAGGACAATCCGATTAGCTCCGATTCGTCCGATTATCCGATTTCCGATTCGGTTCCGATTATCCGATTTCCGATTCTGCCTCCGATTAAGTGTCCGATTTCTTCCGATTTGGCTCATCCGATTTCCGATTAGCAAATTATATCCGATTTCCGATTCTCCGATTCTTCCGATTGCGAGTCCGATTACTCCGATTTCCGATTCAGCTTTCTCCGATTTCCGATTTCAATCCGATTTTCCGATTCGACTCAGCGTCCGATTCATCCGATTGCCATCCGATTAAGTCCGATTCCTCCGATTGATCCTTCCGATTGCTCCGATTTCCGATTCGTGGTCCGATTTCCGATTGTTTCCGATTTTCCGATTCATCCGATTGTCCGATTATCCGATTTCCGATTCCAGGCTGTCCGATTAGGATCCGATTATCCGATTGGTCCGATTCTCCGATTTATCCGATTAAAATCCGATTGTCCGATTTCCGATTCCGTAAATCCCGATTTTCCGATTTCTCCGATTTTGCTCCGATTTCCGATTGTCCGATTCTCTTATGTCCGATTCAATCCGATTTCCGATTATCCGATTTGATCCGATTGCCGCGTCCGATTGGCTCCGATTGGTCCGATTTCCGATTGTCCGATTTCCGATTTCCATCCGATTAGTCCGATTCTCCGATTCACTCCGATTACGCCGGAATTCCGATTTCAGGGCTCAATTCCGATTGCGATCCGATTCTTCCGATTTTCCGATTTCCGATTTCCGATTTCCGATTTAATCCGATTGGGGTTGTCCTCCGATTAAGTAATTCCGATTATCCTACAATCCGATTTCCGATTTCCGATTGTTCCGATTAGGTCCGATTTCCGATTTCCGATTCTAAGTGTCCGATTTCCGATTGTCCGATTTCCGATTTATGTTATCCGATTATATCCGATTTGTCCGATTTCTCCGATTCGGCTTTCCGATTATCCGATTGCTGATCCGATTCTCCGATTATCCGATTTACAGCAAAAGGTCCGATTGTCCGATTCATCCGATTTTCCGATTGGTTTCCGATTTCCGATTTCCGATTCTGGTTCCGATTTCTGGTTCCGATTAGACGATCCGATTTCCGATTTTTCTCCGATTCCTCCGATTCAATTGCCATCCGATTGGCGATATCCGATTTCCGATTTCCGATTGTCCGATTGTCCGATTTCCGATTGCTGTCCGATTGTCCGATTTAACCGATCCGATTAAGGCTTCCGATTAGATCCGATTTCCGATTTCCGATTTCCGATTATCCGATTTGTGGTCTTCCGATTAACTCCGATTTTCCGATTTCCGATTAACCACTCCGATTATCCGATTGTTCCGATTTGGCTCCGATTAATCCGATTTCCGATTCAGACGTCCGATTAGAAGCTGGCGATCCGATTGTGGGGTCCGATTTCCGATTTATCGCCTCCGATTGCTCCGATTATTCCGATTGGTTGCTCCGATTTTCCGATTTGTCCGATTGTCATCCGATTTCCGATTACTATTTTCCGATTCATTTTCCGATTTCCGATTGCCCTGTCCGATTGGTCCGATTCAGTTCCGATTTCCGATTTGCATGAGTTGATCCGATTCCTTACGGAATCCCTCCGATTTCTTCCGATTCTCCGATTTTCCGATTTCCGATTACGAAGGTCCGATTAGTCCGATTCAGTATTCGAGTTGTCCGATTTGTCCGATTCCGCCGGTTCAGTCCGATTGTCCGATTACGCAACTATCCGATTCTCCGATTTTCCCCTTATCCGATTACTCCGATTCCCAGGTTCCGATTATTCCGATTTCCGATTCTGGATTCCGATTCTGTCCGATTTCCGATTTTCCGATTCCTTCCGATTGCCTCCGATTTCCGATTGTTCCGATTCCTCCGATTGCGTCCGATTTAGTTAACTATCCGATTCATCCGATTGTTCCGATTCCCTAGACATTCCGATTTCCGATTAAGGCGTCCGATTTCCGATTTTCCGATTTGTCTTCCGTTCCGATTCGTCCGATTGAGTTCGACTATCCGATTGTCCGATTGTGCATCCGATTCTAGCTAGTCCGATTTTCCGATTTTCGTCCGATTTTCCCACACATCCGATTAGAGCCATCCGATTTACATGTCCGATTTCCGATTTCCGATTTCAAATCCGATTTCCGATTGTGTTCCGATTGACTCCGATTTCCGATTTCCGATTCAATCTCCGATTTCCGATTGGCCTCCGATTCTTCCGATTCACTGCCGTCCGATTTTCCGATTTAAACTCTCCGATTCGGTTGAAGACTCCGATTTTCCGATTTGCGTATGGTGTCCGATTTCGTCCGATTTCCGATTATAGTCCGATTTTTTCCGATTTCCGATTAGGAGAATATCCGATTATCCGATTGATCCGATTAGTAAGCGGCATCCGATTTTCCGATTTCCGATTCCCGGCGTCCGATTGGGTCCGATTACGGATCCGATTTTCCGATTATCTCCGATTAATCCGATTAGTCCGATTCGATCACTCCGATTTCCGATTTCCGATTGTTCCGATTCACTCCGATTCCCTCCGATTTTTCCGATTTTCCGATTGGTCCGATTATAGATCCGATTTCCGATTTTCCGATTAGACATCCGATTTCCGATTGACCCCCCTCCGATTTTCCGATTTTCTCCGATTGTCCATACGACTCCGATTAAGGTATTCCGATTATACGTCCGATTCTCCGATTAATCCGATTTCCGATTAATCCGATTGATTCCGATTTCCGATTGTCCGATTTCCGATTCATATCCGATTGTCCGATTTCCGATTGTTCCGATTGTAGTCCGATTGTGTCTACGTGCTACTATCCGATTCTTCCGATTATCCGATTCTGTAAGTATCCGATTCCTCCGATTCATTCGTTCTGCTCCGATTTGCATGTCCGATTTCCGATTGCCCCTCCGATTAATTCCGATTGCTCCGATTTCCGATTCTCCGATTGCTCCGATTGTCGATTTCCGATTCTCCGATTCTCCGATTGAGGCCGCGTTCCGATTTCCGATTCTTACTTCCGATTGATTTTCCGATTTCCGATTATCCGATTAGAGCCGTCCGATTCACATCCGATTACTAGGTTGCTTCCGATTGTCCGATTCCATATCCGATTTCCGATTACTCCGATTTATCCGATTTGTCCGATTATCTCCTAGTTTTTCCGATTATCCGATTGCTCCGATTTCCGATTAGTCCGATTTTTTTCCGATTCATCCGATTTTGGTCCGATTGTTACTTCCGATTTCCGATTGATTTCCGATTAACTTCCGATTATGTATCCGATTTCCGATTTCCGATTGCGAGCTGTCCGATTCGTCCGATTTCCGATTCTGTCCGATTCTCCGATTGTAATGGCGAGTTCTCCGATTGGTCCGATTTCCGATTTCCGATTTGCATCTCCGATTTCCACTCTTCCGATTATTCCGATTTCCGATTTCCGATTTGGTCCGATTAATCCGATTTCCGATTTTTCCGATTTTCCGATTTTCCGATTCATTCCGATTTCCGATTCACTGAGAGATCCGATTATTCCGATTGCCTCCGATGTTCTTACCATCCGATTAACTTCCGATTTCCGATTGAATTACTCCGATTTCCGATTTATTCCTCCCTCCGATTTCCGATTTCCGATTTGGAAGCTCCGATTTTCCGATTTCCGATTACACTCCGATTCTATGTATGTCCGATTCAGTCCGATTCTGTCACGTCCGATTATCCGATTATCCGATTTTCCGATTTTATTAGAATCCGATTGTTCCGATTTGTCCGATTAAATTCCGATTCTGGTCGCCTCCGATTCATGCTCCGATTACAAATCCGATTGCGTCCGATTATCCGATTTCCGATTTGTACGTTTCCGATTTCCGATTTCCGATTATCCGATTTACTTTCCGATTTCCGATTATCGGCTTTCCGATTATTCCGATTCGAGCATCCGATTCTCGTCCGATTGTCTAAGAATCCGATTAGTTTTACGCCCTCCGATTTTCCGATTTCCGATTCTCCGATTTCCGATTTTCCGATTTTCCGATTATCCGATTTGTCCGATTCAGATTCCGATTAACTCCGATTGGCTAGTCCGATTGGTTCCGATTTTCCGATTATCCGATTGACCTCCGATTGATCCGATTGTCCGATTTCATCCGATTTGCATCCGATTTTCCGATTGAGTCCGATTAATCCGATTCTTATCCGATTTCCGATTTCCGATTTTCCGATTACTATGAACATCCGATTGAGCTGCCTGTTCCGATTGACCTATCCGATTGGCCTGCGGATCAGAGCATGTCCGATTGGGCGTAAGATCCGATTGTCCGATTTCCGATTTCCGATTATTATGTCCGATTAACTGCATCCGATTTCATCCGATTATCCGATTTCCGATTTGGTCCGATTTCCGATTCATCTCCGATTTCCGATTACCTCTCGTTCCGATTTTTTTCCGATTGATCTCCGATTGATCCGATTTTCCGATTGGTGTCCGATTTCCGATTTCCGATTGGTCCGATTCAGCAATCCGATTCTGAGGTCCGATTGTTTCCGATTAAGAGTCCGATTAGAAAATCCGATTTTTCCTTGAGTCCGATTTCCGATTATACATCCACGTCCGATTTCCGATTTCCGATTTCCTTCCGATTCTCCGATTCTCCGATTAGTGTTTTGATATAGTCCGATTTCCGATTATTTTCCGATTAGGGTAGTCCGATTGGTCCGATTCAGCTCCGATTATAGTTGTCCGATTGTTATTTGTCCGATTTTCCGATTGGCTTCCGATTTTCCGATTTCCGATTCGTCCGATTGCTCCGATTCTGGCCTGGTCCGATTTCCGATTGTTATCCGATTCCCATGTCCGATTAGCGTCCGATTTTACTTCCGATTACTCCGATTGGAGTCCGATTGAGGGGTTATCCGATTTCCGATTTCCGATTGGATCCGATTAGGTTCCGATTTCCGATTTACGGGATTCCGATTTTTCCGATTCTCCGATTTCCGATTAACTCCGATTTCCGATTTGCGCGCCTCCGATTGATCCGATTCTTCCGATTTGCATCCGATTGCCGTGTCCGATTGTCCGATTTCTCCGATTCCGAATCCGATTAATGTAAGTCCGATTCCCGGTCCGATTTCCACTCCGATTTCCGATTATCCGATTTGCTATAGATCGTCCGATTATCCGGTACGAATTCCGATTAACTCCGATTGCTCCGATTGTCCGATTCTCCGATTGTCCGATTGTCGTTTCCGATTTCCGATTTCCGATTTCCGATTATCCGATTTCCGATTATGTTCCGATTTCCGATTCGGCCTCCGATTTCCGATTTCCGATTGTCTCCGATTCGTTTCCGATTGTCCGATTCCGGTCGTAACTCCGATTTCCGATTGTCCGATTAGTTAACGACTCCGATTTCCGATTTCCGATTTTTTTCCGATTTTCCGATTGAGAGTCCGATTTCCGATT'
    pattern = 'TCCGATTTC'
    res = hm.approximate_pattern_matching(text=text, pattern=pattern, d=0)
    assert res == [3, 10, 17, 61, 68, 75, 95, 119, 126, 133, 161, 183, 200, 247, 262, 295, 355, 362, 432, 439, 446, 498,
                   563, 587, 621, 628, 643, 663, 701, 708, 745, 771, 812, 829, 876, 883, 959, 1044, 1062, 1069, 1106,
                   1197, 1235, 1260, 1284, 1304, 1416, 1432, 1469, 1493, 1500, 1507, 1542, 1610, 1660, 1730, 1778, 1848,
                   1876, 1883, 1969, 2045, 2097, 2121, 2196, 2220, 2229, 2272, 2292, 2299, 2306, 2461, 2480, 2514, 2552,
                   2693, 2709, 2738, 2841, 2867, 2898, 2922, 2947, 2992, 3014, 3021, 3118, 3137, 3187, 3273, 3303, 3323,
                   3363, 3420, 3435, 3442, 3497, 3544, 3551, 3558, 3622, 3629, 3655, 3662, 3683, 3698, 3738, 3851, 3858,
                   3877, 3903, 3960, 3967, 3990, 4053, 4060, 4067, 4115, 4172, 4224, 4304, 4337, 4377, 4424, 4450, 4609,
                   4639, 4674, 4761, 4781, 4938, 4945, 4952, 4964, 4992, 4999, 5018, 5126, 5136, 5164, 5230, 5323, 5330,
                   5404, 5431, 5532, 5558, 5573, 5599, 5724, 5762, 5833, 5865, 5943, 6014, 6074, 6115, 6122, 6153, 6208,
                   6215, 6235, 6259, 6266, 6292, 6334, 6404, 6425, 6449, 6456, 6485, 6671, 6693, 6700, 6727, 6825, 6840,
                   6969, 7028, 7035, 7150, 7157, 7191, 7209, 7226, 7244, 7317, 7324, 7418, 7443, 7450, 7457, 7506, 7614,
                   7655, 7745, 7752, 7780, 7819, 7836, 7908, 7956, 7968, 8078, 8085, 8092, 8107, 8125, 8144, 8151, 8205,
                   8237, 8244, 8282]

    text = 'CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT'
    pattern = 'ATTCTGGA'
    d = 3
    res = hm.approximate_pattern_matching(text=text, pattern=pattern, d=d)
    assert res == [6, 7, 26, 27]


def test_compute_frequencies():
    freq = hm.compute_frequencies('ACGCGGCTCTGAAA', k=2)
    assert freq == {'AC': 1, 'CG': 2, 'GC': 2, 'GG': 1, 'CT': 2, 'TC': 1, 'TG': 1, 'GA': 1, 'AA': 2}

    freq = hm.compute_frequencies('CGCCTAAATAGCCTCGCGGAGCCTTATGTCATACTCGTCCT', k=3)
    assert freq == {'CGC': 2, 'GCC': 3, 'CCT': 4, 'CTA': 1, 'TAA': 1, 'AAA': 1, 'AAT': 1, 'ATA': 2, 'TAG': 1, 'AGC': 2,
                    'CTC': 2, 'TCG': 2, 'GCG': 1, 'CGG': 1, 'GGA': 1, 'GAG': 1, 'CTT': 1, 'TTA': 1, 'TAT': 1, 'ATG': 1,
                    'TGT': 1, 'GTC': 2, 'TCA': 1, 'CAT': 1, 'TAC': 1, 'ACT': 1, 'CGT': 1, 'TCC': 1}


def test_find_clumps():
    genome = 'CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA'
    kmers = hm.find_clumps(genome=genome, k=5, l=50, t=4)
    assert sorted(kmers) == ['CGACA', 'GAAGA']

    genome, k, l, t = fetch_clumps_finding_input(Path('test/testcase01.txt'))
    kmers = hm.find_clumps(genome, k, l, t)
    assert kmers == ['AAACCAGGTGG']

    genome, k, l, t = fetch_clumps_finding_input(Path('test/testcase02.txt'))
    kmers = hm.find_clumps(genome, k, l, t)
    assert kmers == []

    genome, k, l, t = fetch_clumps_finding_input(Path('test/testcase03.txt'))
    kmers = hm.find_clumps(genome, k, l, t)
    assert sorted(kmers) == ['GAAGAACGG', 'GGGGGGGGG']


def test_hamming_distance():
    dist = hm.hamming_distance('GGGCCGTTGGT', 'GGACCGTTGAC')
    assert dist == 3

    dist = hm.hamming_distance('ACCGCTCGA', 'ACCGCTCGA')
    assert dist == 0

    dist = hm.hamming_distance('CTACAGCAATACGATCATATGCGGATCCGCAGTGGCCGGTAGACACACGT',
                               'CTACCCCGCTGCTCAATGACCGGGACTAAAGAGGCGAAGATTATGGTGTG')
    assert dist == 36


def test_approximate_pattern_count():
    text = 'TTTAGAGCCTTCAGAGG'
    pattern = 'GAGG'
    d = 2
    count = hm.approximate_pattern_count(text=text, pattern=pattern, d=d)
    assert count == 4

    text = 'TTGGAATGACAATTTGACTGTAGTTCCCGTGCGAACTAGCGCGTGGTAACAATCCCGTCGGGGGTAGCAGTTGTCTTATCTCAGCCACATTGCACACTGCCGGCCCGTGGAGGCTG'
    pattern = 'AACGAG'
    d = 3
    count = hm.approximate_pattern_count(text=text, pattern=pattern, d=d)
    assert count == 15

    text = 'CATGCCATTCGCATTGTCCCAGTGA'
    pattern = 'CCC'
    d = 2
    count = hm.approximate_pattern_count(text=text, pattern=pattern, d=d)
    assert count == 15


def test_frequent_words_with_complements():
    text = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
    k = 4
    d = 1
    kmers = hm.frequent_words_with_complements(text, k, d)
    assert sorted(kmers) == sorted(['ATGT', 'ACAT'])

    text = 'AAGCGGCGGACGTGGGGACAAGAAGGTGTCGGGTAAGCGGGGACACCGGACACGGCGGCGGGGGGCGGACACCGGCGGCGGAAGACGTACCGGGGACACGGAAGAAGGGGGCGGCGGGTACGGGGACGTGTGGAAGCGGACGGGTAAGGGGGGTGGGTGTAAGACCGGGGACCGGCGGGGAAGGGGGGTAAGGTACACAAGGGGGACGTAAGGTCGGGGGTCGGGGGGACACGGAAGACAAG'
    k = 5
    d = 2
    kmers = hm.frequent_words_with_complements(text, k, d)
    assert sorted(kmers) == ['CCCCC', 'GGGGG']
