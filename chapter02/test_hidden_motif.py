import hidden_motif
import stepik_hidden_motif
from stepik_hidden_motif import fetch_greedy_motif_input, fetch_kmer_to_dna_dist_input
import pytest
from math import isclose
from random import seed
import numpy as np
from pathlib import Path


def test_kmers_from_dna():
    dna = 'GATTACA'
    kmers = [nucleotide for nucleotide in hidden_motif.kmers_from_dna(dna, 1)]
    assert dna == ''.join(kmers)

    dna2 = 'GAAA'
    kmers2 = list(hidden_motif.kmers_from_dna(dna2, 4))
    assert kmers2 == [dna2]

    dna3 = 'GAGAC'
    kmers3 = list(hidden_motif.kmers_from_dna(dna3, 3))
    assert kmers3 == ['GAG', 'AGA', 'GAC']

    dna4 = 'TACG'
    res = []
    for item in hidden_motif.kmers_from_dna(dna4, 2):
        res.append(item)
    assert res == ['TA', 'AC', 'CG']


def test_motif_enumeration():
    k, d = 3, 1
    dna = ['ATTTGGC', 'TGCCTTA', 'CGGTATC', 'GAAAATT']
    res = stepik_hidden_motif.motif_enumeration(dna, k, d)
    assert sorted(res) == sorted(['ATA', 'ATT', 'GTT', 'TTT'])
    assert type(res) is list


def test_hamming_distance():
    dist1 = hidden_motif.hamming_distance('GATTACA', 'CATTAGA')
    assert dist1 == 2

    dist2 = hidden_motif.hamming_distance('GATTACA', 'GATTACA')
    assert dist2 == 0

    dist3 = hidden_motif.hamming_distance('G', 'C')
    assert dist3 == 1

    dist4 = hidden_motif.hamming_distance('', '')
    assert dist2 == 0


def test_generalised_hamming_distance():
    kmer = 'GATTACA'
    dna = 'ACGTGATTACACGTGATTACAA'
    dist1 = hidden_motif.generalised_hamming_distance(kmer, dna)
    assert dist1 == 0

    dist2 = hidden_motif.generalised_hamming_distance(dna, kmer)
    assert dist2 == 0

    dist3 = hidden_motif.generalised_hamming_distance(kmer, kmer)
    assert (dist3 == 0)

    dist4 = hidden_motif.generalised_hamming_distance('CGTCGT', 'A')
    assert dist4 == 1

    dist5 = hidden_motif.generalised_hamming_distance('CGTCGT', 'CA')
    assert dist5 == 1

    dist6 = hidden_motif.generalised_hamming_distance('CGTCGT', 'CAA')
    assert dist6 == 2

    dist7 = hidden_motif.generalised_hamming_distance('GATTCTCA', 'GCAAAGACGCTGACCAA')
    assert dist7 == 3


def test_kmer_to_dna_distance():
    kmer = 'GATTACA'
    motif = ['ACGTGATTACACGTGATTACAA', 'GATTACC', 'AATTACG']
    dist = hidden_motif.kmer_to_dna_distance(kmer, motif)
    assert dist == 3

    with pytest.raises(TypeError):
        _ = hidden_motif.kmer_to_dna_distance(kmer, 42)

    kmer = 'AAA'
    motif = ['TTACCTTAAC', 'GATATCTGTC', 'ACGGCGTTCG', 'CCCTAAAGAG', 'CGTCAGAGGT']
    dist = hidden_motif.kmer_to_dna_distance(kmer, motif)
    assert dist == 5

    kmer, motif = fetch_kmer_to_dna_dist_input(Path('test/testcase03.txt'))
    dist = hidden_motif.kmer_to_dna_distance(kmer, motif)
    assert dist == 76

    kmer, motif = fetch_kmer_to_dna_dist_input(Path('test/testcase04.txt'))
    dist = hidden_motif.kmer_to_dna_distance(kmer, motif)
    assert dist == 70


def test_flattened():
    flat1 = hidden_motif.flattened([('',), ('a', 'b', 'c'), [1, 2]])
    assert flat1 == ['', 'a', 'b', 'c', 1, 2]


def test_dna_to_number():
    for k in range(1, 8):
        for n in range(0, 4 ** k):
            kmer = hidden_motif.number_to_kmer(n, k)
            number = hidden_motif.kmer_to_number(kmer)
            assert number == n


def test_all_kmers():
    kmers2 = hidden_motif.all_kmers(1)
    assert sorted(kmers2) == sorted(['G', 'A', 'T', 'C'])

    kmers3 = hidden_motif.all_kmers(2)
    assert sorted(kmers3) == sorted(['AA', 'AT', 'AG', 'AC',
                                     'GA', 'GT', 'GG', 'GC',
                                     'TA', 'TT', 'TG', 'TC',
                                     'CA', 'CT', 'CG', 'CC'])

    kmers4 = hidden_motif.all_kmers(8)
    assert len(list(kmers4)) == 4 ** 8


def test_profile_most_probable_kmer():
    dna = 'ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT'
    k = 5

    profile = {
        'A': [0.2, 0.2, 0.3, 0.2, 0.3],
        'C': [0.4, 0.3, 0.1, 0.5, 0.1],
        'G': [0.3, 0.3, 0.5, 0.2, 0.4],
        'T': [0.1, 0.2, 0.1, 0.1, 0.2]}

    res = hidden_motif.profile_most_probable_kmer(dna, k, profile)
    assert res == 'CCGAG'

    dna = 'TGCCCGAGCTATCTTATGCGCATCGCATGCGGACCCTTCCCTAGGCTTGTCGCAAGCCATTATCCTGGGCGCTAGTTGCGCGAGTATTGTCAGACCTGATGACGCTGTAAGCTAGCGTGTTCAGCGGCGCGCAATGAGCGGTTTAGATCACAGAATCCTTTGGCGTATTCCTATCCGTTACATCACCTTCCTCACCCCTA'
    k = 6
    profile = {'A': [0.364, 0.333, 0.303, 0.212, 0.121, 0.242],
               'C': [0.182, 0.182, 0.212, 0.303, 0.182, 0.303],
               'G': [0.121, 0.303, 0.182, 0.273, 0.333, 0.303],
               'T': [0.333, 0.182, 0.303, 0.212, 0.364, 0.152]}
    res = hidden_motif.profile_most_probable_kmer(dna, k, profile)
    assert res == 'TGTCGC'

    dna = 'GCTGTGTAATCCTTGATGCCAATTCTTGGAGGTCCCCCGCGACTTACTCTATTGAGTGAGCAGCGACTTTTGTTTTGCCGAATATATCCACGCTAATTGTCTGTGGACGTCAAAAGAGATTGAGTTTACTTTTAGCACGAGAGTGCCCAGGTTAGATCGTGGGGGGACTGAACACTAGGTTCAGCCAGGTCGCAATGAAAGTTAGAATTACAACAATGGTGTAGGTGCCCACCCTTGATGGATCTAGACCGTTGTGTCCCGCCTCTGTTCTAAAGTGCATCGCCACCTCAGCTTCATTCACACGGTAGCGATGAATGGGCTGCCGAATCCCCTGTGAACGTACACAAGGGAACGCTGAAAGGGGTACCTCCTGACAATGAAACCGATATTGTCTATGAGCCGGATTGAACATACAGCTTTCTTTGATGGGGCTGGTTAATCTTCTTTTACTTCCTTCTAATACAATGTACATTTATCATTAAGTCTGGCAGCTTTAAACGACCTACCTATGAGTCGAGCCACAGGTTTTTAAAGCCTTCCTTGGCCTACCGTTGGGGCTTAACCCTGGGGGAGAACATATGACTCCCACATATCGAGGGACAGCAGATCTTAAACCGAAACGCCCGATCAGTATAAACAAATTAAGGGAACGACTGCTATTGGCGACGATTACAGAACTCATCTAGAGTCAATGCATGTGTCCCTCTATTGCGTACAACTTTCCGTGAGTTCGTCTAGTACCGGCCTTTGGAATTGGGTTCCAGCGAGGCTACCCGGCCGCTAATACTCTGGGTGAGTCTCTAAGCAGTCATGGGCAGCCGCTTATCGACATCTACCCGCGAAAGTACGCTCGGAGATAGGTCACCTGCCGTTGTCCGGCGCGGTTTTCTTTTAATCAAACACAGTGGGCGTCGACTACCAGGGGTTTTGCAGACCTCGCTACATGCCTCTGCCCTAATGTTTTGGCCCAGTGTGGTTCTTAAACAACCGGAACCTGAGA'
    k = 15
    profile = {
        'A': [0.364, 0.212, 0.318, 0.136, 0.288, 0.227, 0.212, 0.242, 0.303, 0.242, 0.318, 0.348, 0.242, 0.258, 0.121],
        'C': [0.212, 0.212, 0.182, 0.288, 0.182, 0.242, 0.273, 0.273, 0.258, 0.273, 0.212, 0.182, 0.242, 0.288, 0.379],
        'G': [0.227, 0.258, 0.273, 0.242, 0.167, 0.258, 0.394, 0.273, 0.167, 0.258, 0.258, 0.197, 0.242, 0.227, 0.273],
        'T': [0.197, 0.318, 0.227, 0.333, 0.364, 0.273, 0.121, 0.212, 0.273, 0.227, 0.212, 0.273, 0.273, 0.227, 0.227]}
    res = hidden_motif.profile_most_probable_kmer(dna, k, profile)
    assert res == 'ATTCTTGGAGGTCCC'


def test_median_string():
    dna = ['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA', 'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG', 'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
           'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC', 'AATCCACCAGCTCCACGTGCAATGTTGGCCTA']
    k = 6
    res, _ = hidden_motif.median_string(dna, k)
    assert res == 'AGGTGC'

    dna2 = ['AAACCGGTT', 'CCTTGGAAA', 'CTGAAACAA', 'TTTAAAGGG']
    k2 = 3
    res2, _ = hidden_motif.median_string(dna2, k2)
    assert res2 == 'AAA'

    dna = ['AAATTGACGCAT',
           'GACGACCACGTT',
           'CGTCAGCGCCTG',
           'GCTGAGCACCGG',
           'AGTTCGGGACAG']
    k = 3
    res, _ = hidden_motif.median_string(dna, k)
    assert res == 'GAC'

    dna = ['TGATGATAACGTGACGGGACTCAGCGGCGATGAAGGATGAGT',
           'CAGCGACAGACAATTTCAATAATATCCGCGGTAAGCGGCGTA',
           'TGCAGAGGTTGGTAACGCCGGCGACTCGGAGAGCTTTTCGCT',
           'TTTGTCATGAACTCAGATACCATAGAGCACCGGCGAGACTCA',
           'ACTGGGACTTCACATTAGGTTGAACCGCGAGCCAGGTGGGTG',
           'TTGCGGACGGGATACTCAATAACTAAGGTAGTTCAGCTGCGA',
           'TGGGAGGACACACATTTTCTTACCTCTTCCCAGCGAGATGGC',
           'GAAAAAACCTATAAAGTCCACTCTTTGCGGCGGCGAGCCATA',
           'CCACGTCCGTTACTCCGTCGCCGTCAGCGATAATGGGATGAG',
           'CCAAAGCTGCGAAATAACCATACTCTGCTCAGGAGCCCGATG']
    k = 6
    res, _ = hidden_motif.median_string(dna, k)
    assert res == 'CGGCGA'

    dna = ['TGCGTTCTGCAACTTCTTCTCACAACACTGGACTCGCGCTGT',
           'GACGGGTCGGCCTATGTTCAGTTCACACCGGGCTCAGAACGT',
           'AAAATCTATAAGGACGACGTAGTGACACTGGGCTCGGAGTCA',
           'TCTCCAAAGAGGCAGAAAGGGTTTTTCGCTACACCGGTGAGT',
           'CCCCCCATAGGACCATTAGTTACGTAGGCACGATAGACACAG',
           'CGACGTTTTATAAAGGAGACACTGTCCCTGTAATAATCTCCC',
           'ACACCGAGACAAAAGATTGAGGGAAGGCACCGTCCCCCTGAG',
           'GTCGAGGTTTTAGCAAAATCCTACCCGAGTACACCGTTAAGT',
           'AGAAATCGTACAAGCACTTTCATGCTAAGACAAGATACACTG',
           'ACACAGGAGTCAGTCCGAGGGCCGACGTCATATTCCTGTTAG']
    k = 6
    res, _ = hidden_motif.median_string(dna, k)
    assert res == 'ACACCG'


def test_motifs_profile():
    dna = ['GGCGTTCAGGCA', 'AAGAATCAGTCA', 'CAAGGAGTTCGC', 'CACGTCAATCAC', 'CAATAATATTCG']
    profile = hidden_motif.motifs_profile(dna, pseudocount=0)
    answer = {
        'A': [0.2, 0.8, 0.4, 0.2, 0.4, 0.4, 0.2, 0.8, 0, 0, 0.2, 0.4],
        'C': [0.6000000000000001, 0, 0.4, 0, 0, 0.2, 0.4, 0, 0, 0.4, 0.6000000000000001, 0.4],
        'G': [0.2, 0.2, 0.2, 0.6000000000000001, 0.2, 0, 0.2, 0, 0.4, 0.2, 0.2, 0.2],
        'T': [0, 0, 0, 0.2, 0.4, 0.4, 0.2, 0.2, 0.6000000000000001, 0.4, 0, 0]}
    for (key1, value1), (key2, value2) in zip(profile.items(), answer.items()):
        assert key1 == key2
        for item1, item2 in zip(value1, value2):
            assert isclose(item1, item2)


def test_score_motif():
    dna = ['GGCGTTCAGGCA', 'AAGAATCAGTCA', 'CAAGGAGTTCGC', 'CACGTCAATCAC', 'CAATAATATTCG']
    assert hidden_motif.score_motif(dna) == 28


def test_greedy_motifs_search():
    dna = ['GGCGTTCAGGCA', 'AAGAATCAGTCA', 'CAAGGAGTTCGC', 'CACGTCAATCAC', 'CAATAATATTCG']
    k = 3
    motif = hidden_motif.greedy_motifs_search(dna, k)
    assert sorted(motif) == sorted(['CAG', 'CAG', 'CAA', 'CAA', 'CAA'])

    k, dna = fetch_greedy_motif_input(Path('test/testcase01.txt'))
    motif = hidden_motif.greedy_motifs_search(dna, k)
    expected = ['AGTGGGTATCTC', 'TAAAAAGGTATA', 'AACCACGAGTAC', 'TGTCATGTGCGG', 'AACCTAAACCCT', 'AGTCGTTATCCC',
                'AGTAATATGTAC', 'AGTGGTTATCAC', 'AGTGGTTATCCC', 'AGTGGCTATCGC', 'AGTGGATATCCC', 'AGTGAGAAGCAA',
                'AGTGACTAGACA', 'TAAGACTAGTTA', 'TATGAAGGGTGA', 'AGTCGGGATAAC', 'AGTGGGTATCTC', 'AGCGGTTAGTCA',
                'AGTGAAATTCCT', 'TGTGGATGGCTT', 'TGTAGGTATCAC', 'TGCAGATATCCA', 'TGTGGTTATCAC', 'TGTCATTATTCA',
                'TGCGTAGATCAA']

    assert sorted(motif) == sorted(expected)
    k, dna = fetch_greedy_motif_input(Path('test/testcase02.txt'))
    motif = hidden_motif.greedy_motifs_search(dna, k)
    expected = ['CACCCTTGAGTT', 'ATTGGAGCATGA', 'TAATAAGTATTG', 'GGAGGCTCAATG', 'AGGCAGGGTCTT', 'AGGGATGCTGTG',
                'AACGAAGGAATG', 'CGTGAAGGAATA', 'AACCGAGCAAGA', 'AGACGGGCAGTG', 'CACTGGGAATGG', 'TGTCGTTGAATT',
                'CAATCAGAAATG', 'AATGAATCAATG', 'TGCGCAGAAATG', 'AATGACGCAATG', 'AGACCAGAAATG', 'CGCGCAGAAATG',
                'AGACCAGAAGTG', 'TATGCAGAAATG', 'ATCGGATAAATT', 'ATCCCAGAACTG', 'CAGGGTTAAAGG', 'TAACCAGAAGTG',
                'AAGTCAGAAGTG']
    assert sorted(motif) == sorted(expected)


def test_laplace_profile_matrix():
    dna = ['GGCGTTCAGGCA', 'AAGAATCAGTCA', 'CAAGGAGTTCGC', 'CACGTCAATCAC', 'CAATAATATTCG']
    profile = hidden_motif.motifs_profile(dna, pseudocount=1)
    answer = {'A': [0.2222222222222222, 0.5555555555555556, 0.3333333333333333, 0.2222222222222222, 0.3333333333333333,
                    0.3333333333333333, 0.2222222222222222, 0.5555555555555556, 0.1111111111111111, 0.1111111111111111,
                    0.2222222222222222, 0.3333333333333333],
              'C': [0.4444444444444444, 0.1111111111111111, 0.3333333333333333, 0.1111111111111111, 0.1111111111111111,
                    0.2222222222222222, 0.3333333333333333, 0.1111111111111111, 0.1111111111111111, 0.3333333333333333,
                    0.4444444444444444, 0.3333333333333333],
              'G': [0.2222222222222222, 0.2222222222222222, 0.2222222222222222, 0.4444444444444444, 0.2222222222222222,
                    0.1111111111111111, 0.2222222222222222, 0.1111111111111111, 0.3333333333333333, 0.2222222222222222,
                    0.2222222222222222, 0.2222222222222222],
              'T': [0.1111111111111111, 0.1111111111111111, 0.1111111111111111, 0.2222222222222222, 0.3333333333333333,
                    0.3333333333333333, 0.2222222222222222, 0.2222222222222222, 0.4444444444444444, 0.3333333333333333,
                    0.1111111111111111, 0.1111111111111111]}
    for (key1, value1), (key2, value2) in zip(profile.items(), answer.items()):
        assert key1 == key2
        for item1, item2 in zip(value1, value2):
            assert isclose(item1, item2)


def test_greedy_motif_search_with_pseudocounts():
    dna = ['GGCGTTCAGGCA', 'AAGAATCAGTCA', 'CAAGGAGTTCGC', 'CACGTCAATCAC', 'CAATAATATTCG']
    k = 3
    t = len(dna)
    pseudocount = 1
    motif = hidden_motif.greedy_motifs_search(dna=dna, k=k, pseudocount=pseudocount)
    assert sorted(motif) == sorted(['TTC', 'ATC', 'TTC', 'ATC', 'TTC'])

    with open('dataset_160_9.txt') as genome_file:
        genome = genome_file.readlines()
    genome = [str.rstrip(line) for line in genome]
    k = 12
    t = 25
    assert t == len(genome)
    result = ['GGTGCGTTAGTT', 'GGTGCGCGAGCT', 'GGTGCGTTAGTT', 'GGTGCGACAGTT', 'GGTGAGCAAGGT', 'GGTGAGTTAGCT',
              'GGTGGGGCAGGT', 'GGTGTGCAAGTT', 'GGTGCGACAGAT', 'GGTGTGATAGCT', 'GGTGGGAGAGTT', 'GGTGTGAAAGAT',
              'GGTGAGGAAGGT', 'GGTGAGCTAGGT', 'GGTGGGCAAGTT', 'GGTGAGCGAGGT', 'GGTGAGTAAGCT', 'GGTGCGACAGAT',
              'GGTGCGAAAGTT', 'GGTGTGTTAGTT', 'GGTGCGCCAGCT', 'GGTGTGCAAGTT', 'GGTGCGTTAGTT', 'GGTGAGTGAGAT',
              'GGTGTGATAGAT']
    motif = hidden_motif.greedy_motifs_search(dna=genome, k=k, pseudocount=1)
    assert sorted(motif) == sorted(result)


def test_mc_test_randomized_motif_search():
    seed(42)
    dna = ['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA', 'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG', 'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
           'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC', 'AATCCACCAGCTCCACGTGCAATGTTGGCCTA']
    k = 8
    motif = hidden_motif.mc_randomized_motif_search(dna, k, times=1000, seed=42)
    assert sorted(motif) == sorted(['CCAAGGTG', 'TACAGGCG', 'TCCACGTG', 'TCTCGGGG', 'TTCAGGTG'])


def test_gibbs_sampler():
    seed(42)
    dna = ['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA', 'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG', 'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
           'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC', 'AATCCACCAGCTCCACGTGCAATGTTGGCCTA']
    k = 8
    motif, _ = hidden_motif.gibbs_sampler(dna, k, n=100)
    # assert sorted(motif) == sorted(['AAACGGCC', 'CAAGTTTC', 'CACGTGCA', 'CGAGACCG', 'TAAGTGCC'])
    assert sorted(motif) == sorted(['AAACGGCC', 'ATACAGGC', 'CAAGGTGC', 'CAAGTTTC', 'CCACGTGC'])


def test_consensus_from_motifs():
    motifs = ['GGCGTTCAGGCA', 'AAGAATCAGTCA', 'CAAGGAGTTCGC', 'CACGTCAATCAC', 'CAATAATATTCG']
    consensus = hidden_motif.consensus_from_motifs(motifs)
    assert consensus == 'CACGTTCATTCA'


def test_sample_random_relative_entropy():
    dna = ['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA', 'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG', 'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
           'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC', 'AATCCACCAGCTCCACGTGCAATGTTGGCCTA']
    samples = hidden_motif.sample_random_relative_entropy(dna=dna, k=7, n=5, seed=42)
    assert np.isclose(samples, [-4.361310120017858, -3.8374744768060447, -3.3695044762291197, -3.3698300963604413,
                                -3.237474476806044]).all()
