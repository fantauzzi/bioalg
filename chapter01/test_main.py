from main import compute_skew, argsmin, DNA_complement
from random import choices, seed

def test_skew():
    genome = 'CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAACTAGTTCACACGGCTTGATGGCAAATGGTTTTTCCGGCGACCGTAATCGTCCACCGAG'
    running_skew = [0] + compute_skew(genome)
    mins = argsmin(running_skew)
    assert mins == [53, 97]

    with open('skew_dataset.txt') as genome_file:
        genome = genome_file.read()
    running_skew = [0] + compute_skew(genome)
    mins = argsmin(running_skew)
    assert sorted(mins) == sorted([89969, 89970, 89971, 90345, 90346])

def test_DNA_Complement():
    res = DNA_complement('GATTACA')
    assert res == 'TGTAATC'

    seed(42)
    dna = ''.join(choices(['A', 'T', 'C', 'G'], k=100))
    compl = DNA_complement(dna)
    assert dna == DNA_complement(compl)