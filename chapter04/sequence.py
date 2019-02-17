from copy import deepcopy


def translate_RNA(rna):
    """
    Returns the protein translated from mRNA. Translation continues until a stop codon, or the end of the mRNA, whichever comes first.
    :param rna: The RNA to be translated, a string with length a multiple of 3.
    :return: The protein, a string.
    """

    def read_codons(polymer):
        assert len(polymer) % 3 == 0
        for i in range(0, len(polymer), 3):
            codon = polymer[i: i + 3]
            yield codon

    translation = {'AAA': 'K',
                   'AAC': 'N',
                   'AAG': 'K',
                   'AAU': 'N',
                   'ACA': 'T',
                   'ACC': 'T',
                   'ACG': 'T',
                   'ACU': 'T',
                   'AGA': 'R',
                   'AGC': 'S',
                   'AGG': 'R',
                   'AGU': 'S',
                   'AUA': 'I',
                   'AUC': 'I',
                   'AUG': 'M',
                   'AUU': 'I',
                   'CAA': 'Q',
                   'CAC': 'H',
                   'CAG': 'Q',
                   'CAU': 'H',
                   'CCA': 'P',
                   'CCC': 'P',
                   'CCG': 'P',
                   'CCU': 'P',
                   'CGA': 'R',
                   'CGC': 'R',
                   'CGG': 'R',
                   'CGU': 'R',
                   'CUA': 'L',
                   'CUC': 'L',
                   'CUG': 'L',
                   'CUU': 'L',
                   'GAA': 'E',
                   'GAC': 'D',
                   'GAG': 'E',
                   'GAU': 'D',
                   'GCA': 'A',
                   'GCC': 'A',
                   'GCG': 'A',
                   'GCU': 'A',
                   'GGA': 'G',
                   'GGC': 'G',
                   'GGG': 'G',
                   'GGU': 'G',
                   'GUA': 'V',
                   'GUC': 'V',
                   'GUG': 'V',
                   'GUU': 'V',
                   'UAC': 'Y',
                   'UAU': 'Y',
                   'UCA': 'S',
                   'UCC': 'S',
                   'UCG': 'S',
                   'UCU': 'S',
                   'UGC': 'C',
                   'UGG': 'W',
                   'UGU': 'C',
                   'UUA': 'L',
                   'UUC': 'F',
                   'UUG': 'L',
                   'UUU': 'F'}

    stop_codons = {'UAA', 'UAG', 'UGA'}

    protein = []
    for codon in read_codons(rna):
        if codon in stop_codons:
            break
        protein.append(translation[codon])
    return ''.join(protein)


def reverse_complement(dna):
    """
    Returns the reverese complement of a given DNA segment.
    :param dna: The DNA segment, a string.
    :return: The reverse complement, a string.
    """
    complements = {'A': 'T',
                   'T': 'A',
                   'G': 'C',
                   'C': 'G'}
    complement = [complements[nucleotide] for nucleotide in reversed(dna)]
    as_string = ''.join(complement)
    return as_string


def peptide_encoding(dna, protein):
    """
    Returns all segments within a piece of cDNA that encode a given protein. A segment is returned if it, or its reverse complement, encodes the protein.
    :param dna: The piece of cDNA, a string.
    :param protein: The given protein, a string.
    :return: A list of all found segments, a list of strings.
    """
    rc_dna = reverse_complement(dna)
    rna = dna.replace('T', 'U')
    rc_rna = rc_dna.replace('T', 'U')
    res = []
    len_encoding_rna = len(protein) * 3
    for i in range(0, len(dna) - len_encoding_rna):
        rna_segment = rna[i:i + len_encoding_rna]
        translation = translate_RNA(rna_segment)
        if translation == protein:
            res.append(dna[i:i + len_encoding_rna])
        rc_rna_segment = rc_rna[i:i + len_encoding_rna]
        rc_translation = translate_RNA(rc_rna_segment)
        if rc_translation == protein:
            res.append(dna[len(dna) - i - len_encoding_rna: len(dna) - i])
    return res


def get_ammino_mass():
    ammino_mass = {
        'G': 57,
        'A': 71,
        'S': 87,
        'P': 97,
        'V': 99,
        'T': 101,
        'C': 103,
        'I': 113,
        'L': 113,
        'N': 114,
        'D': 115,
        'K': 128,
        'Q': 128,
        'E': 129,
        'M': 131,
        'H': 137,
        'F': 147,
        'R': 156,
        'Y': 163,
        'W': 186
    }
    return ammino_mass


def get_ammino_mass_red():
    ammino_mass = get_ammino_mass()
    del ammino_mass['L']
    del ammino_mass['Q']
    return ammino_mass


def peptide_spectrum(peptide, cyclic = True):
    """
    Returns the theoretical peptide_spectrum for a peptide.
    :param peptide: The peptide, a string.
    :return: The peptide_spectrum, a list of numbers sorted in ascending order.
    """

    ammino_mass = get_ammino_mass()

    prefix_mass = [0]
    for ammino in peptide:
        prefix_mass.append(prefix_mass[-1] + ammino_mass[ammino])

    peptide_mass = prefix_mass[-1]
    spectrum = [0]
    for i in range(0, len(peptide)):
        for j in range(i + 1, len(peptide) + 1):
            spectrum.append(prefix_mass[j] - prefix_mass[i])
            if cyclic and i > 0 and j < len(peptide):
                spectrum.append(peptide_mass - (prefix_mass[j] - prefix_mass[i]))

    return list(sorted(spectrum))


def peptide_mass(peptide):
    ammino_mass = get_ammino_mass()
    mass = sum([ammino_mass[ammino] for ammino in peptide])
    return mass


from itertools import chain
from collections import Counter


def flatten(seq_of_seq):
    return list(chain(*seq_of_seq))


def sequence_cyclopeptide(spectrum):
    def expand(peptide):
        ammino_mass = get_ammino_mass_red()
        expanded = [peptide + ammino for ammino in ammino_mass.keys()]
        return expanded

    def consistent(peptide, spectrum):
        peptide_spect = peptide_spectrum(peptide, cyclic=False)
        peptide_count = Counter(peptide_spect)
        spectrum_count = Counter(spectrum)
        it_is_consistent = True if peptide_count & spectrum_count == peptide_count else False
        return it_is_consistent

    def consistent2(peptide, spectrum):
        ammino_mass = get_ammino_mass_red()
        peptide_masses = [ammino_mass[ammino] for ammino in peptide]
        peptide_count = Counter(peptide_masses)
        spectrum_count = Counter(spectrum)
        it_is_consistent = True if peptide_count & spectrum_count == peptide_count else False
        return it_is_consistent

    candidate_peptides = ['']
    final_peptides = []
    while candidate_peptides:
        candidate_peptides = flatten([expand(peptide) for peptide in candidate_peptides])
        new_candidates = []
        for peptide in candidate_peptides:
            the_peptide_spectrum = peptide_spectrum(peptide, cyclic=True)
            # if sum(the_peptide_spectrum) == sum(spectrum):
            if the_peptide_spectrum == spectrum:
                    final_peptides.append(peptide)
            elif consistent(peptide, spectrum):
                new_candidates.append(peptide)
        candidate_peptides = new_candidates

    return final_peptides
