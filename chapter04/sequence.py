from copy import deepcopy
from itertools import chain
from collections import Counter


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
    """
    Returns an association between amminoacids and their mass.
    :return: A dictionary, associating the ammino acid one letter code in the key, a string, with its mass, an integer number.
    """
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


def get_reverse_ammino_mass():
    ammino_mass = get_ammino_mass()
    reversed = {value: key for key, value in ammino_mass.items()}
    reversed[113] = 'I/L'
    reversed[128] = 'K/Q'
    return reversed


def get_ammino_mass_red():
    ammino_mass = get_ammino_mass()
    del ammino_mass['L']
    del ammino_mass['Q']
    return ammino_mass


def peptide_spectrum(peptide, cyclic=True):
    """
    Returns the theoretical peptide_spectrum for a peptide.
    :param peptide: The peptide, a string.
    :return: The peptide spectrum, a list of numbers sorted in ascending order.
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
    """
    Returns the overall mass for a peptide, computed from its teoretical spectrum.
    :param peptide: The peptide, a string.
    :return: The mass, an integer.
    """
    ammino_mass = get_ammino_mass()
    mass = sum([ammino_mass[ammino] for ammino in peptide])
    return mass


def flatten(seq_of_seq):
    """
    Returns a sequence obtained by flattening a sequence of sequences. E.g. [[1, 2, 3], [4, 5]] becomes [1, 2, 3, 4, 5]
    :param seq_of_seq: The sequence of sequences.
    :return: The sequence after flattening, a list.
    """
    return list(chain(*seq_of_seq))


def expand(peptide):
    """
    Expands a given peptide into the list of peptides that have it for prefix,
    :param peptide: The given peptide, a string.
    :return: The list of peptides, a list of strings.
    """
    ammino_mass = get_ammino_mass_red()
    expanded = [peptide + ammino for ammino in ammino_mass.keys()]
    return expanded


def sequence_cyclopeptide(spectrum):
    """
    Returns all cyclopeptides with a given mass spectrum. The implementation proceeds by branch and bound: it generates candidate peptides of increasing length, discarding those that cannot be, based on their spectrum, an infix of the peptide whose spectrum is given; iterations continue until all circular peptides with the given spectrum have been generated.
    :param spectrum: The given spectrum, a list of integers in non-decreasing order; the same value may be repeated in the list multiple times as needed.
    :return: The resulting cyclopeptides, a list of lists; each element in the list is a cyclopeptide, which is represented as the list of atomic weights of its component ammino acids (a list of integers).
    """

    def consistent(peptide, spectrum):
        """
        Checks and returns if the linear spectrum of a given peptide is a subset of a certain spectrum.
        :param peptide: The given peptide, a string.
        :param spectrum: The spectrum to be checked against, a list of integers in non-decreasing order.
        :return: True if and only if the checked condition is met, False otherwise.
        """
        peptide_spect = peptide_spectrum(peptide, cyclic=False)
        peptide_count = Counter(peptide_spect)
        spectrum_count = Counter(spectrum)
        it_is_consistent = True if peptide_count & spectrum_count == peptide_count else False
        return it_is_consistent

    # List of candidate peptides, to be grown
    candidate_peptides = ['']
    # List of peptides found to match the given spectrum
    final_peptides = []
    while candidate_peptides:
        # Expand the candidate peptides; expanded candidates are all the peptides with current candidates as prefix
        candidate_peptides = flatten([expand(peptide) for peptide in candidate_peptides])
        # List of candidates to be further explored at the next iteration
        new_candidates = []
        for peptide in candidate_peptides:
            # If the candidate peptide has a spectrum matching the given spectrum, then list it in the solution
            if peptide_spectrum(peptide, cyclic=True) == spectrum:
                final_peptides.append(peptide)
            elif consistent(peptide, spectrum):  # Otherwise, could it still be expanded into a solution?
                new_candidates.append(peptide)  # If yes, list it to be expanded at the next iteration
        candidate_peptides = new_candidates

    ammino_mass = get_ammino_mass_red()
    masses = [[ammino_mass[ammino] for ammino in peptide] for peptide in final_peptides]

    return masses


def score_peptide(peptide, spectrum, cyclic=True):
    """
    Return the score for a peptide, linear or cyclic, against a given spectrum. The score is given by the number of masses that belong to either the peptide spectrum or the given spectrum but not both.
    :param peptide: The peptide, a string.
    :param spectrum: The spectrum, a list of integer numbers in non-descending order.
    :param cyclic: True to indicate that the peptide is cyclic, False to indicate it is linear.
    :return: The score, an integer.
    """
    pept_spectr = peptide_spectrum(peptide, cyclic=cyclic)
    pept_count = Counter(pept_spectr)
    spectr_count = Counter(spectrum)
    # Take the max (union) of the two Counters
    totals_count = pept_count | spectr_count
    total = sum([value for value in totals_count.values()])
    diff_count = deepcopy(pept_count)
    diff_count.subtract(spectr_count)
    total_diffs = sum([abs(value) for value in diff_count.values()])
    score = total - total_diffs
    return score


def leaderboard_peptide_sequence(spectrum, n):
    def trim(scored_peptides, n):
        trimmed = Counter({})
        best_n = scored_peptides.most_common(n)
        lowest_passable = min([value for (_, value) in best_n])
        for peptide, score in scored_peptides.items():
            if score >= lowest_passable:
                trimmed[peptide] = score
        return trimmed

    parent_mass = max(spectrum)
    leader_peptide = ''
    leader_score = score_peptide(leader_peptide, spectrum, cyclic=False)
    leaderboard = Counter({'': leader_score})
    while leaderboard:
        expanded_peptides = flatten(
            [expand(peptide) for peptide in leaderboard.keys()])  # TODO move into its own function as common code
        leaderboard = Counter(
            {peptide: score_peptide(peptide, spectrum, cyclic=False) for peptide in expanded_peptides})
        new_leaderboard = Counter({})
        for peptide, peptide_score in leaderboard.items():
            if peptide_mass(peptide) == parent_mass:
                new_leaderboard[peptide] = peptide_score
                if peptide_score > leader_score:
                    leader_score = peptide_score
                    leader_peptide = peptide
            elif peptide_mass(peptide) < parent_mass:
                new_leaderboard[peptide] = peptide_score
        leaderboard = trim(new_leaderboard, n) if new_leaderboard else new_leaderboard

    ammino_mass = get_ammino_mass_red()
    masses = [ammino_mass[ammino] for ammino in leader_peptide]

    return masses


def spectral_convolution(spectrum):
    res = [spectrum[j] - spectrum[i] for j in range(0, len(spectrum)) for i in range(0, j)]
    zeros_stripped = []
    for item in res:
        if item > 0:
            zeros_stripped.append(item)
    return zeros_stripped
