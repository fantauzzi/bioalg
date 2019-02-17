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
            res.append(dna[len(dna)-i-len_encoding_rna: len(dna)-i])
    return res
