import hidden_motif


def immediate_neighbors(kmer):  # TODO should it be a generator?
    """
    All k-mers that dist no more than one from a given k-mer (by Hamming distance).
    :param kmer: The given k-mer.
    :return:A list of k-mers; it includes the one given as input parameter.
    """
    res = []
    lookup = {'G': ('T', 'A', 'C'),
              'T': ('G', 'A', 'C'),
              'A': ('G', 'T', 'C'),
              'C': ('G', 'T', 'A')}
    for i in range(0, len(kmer)):
        choices = lookup[kmer[i]]
        neighbors = [kmer[:i] + nucleotide + kmer[i + 1:] for nucleotide in choices]
        res += neighbors
    return res


def neighbors(kmer, d):
    """
    All k-mers that dist no more than a certain amount from a given one (by Hamming distance).
    :param kmer:The given k-mer.
    :param d: The distance not to be exceeded.
    :return: A list of k-mers, including teh given one.
    """
    kmers = {kmer}
    processed = set()
    for _ in range(1, d + 1):
        for item in kmers:
            if item in processed:
                continue
            imm_neighbors = immediate_neighbors(item)
            kmers = kmers.union(set(imm_neighbors))
            processed.add(item)
    return kmers


def motif_enumeration(dna, k, d):
    """
    Produce all (k, d)-motifs in a given DNA segment.
    :param dna: The given DNA segment.
    :param k: The length (in nucleotides) of each motif to be produced.
    :param d: The maximum allowed Hamming distance between a motif and a corresponding k-mer in the DNA segment.
    :return: A list of motifs, each a string of length k.
    """

    def is_kmer_in_dna(kmer, dna, d):
        k = len(kmer)
        for current_kmer in hidden_motif.kmers_from_dna(dna, k):
            if hidden_motif.hamming_distance(current_kmer, kmer) <= d:
                return True
        return False

    patterns = set()
    for i_dna, dna_piece in enumerate(dna):
        for kmer in hidden_motif.kmers_from_dna(dna_piece, k):
            kmer_neighbors = neighbors(kmer, d)
            for kmer_neighbor in kmer_neighbors:
                found = True
                for j_dna in range(0, len(dna)):
                    if not is_kmer_in_dna(kmer_neighbor, dna[j_dna], d):
                        found = False
                        break
                if found:
                    patterns.add(kmer_neighbor)
    return list(patterns)


"""
A bunch of adapters to facilitate running the stepik challenges
"""


def MotifEnumeration(dna, k, d):
    return motif_enumeration(dna, k, d)


def DistanceBetweenPatternAndStrings(pattern, dna):
    return hidden_motif.kmer_to_dna_distance(pattern, dna)


def MedianString(dna, k):
    median, _ = hidden_motif.median_string(dna, k)
    return median


def ProfileMostProbableKmer(text, k, profile):
    return hidden_motif.profile_most_probable_kmer(text, k, profile)


def GreedyMotifSearch(dna, k, t):
    assert t == len(dna)
    return hidden_motif.greedy_motifs_search(dna, k, pseudocount=0)


def GreedyMotifSearchWithPseudocounts(dna, k, t, pseudocount=1):
    assert t == len(dna)
    return hidden_motif.greedy_motifs_search(dna, k, pseudocount)
