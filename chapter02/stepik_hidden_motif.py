import hidden_motif

"""
A bunch of adapters to facilitate running the stepik challenges
"""

def MotifEnumeration(dna, k, d):
    return hidden_motif.motif_enumeration(dna, k, d)

def DistanceBetweenPatternAndStrings(pattern, dna):
    return hidden_motif.kmer_to_motif_distance(pattern, dna)

def MedianString(dna, k):
    return hidden_motif.median_string(dna, k)

def ProfileMostProbableKmer(text, k, profile):
    return hidden_motif.profile_most_probable_kmer()

def GreedyMotifSearch(dna, k, t):
    return hidden_motif.greedy_motif_search()

def GreedyMotifSearchWithPseudocounts(dna, k, t, pseudocount=1):
    return hidden_motif.greedy_motif_search_with_pseudocounts(dna, k, t, pseudocount)

