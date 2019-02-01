import functools
import operator
from collections import Counter
from random import randint, choices
from copy import deepcopy
from random import seed


def hamming_distance(a, b):
    """
    Compute the Hamming distance between two strings of the same length (case sensitive).
    :param a: The first string.
    :param b: The second string.
    :return: The Hamming distance between the two given strings.
    """
    assert len(a) == len(b)
    dist = sum(item_a != item_b for item_a, item_b in zip(a, b))
    return dist


def kmers_from_dna(dna, k):
    """
    Generate the list of k-mers in DNA segment, one k-mer at a time, in the same order as they appear in the segment.
    :param dna: The DNA segment.
    :param k: The k-mer size (number of nucleotides).
    :return: Produces on k-mer at a time.
    """
    assert k >= 1
    assert len(dna) >= k

    assert len(dna) >= k
    for i in range(0, len(dna) - k + 1):
        kmer = dna[i:i + k]
        yield kmer


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
    kmers = set((kmer,))
    processed = set()
    for _ in range(1, d + 1):
        for item in kmers:
            if item in processed:
                continue
            imm_neighbors = immediate_neighbors(item)
            kmers = kmers.union(set(imm_neighbors))
            processed.add(item)
    return kmers


"""     Input: Integers k and d, followed by a collection of strings Dna.
     Output: All (k, d)-motifs in Dna.
     """


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
        for current_kmer in kmers_from_dna(dna, k):
            if hamming_distance(current_kmer, kmer) <= d:
                return True
        return False

    patterns = set()
    for i_dna, dna_piece in enumerate(dna):
        for kmer in kmers_from_dna(dna_piece, k):
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


def generalised_hamming_distance(dna1, dna2):
    """
    The distance between two DNA segments not necessarily of the same length. If the two segments are of the same length, then it is the same as the Hamming distance.
    :param dna1: A DNA segment.
    :param dna2: A DNA segment.
    :return: The distance, and integer number.
    """
    if len(dna1) == len(dna2):
        return hamming_distance(dna1, dna2)
    if len(dna1) > len(dna2):
        dna = dna1
        kmer = dna2
    else:
        dna = dna2
        kmer = dna1
    k = len(kmer)

    dist = min([hamming_distance(kmer, kmer2) for kmer2 in kmers_from_dna(dna, k)])
    return dist


def kmer_to_motif_distance(kmer, motif):
    """
    The distance between a k-mer and a motif (i.e. a collection of DNA segments). It is the minimum distance between the k-mer and each segment in DNA.
    :param kmer: The k-mer.
    :param motif: The motif, an iterable container.
    :return: The distance, an integer number.
    :raises: TypeError, if motif is not an iterable container
    """
    try:
        iter(motif)
    except TypeError:
        raise
    dist = sum([generalised_hamming_distance(kmer, dna) for dna in motif])
    return dist


def flattened(list_of_lists):
    """
    Flatten a list of lists into a list. E.g., [[], [1], [1, 3]] is flattened into [1, 1, 3].
    :param list_of_lists: The list to be flattened.
    :return: The flattened list.
    """
    res = functools.reduce(operator.iconcat, list_of_lists, [])
    return res


def dna_to_number(dna):
    """
    Convert a DNA string into an integer number, such that the conversion and its inversion make a bijection. The inversion of the conversion in implemented by number_to_kmer().
    :param dna: The DNA string.
    :return: The corresponding integer number.
    """
    nucleotide_to_number = {'A': 0, 'C': 1, 'G': 2, 'T': 3}  # If you change this, also change it in number_to_kmer()
    if dna == '':
        return 0
    last_nucleotide = dna[-1]
    prefix = dna[:len(dna) - 1]
    res = 4 * dna_to_number(prefix) + nucleotide_to_number[last_nucleotide]
    return res


def number_to_kmer(n, k):
    """
    Convert a number into a k-mer, such that the conversion and its inversion make a bijection. The inversion of the conversion in implemented by dna_to_number().
    :param n: The integer number to be converted.
    :param k: The length of the desired k-mer.
    :return: A k-mer.
    """
    assert k >= 1
    number_to_nucleotide = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}  # If you change this, also change it in dna_to_number()
    if k == 1:
        return number_to_nucleotide[n]
    prefix_number = n // 4
    r = n % 4
    nucleotide = number_to_nucleotide[r]
    prefix_dna = number_to_kmer(prefix_number, k - 1)
    res = prefix_dna + nucleotide
    return res


def all_kmers(k):
    """
    Generate all k-mers with a given length.
    :param k: The desired k-mer length.
    :return: One k-mer at a time, all the k-mers.
    """
    for i in range(0, 4 ** k):
        res = number_to_kmer(i, k)
        yield res


def median_string(dna, k):  # TODO add unit test
    """
    Among all possible k-mers, find the one with minimum distance from a given collection of DNA segments. If multiple k-mers satisfy the requirement, choose any of them.
    :param dna: The collection of DNA segments (each segment is a string).
    :param k: The length of the k-mer to be found.
    :return: A k-mer, a string.
    """
    min_dist = float('inf')
    best_kmer = None
    for kmer in all_kmers(k):
        dist = kmer_to_motif_distance(kmer, dna)
        if dist < min_dist:
            min_dist = dist
            best_kmer = kmer
    return best_kmer


def profile_most_probable_kmer(text, k, profile):
    """
    Choose the k-mer within a DNA segment that is most probable given a certain probability profile.
    :param text: The DNA segment, a string.
    :param k: The size of the wanted k-mer.
    :param profile: The probability profile; a dictionary with keys 'A', 'C', 'G' and 'T' which associates every nucleotide with a list of probabilities of length k.
    :return: The most probable k-mer, a string.
    """

    most_prob_kmer = None
    highest_prob = -1
    for kmer in kmers_from_dna(text, k):
        prob = functools.reduce(operator.mul,
                                [profile[nucleotide][i] for nucleotide, i in zip(kmer, range(0, len(kmer)))], 1)
        if prob > highest_prob:
            highest_prob = prob
            most_prob_kmer = kmer
    return most_prob_kmer


def profile_matrix(motif):
    """
    Determine the probability profile matrix for a motif.
    :param motif: The given motif.
    :return: The probability matrix, a dictionary with keys 'A', 'C', 'G' and 'T' which associates every nucleotide with a list of probabilities of length k.
    """
    n_rows = len(motif)
    n_cols = len(motif[0])
    profile = {'A': [0] * n_cols, 'C': [0] * n_cols, 'G': [0] * n_cols, 'T': [0] * n_cols}
    for column in range(0, n_cols):
        for row in range(0, n_rows):
            nucleotide = motif[row][column]
            profile[nucleotide][column] += 1 / n_rows
    return profile


def laplace_profile_matrix(motif, pseudocount=1):
    """
    Determine the probability profile matrix for a motif using Laplace's Rule of Succession.
    :param motif: The given motif.
    :param pseudocount: The constant value added to the count of every nucleotide in every position, before normalising the counts into probabilities.
    :return: The probability matrix, a dictionary with keys 'A', 'C', 'G' and 'T' which associates every nucleotide with a list of probabilities of length k.
    """
    n_rows = len(motif)
    n_cols = len(motif[0])
    profile = {'A': [pseudocount] * n_cols, 'C': [pseudocount] * n_cols, 'G': [pseudocount] * n_cols,
               'T': [pseudocount] * n_cols}
    for column in range(0, n_cols):
        for row in range(0, n_rows):
            nucleotide = motif[row][column]
            profile[nucleotide][column] += 1
    # This sucks, redo it with pandas or numpy
    for col in range(0, n_cols):
        col_total = profile['A'][col] + profile['C'][col] + profile['G'][col] + profile['T'][col]
        profile['A'][col] /= col_total
        profile['C'][col] /= col_total
        profile['G'][col] /= col_total
        profile['T'][col] /= col_total
    return profile


def score_motif(motif):
    """
    Determine the score of a motif, i.e. a set of k-mers all of the same length. The more dissimilar the k-mers are in every position (from 0 to k-1), the higher the score. A score of 0 indicates that all k-mers in the motif are identical.
    :param motif: The motif to be scored.
    :return: The score, an integer number.
    """
    n_rows = len(motif)
    n_cols = len(motif[0])
    scores = []
    for col in range(0, n_cols):
        # ps = pandas.Series([motif[row][col] for row in range(0, n_rows)])
        # counts = ps.value_counts()
        counts = Counter([motif[row][col] for row in range(0, n_rows)])
        scores.append(n_rows - max(counts.values()))
    total_score = sum(scores)
    return total_score


# Your function should return a list of strings.
def greedy_motif_search(dna, k,
                        t):  # TODO make the function with snake_case signature and omitting the redundant t parameter
    """
    Find the motifs that best fit a collection of DNA segments (have the lowest score), using a greedy search. If multiple motifs have the same score, then provide just one of them.
    :param dna: The collection of DNA segments (strings), not necessarily of the same length.
    :param k: The size of motifs to be discovered; the same for every motif, an integer number.
    :param t: The number of motifs to be found, must be the same as the number of segments in dna, an integer number.
    :return: The discovered motifs, a list of strings.
    """
    best_score = float('inf')
    best_motif = None
    for kmer in kmers_from_dna(dna[0], k):
        motif = [kmer]
        for i in range(1, t):
            profile = profile_matrix(motif)
            most_prob_kmer = profile_most_probable_kmer(dna[i], k, profile)
            motif.append(most_prob_kmer)
        score = score_motif(motif)
        if score < best_score:
            best_motif = motif
            best_score = score
    return best_motif


def greedy_motif_search_with_pseudocounts(dna, k, t,
                                          pseudocount=1):  # TODO this can be the same as the function above, when pseudocounts=0
    best_score = float('inf')
    best_motif = None
    for kmer in kmers_from_dna(dna[0], k):
        motif = [kmer]
        for i in range(1, t):
            profile = laplace_profile_matrix(motif, pseudocount)
            most_prob_kmer = profile_most_probable_kmer(dna[i], k, profile)
            motif.append(most_prob_kmer)
        score = score_motif(motif)
        if score < best_score:
            best_motif = motif  # TODO Should I do a deep copy here?
            best_score = score
    return best_motif


def randomized_motif_search(dna, k):
    """
    Find the motifs that best fit a collection of DNA segments (have the lowest score), using a randomised search. Setting the same seed with random.seed() before calling the function, will produce the same result every time.
    :param dna: The collection of DNA segments (strings), not necessarily of the same length.
    :param k: The size of motifs to be discovered; the same for every motif, an integer number.
    :return: The discovered motifs, a list of strings.
    """
    seed(42)  # TODO remove from here. Make it a parameter.
    t = len(dna)
    s = len(dna[0])
    assert s >= k
    motif = [dna[row][i:i + k] for row, i in zip(range(0, t), [randint(0, s - k - 1) for _ in range(0, t)])]
    best_motif = deepcopy(motif)
    best_motif_score = score_motif(best_motif)
    while True:
        profile = laplace_profile_matrix(motif)
        motif = [profile_most_probable_kmer(text, k, profile) for text in dna]
        score = score_motif(motif)
        if score < best_motif_score:
            best_motif = deepcopy(motif)
            best_motif_score = score
        else:
            return best_motif, best_motif_score


def mc_randomized_motif_search(dna, k, times):
    """
    Repeat the randomised motif search in DNA segments a given number of times, and select the best result. Setting the same seed with random.seed() before calling the function, will produce the same result every time.
    :param dna: The collection of DNA segments (strings), not necessarily of the same length.
    :param k: The size of motifs to be discovered; the same for every motif, an integer number.
    :param times: The number of time the motifs search has to be repeated
    :return: The discovered motifs, a list of strings.
    """
    best_score = float('inf')
    best_motif = None
    for _ in range(0, times):
        motif, score = randomized_motif_search(dna, k)
        if score < best_score:
            best_score = score
            best_motif = deepcopy(motif)
    return best_motif


def gibbs_sampler(dna, k, n):
    """
    Find the motifs that best fit a collection of DNA segments (have the lowest score), using a Gibbs sampler. Setting the same seed with random.seed() before calling the function, will produce the same result every time.
    :param dna: The collection of DNA segments (strings), not necessarily of the same length.
    :param k: The size of motifs to be discovered; the same for every motif, an integer number.
    :param n: The number of iterations to be performed, an integer number.
    :return: The discovered motifs, a list of strings.
    """
    t = len(dna)  # No. of strings in dna
    s = len(dna[0])  # Length of every string in dna
    assert s >= k
    # Randomly select one motif per dna string, each of length  k, as a starting set of motifs
    motifs = [dna[row][i:i + k] for row, i in zip(range(0, t), [randint(0, s - k - 1) for _ in range(0, t)])]
    best_motifs = deepcopy(motifs)
    best_motif_score = score_motif(best_motifs)
    for _ in range(0, n):  # TODO find better stop criteria
        # Randomly choose one of the motifs
        i = randint(0, t - 1)
        # Find the probability profile of the motifs, without the randomly choosen one
        motif_ex_i = motifs[:i] + motifs[i + 1:]
        profile = laplace_profile_matrix(motif_ex_i)
        # Determine the likelihood of every k-mer in the i-th string of dna, based on the just calculated profile
        # proportions = [profile[nucleotide][col] for nucleotide, col in zip(dna[i], range(0, s))]
        proportions = []
        for kmer in kmers_from_dna(dna[i], k):
            prop = functools.reduce(operator.mul,
                                    [profile[nucleotide][col] for nucleotide, col in zip(kmer, range(0, len(kmer)))], 1)
            proportions.append(prop)
        sum_proportions = sum(proportions)
        prob_distr = [prop / sum_proportions for prop in proportions]
        assert len(prob_distr) == s - k + 1
        # Sample one kmer from the i-th string in dna, based on the just calculated probability distribution
        drafted_kmer_i = choices(range(0, len(prob_distr)), weights=prob_distr)[0]
        # Replace the i-th motif in the current set of motifs with the sampled one
        motifs[i] = dna[i][drafted_kmer_i:drafted_kmer_i + k]
        # Check if the obtaiend motif is the best so far
        score = score_motif(motifs)
        if score < best_motif_score:
            best_motif_score = score
            best_motifs = deepcopy(motifs)
    return best_motifs


def main_for_mc_randomized_motif_search():
    s1 = input()
    k, t = s1.split(sep=' ')
    k, t = int(k), int(t)
    dna = []
    for _ in range(0, t):
        s2 = input()
        dna.append(s2)

    res = mc_randomized_motif_search(dna, k, times=1000)

    for item in res:
        print(item)


def main_for_gibbs_sampler():
    seed(2)
    s1 = input()
    k, t, n = s1.split(sep=' ')
    k, t, n = int(k), int(t), int(n)
    dna = []
    for _ in range(0, t):
        s2 = input()
        dna.append(s2)
    assert len(dna) == t

    res = gibbs_sampler(dna, k, n)

    for item in res:
        print(item)


if __name__ == '__main__':
    pass
