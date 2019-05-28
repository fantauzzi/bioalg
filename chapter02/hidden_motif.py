import functools
import operator
from collections import Counter
from random import randint, choices, shuffle
from copy import deepcopy
import random
import matplotlib.pyplot as plt
import matplotlib
from tqdm import tqdm
from Bio import SeqIO
from numpy import argmin
from math import log2, isclose

plt.ion()

"""
Transcription factors are proteins that regulate the expression of genes, and turn them on and off. They do so by binding to specific short segments of DNA, tipically upstream of genes, called regulatory motifs. Here we try to identify regulatory motifs via computational analysis.

We assume that the regions upstream of the genes of interest have been sequenced, for around a thousand nucleotides, and we search them for recurring, short patterns, in the order of 10-20 nucleotides. The search is complicated by the fact that motifs upstream of different genes have variations between each other, in spite of the genese being related in function.

Let's begin to formalise the problem. The sequenced regions are strings of nucleotides of length of up to 1000. The motifs we are looking for are k-mers, that is strings of nucleotides of length exactly k.

To match k-mers with each other and then with longer segments, we define a distance. The choice falls on the Hamming distance: given two strings of the same length, the Hamming distance is the number of respective positions where the strings differ. For instance, the Hamming distance between 'ACGTT' and'ATGAT' is 2.     
"""


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


"""
We need to generalise the above definition to strings of different length. To facilitate that, we first write a generator that, given a DNA segment, produces all the k-mers in it, for an assigned value of k. For example, all 3-mers in GATTACA are: GAT, ATT, TTA, TAC, ACA.
"""


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


"""
We are ready to define the distance between strings of different length. Given a k-mer and a string of length greater than k, we take the *minimum* of the distances between the given k-mer and every k-mer from the longer string.
"""


def generalised_hamming_distance(a, b):
    """
    The distance between two string, not necessarily of the same length. If the two strings are of the same length, then it is the same as the Hamming distance.
    :param a: A string.
    :param b: Another string.
    :return: The distance, an integer number.
    """
    if len(a) == len(b):
        return hamming_distance(a, b)
    if len(a) > len(b):
        dna = a
        kmer = b
    else:
        dna = b
        kmer = a
    k = len(kmer)

    dist = min([hamming_distance(kmer, kmer2) for kmer2 in kmers_from_dna(dna, k)])
    return dist


"""
Still working toward finding motifs in a collection of DNA segments. Those motifs tend to have commonalities (same nucleotides in the same positions) but also variations. A possible strategy could be to set a length k, then compare all k-mers of every DNA segment with every k-mer of every other DNA segment, looking for motifs that are within a set distance d. Problem with this strategy is, it doesn't work. Beside being computationally unsustainable, to find the motifs it requires us to relax the value of d to the point that also many false positives are detected.

This is what we are going to do instead: we put together a string of nucleotides that is representative of every motif. We will set a length k, then look for a k-mer that doesn't necessarily appear in any of the DNA segments, but that minimises the sum of its distances from the DNA segments. Such k-mer is sort of "blueprint" for the motifs we are looking for. It is called, the *median string* for the set of DNA segments.

The algorithm is "brute force", where we enumerate all possible k-mers, and compute their distances from every DNA segment.
"""


def kmer_to_dna_distance(kmer, dna):
    """
    The distance between a k-mer and a collection of DNA segments. It is the sum of the distances between the k-mer and each segment.
    :param kmer: The k-mer.
    :param dna: The motif, an iterable container.
    :return: The distance, an integer number.
    :raises: TypeError, if motif is not an iterable container
    """
    try:
        iter(dna)
    except TypeError:
        raise
    dist = sum([generalised_hamming_distance(kmer, dna) for dna in dna])
    return dist


def flattened(list_of_lists):
    """
    Flatten a list of lists into a list. E.g., [[], [1], [1, 3]] is flattened into [1, 1, 3].
    :param list_of_lists: The list to be flattened.
    :return: The flattened list.
    """
    res = functools.reduce(operator.iconcat, list_of_lists, [])
    return res


"""
To enumerate al k-mers for a given k, we establish a bijection between those k-mers and an interval of natural numbers. Then, enumerating all k-mers will be as easy as counting natural numbers in that interval, and mapping each number to the corresponding k-mer.
"""


def nucleotide_numbering():
    """
    Establish a mapping between the four DNA nucleotides and numbers.
    :return: a pair of dictionaries, the first one mapping nucleotides (strings) to numbers, the second mapping numbers to nucleotides.
    """
    nucleotide_to_number = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    number_to_nucleotide = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    return nucleotide_to_number, number_to_nucleotide


def kmer_to_number(kmer):
    """
    Convert a DNA segment into an integer number, such that the conversion and its inverse make a bijection. The inverse is implemented by number_to_kmer().
    :param kmer: The DNA segment.
    :return: The corresponding integer number.
    """
    nucleotide_to_number, _ = nucleotide_numbering()
    if kmer == '':
        return 0
    last_nucleotide = kmer[-1]
    prefix = kmer[:len(kmer) - 1]
    res = 4 * kmer_to_number(prefix) + nucleotide_to_number[last_nucleotide]
    return res


def number_to_kmer(n, k):
    """
    Convert a number into a k-mer, such that the conversion and its inversion make a bijection. The inversion of the conversion in implemented by kmer_to_number().
    :param n: The integer number to be converted.
    :param k: The length of the desired k-mer.
    :return: A k-mer.
    """
    assert k >= 1
    _, number_to_nucleotide = nucleotide_numbering()
    if k == 1:
        return number_to_nucleotide[n]
    prefix_number = n // 4
    r = n % 4
    nucleotide = number_to_nucleotide[r]
    prefix_dna = number_to_kmer(prefix_number, k - 1)
    res = prefix_dna + nucleotide
    return res


"""
Now ready to implement a generator for all possible k-mers, for a given k.
"""


def all_kmers(k):
    """
    Generate all k-mers with a given length.
    :param k: The desired k-mer length.
    :return: One k-mer at a time, all the k-mers.
    """
    for i in range(0, 4 ** k):
        res = number_to_kmer(i, k)
        yield res


"""
And finally, the calculation of the median string of length k.
"""


def median_string(dna, k, progress_bar=False):
    """
    Among all possible k-mers, returns the one with minimum distance from a given collection of DNA segments. If multiple k-mers satisfy the requirement, choose any of them.
    :param dna: The collection of DNA segments (each segment is a string).
    :param k: The length of the k-mer to be found.
    :return: A pair with the median string and its score, a string and an integer respectively.
    """
    best_kmers, min_dist = all_median_strings(dna, k, progress_bar=False)
    return best_kmers[0], min_dist


def all_median_strings(dna, k, progress_bar=False):
    """
    Among all possible k-mers, returns those with minimum distance from a given collection of DNA segments.
    :param dna: The collection of DNA segments (each segment is a string).
    :param k: The length of the k-mers to be found.
    :return: A pair with median strings and the minimum score, respectively a list of strings and an integer.
    """
    min_dist = float('inf')
    best_kmers = None
    kmers_count = 4 ** k
    for kmer in tqdm(iterable=all_kmers(k), total=kmers_count, miniters=64) if progress_bar else all_kmers(k):
        dist = kmer_to_dna_distance(kmer, dna)
        if dist < min_dist:
            min_dist = dist
            best_kmers = [kmer]
        elif dist == min_dist:
            best_kmers.append(kmer)
    return best_kmers, min_dist


"""
The brute force approach finds the median string, but it takes time exponential in the k-mer length, k. Also, we don't know how long the motif is, therefore we need to find the median string for different values of k, and then try and deduce the motif length. Brute force becomes unsuitable as k increases.

We now look into more efficient possibilities, first a greedy algorithm, and then a randomized one. Instead of looking for a median string, we will be looking for motifs in every DNA segment, such that they are as similar to each other as possible. We will do so while enumerating the possible choices of motifs in a "clever" way, as opposed to enumerating all possible choices, which would be computationally intractable.

For starters, we define a profile for k-mers. Given a set of k-mers, the profile tells with what frequency each nucleotide appears in each of the k positions of the k-mers. In the implementation there is an extra parameter, pseudocount, that we will need later, for now it is set to a default of 0.
"""


def motifs_profile(motif, pseudocount=0):
    """
    Determine the frequency profile matrix for a motif using Laplace's Rule of Succession.
    :param motif: The given motif.
    :param pseudocount: The constant value added to the count of every nucleotide in every position, before normalising the counts into frequencies.
    :return: The profile, a dictionary with keys 'A', 'C', 'G' and 'T' which associates every nucleotide with a list of frequencies of length k. Each frequency is a real number between 0 and 1 included, and the sum of frequencies for a given position across all nucleotides is 1.
    """
    n_rows = len(motif)
    n_cols = len(motif[0])
    the_profile = {'A': [pseudocount] * n_cols, 'C': [pseudocount] * n_cols, 'G': [pseudocount] * n_cols,
                   'T': [pseudocount] * n_cols}
    for column in range(0, n_cols):
        for row in range(0, n_rows):
            nucleotide = motif[row][column]
            the_profile[nucleotide][column] += 1

    col_total = pseudocount * 4 + 1 * n_rows
    for col in range(0, n_cols):
        the_profile['A'][col] /= col_total
        the_profile['C'][col] /= col_total
        the_profile['G'][col] /= col_total
        the_profile['T'][col] /= col_total

    return the_profile


'''
Next thing, given a DNA segment, we want to find, among all its k-mers, the most likely based on a motifs.   
'''


def profile_most_probable_kmer(text, k, profile):
    """
    Choose the k-mer within a DNA segment that is most probable given a profile.
    :param text: The DNA segment, a string.
    :param k: The size of the wanted k-mer.
    :param profile: The frequency profile; a dictionary with keys 'A', 'C', 'G' and 'T' which associates every nucleotide with a list of probabilities of length k.
    :return: The most probable k-mer, a string.
    """

    most_prob_kmer = None
    highest_prob = -1
    for kmer in kmers_from_dna(text, k):
        prob = functools.reduce(operator.mul,
                                [profile[kmer[i]][i] for i in range(0, k)], 1)
        if prob > highest_prob:
            highest_prob = prob
            most_prob_kmer = kmer
    return most_prob_kmer


"""
Finally, given a choice of motifs of the same length, we define a way to score it, based on how similar they are to each other. The more similar they are, the lower the score, with 0 the score for motifs that are all identical.
"""


def score_motif(motifs):
    """
    Determine the score for motifs, i.e. a set of k-mers all of the same length. The more dissimilar the k-mers are in every position (from 0 to k-1), the higher the score. A score of 0 indicates that all k-mers in the motif are identical.
    :param motifs: The motif to be scored.
    :return: The score, an integer number.
    """
    n_rows = len(motifs)
    n_cols = len(motifs[0])
    scores = []
    for col in range(0, n_cols):
        counts = Counter([motifs[row][col] for row in range(0, n_rows)])
        scores.append(n_rows - max(counts.values()))
    total_score = sum(scores)
    return total_score


def relative_entropy(motifs, nucleotides_freq):
    profile = motifs_profile(motifs, pseudocount=0)
    k = len(motifs[0])
    total = .0
    for j in range(0, k):
        for r in 'ACGT':
            the_log2 = log2(profile[r][j]) if profile[r][j] > 0 else .0
            total += profile[r][j] * the_log2 - profile[r][j] * log2(nucleotides_freq[r])
    return -total


"""
We are now all set to implement a greedy motifs search. Again, there is an additional parameter pseudocount, which we will use later
"""


def greedy_motifs_search(dna, k, pseudocount=0):
    """
    Find the motifs that best fit a collection of DNA segments (have the lowest score), using a greedy search. If multiple motifs have the same score, then provide just one of them.
    :param dna: A sequence of DNA segments (strings), not necessarily of the same length.
    :param k: The size of motifs to be discovered; the same for every motif, an integer number.
    :param pseudocount: The constant value added to the count of every nucleotide in every position, for computation of the frequency profile.
    :return: The discovered motifs, a list of strings.
    """
    t = len(dna)
    # Keep track of the best set of motifs so far
    best_score = float('inf')
    best_motifs = None
    # For every k-mer in the first DNA segment...
    for kmer in kmers_from_dna(dna[0], k):
        # ... start with list of motifs that only contains that k-mer...
        motifs = [kmer]
        # ... then add to it one k-mer per DNA segment, from the second segment on
        for i in range(1, t):
            ''' The k-mer to be added to the list of motifs is choosen based on which
            is the most probable, giving the profile of the k-mers already in the list '''
            the_profile = motifs_profile(motifs, pseudocount=pseudocount)
            most_prob_kmer = profile_most_probable_kmer(dna[i], k, the_profile)
            motifs.append(most_prob_kmer)
        # Once the list of motifs contain one k-mer per DNA segment, score it
        score = score_motif(motifs)
        # Update the best set of motifs and their score (lower score is better)
        if score < best_score:
            best_motifs = motifs
            best_score = score
    return best_motifs


'''
The greedy motif search performs poorly. Issue is, it constructs and update a profile which is rather sparse, with many values set to 0. As a consequence, the likelyhood assigned to most k-mers in every DNA segment is zero, and most k-mers are not evaluated at all for inclusion among the motifs. K-mers that are good candidate to be motifs are instead skipped altogether. 

To address the issue, we adopt "Laplace's Rule fo Succession": when we should set a value of the profile to 0, we set it to some small value instead, like if that nucleotide was observed a few times in that position (typically one). This is why we have the pseudocount paramter in functions motifs_profile() and greedy_motif_search(): it is a positive integer value to be used instead of 0 in counting nucleotides.   

Let's now explore an alternative kind of algorithms, randomised algorithms. 

We start with a random choice of motifs, taken one from every DNA segment. We build a profile for the current choice of motifs, then choose a new set of motifs (one per DNA segment) based on the profile; then we iterate, as long as we get a choice of motifs with better (lower) score than before. 
'''


def randomized_motif_search(dna, k, seed=None):
    """
    Find the motifs that best fit a collection of DNA segments (have the lowest score), using a randomised search. Setting the same seed with random.seed() before calling the function, will produce the same result every time.
    :param dna: The collection of DNA segments (strings), not necessarily of the same length.
    :param k: The size of motifs to be discovered; the same for every motif, an integer number.
    :param seed: seed for random number generator initialisation: if set to None, no initialisation is performed.
    :return: The discovered motifs, a list of strings.
    """
    if seed is not None:
        random.seed(seed)
    t = len(dna)
    s = len(dna[0])
    assert s >= k
    motif = [dna[row][i:i + k] for row, i in zip(range(0, t), [randint(0, s - k - 1) for _ in range(0, t)])]
    best_motif = motif
    best_motif_score = score_motif(best_motif)
    while True:
        the_profile = motifs_profile(motif, pseudocount=1)
        motif = [profile_most_probable_kmer(text, k, the_profile) for text in dna]
        score = score_motif(motif)
        if score < best_motif_score:
            best_motif = motif
            best_motif_score = score
        else:
            return best_motif, best_motif_score


'''
We then run the randomised search repeatedly, and keep the motifs with the best (lowest) score.
'''


def mc_randomized_motif_search(dna, k, times, seed=None):
    """
    Repeat the randomised motif search in DNA segments a given number of times, and select the best result. Setting the same seed with random.seed() before calling the function, will produce the same result every time.
    :param dna: The collection of DNA segments (strings), not necessarily of the same length.
    :param k: The size of motifs to be discovered; the same for every motif, an integer number.
    :param times: The number of time the motifs search has to be repeated
    :param seed: seed for random number generator initialisation: if set to None, no initialisation is performed.
    :return: The discovered motifs, a list of strings.
    """
    if seed is not None:
        random.seed(seed)
    best_score = float('inf')
    best_motif = None
    for _ in range(0, times):
        motif, score = randomized_motif_search(dna, k, None)
        if score < best_score:
            best_score = score
            best_motif = motif
    return best_motif


def nucleotides_frequency(dna):
    nucleotides_count_by_segm = [Counter(dna[i]) for i in range(0, len(dna))]
    nucleotides_freq = functools.reduce(operator.add, nucleotides_count_by_segm, Counter(A=0, C=0, G=0, T=0))
    total_nucleotides_no = nucleotides_freq['A'] + nucleotides_freq['C'] + nucleotides_freq['G'] + nucleotides_freq['T']
    assert total_nucleotides_no == sum([len(segment) for segment in dna])
    nucleotides_freq['A'] /= total_nucleotides_no
    nucleotides_freq['C'] /= total_nucleotides_no
    nucleotides_freq['G'] /= total_nucleotides_no
    nucleotides_freq['T'] /= total_nucleotides_no
    return nucleotides_freq


'''
The randomised search samples a new set of motifs at every iteration. A more prudent approach is to sample only one new motif per iteration, and replace it in the current list of motifs. 

Algorithm "Gibbs Sampler", at every iteration:
 - discards one motif, randomly chosen, say corresponding to DNA segment i;
 - calculates a profile, based on the remaining list of motifs;
 - based on the profile, chooses the most likely motif in the i-th DNA segment, and insert it into the motifs list;
 - keeps track of the best (lowest score) list of motifs.  
'''


def gibbs_sampler(dna, k, n, seed=None):
    """
    Find the motifs that best fit a collection of DNA segments (have the lowest score), using a Gibbs sampler. Setting the same seed with random.seed() before calling the function, will produce the same result every time.
    :param dna: The collection of DNA segments (strings), not necessarily of the same length.
    :param k: The size of motifs to be discovered; the same for every motif, an integer number.
    :param n: The number of iterations to be performed, an integer number.
    :param seed: seed for random number generator initialisation: if set to None, no initialisation is performed.
    :return: The discovered motifs, a list of strings.
    """
    if seed is not None:
        random.seed(seed)
    t = len(dna)  # No. of strings in dna
    s = len(dna[0])  # Length of every string in dna
    assert s >= k

    nucleotides_freq = nucleotides_frequency(dna)

    # Randomly select one motif per dna string, each of length  k, as a starting set of motifs
    motifs = [dna[row][i:i + k] for row, i in zip(range(0, t), [randint(0, s - k - 1) for _ in range(0, t)])]
    # Initialise the needful to keep track of the best motifs discovered so far
    best_motifs = deepcopy(motifs)  # Deep copy needed to prevent updates to motifs[i] from changing best_motifs[i]
    best_motif_score = score_motif(best_motifs)
    scores = []
    for _ in range(0, n):  # Basic stopping criteria
        # Randomly choose one of the motifs
        i = randint(0, t - 1)
        # Find the probability motifs_profile of the motifs, without the randomly choosen one
        motif_ex_i = motifs[:i] + motifs[i + 1:]
        the_profile = motifs_profile(motif_ex_i, pseudocount=1)
        # Determine the likelihood of every k-mer in the i-th string of dna, based on the just calculated profile
        proportions = []
        for kmer in kmers_from_dna(dna[i], k):
            prop = functools.reduce(operator.mul,
                                    [the_profile[nucleotide][col] for nucleotide, col in
                                     zip(kmer, range(0, len(kmer)))], 1)
            proportions.append(prop)
        sum_proportions = sum(proportions)
        prob_distr = [prop / sum_proportions for prop in proportions]
        assert len(prob_distr) == s - k + 1
        # Sample one kmer from the i-th string in dna, based on the just calculated probability distribution
        drafted_kmer_i = choices(range(0, len(prob_distr)), weights=prob_distr)[0]
        # Replace the i-th motif in the current set of motifs with the sampled one
        motifs[i] = dna[i][drafted_kmer_i:drafted_kmer_i + k]
        # Check if the obtained motif is the best so far
        # score = score_motif(motifs)
        score = relative_entropy(motifs, nucleotides_freq)
        scores.append(score)
        if score < best_motif_score:
            # Deep copy needed, or else code above that updates motifs[i] will overwrite best_motifs[i]
            best_motifs = deepcopy(motifs)
            best_motif_score = score
    return best_motifs, scores


def consensus_from_motifs(motifs):
    n_cols = len(motifs[0])
    n_rows = len(motifs)
    consensus = []
    for col in range(0, n_cols):
        counts = Counter([motifs[row][col] for row in range(0, n_rows)])
        consensus_nucleotide, _ = counts.most_common(1)[0]
        consensus.append(consensus_nucleotide)
    consensus_as_str = ''.join(consensus)
    return consensus_as_str


def sample_random_relative_entropy(dna, k, n, seed=None):
    if seed is not None:
        random.seed(seed)
    nucleotides_freq = nucleotides_frequency(dna)
    n_rows = len(dna)
    n_cols = len(dna[0])
    samples = []
    for _ in range(0, n):
        motifs = [dna[segment][i:i + k] for (segment, i) in
                  zip(range(0, n_rows), [randint(0, n_cols - k - 1) for _ in range(0, n_rows)])]
        rel_entropy = relative_entropy(motifs, nucleotides_freq)
        samples.append(rel_entropy)
    return samples


def generate_dataset(file_name, n_segments, segment_length, motif_length, proportions=None, seed=None):
    if seed is not None:
        random.seed(seed)

    if proportions is not None:
        weights = (proportions['A'], proportions['C'], proportions['G'], proportions['T'])
        assert isclose(sum(weights), 1)
    else:
        weights = [.25] * 4

    motif = ''.join(choices(population='ACGT', weights=weights, k=motif_length))
    with open(file_name, 'w') as fasta_file:
        for i_segment in range(0, n_segments):
            fasta_file.write('>synth' + str(i_segment) + '\n')
            motif_pos = random.randint(0, segment_length - motif_length - 1)
            prefix = ''.join(choices(population='ACGT', weights=weights, k=motif_pos))
            postfix = ''.join(
                choices(population='ACGT', weights=weights, k=segment_length - motif_length - len(prefix)))
            segment = prefix + motif + postfix + '\n'
            fasta_file.write(segment)

    return motif


def main3():
    motif = generate_dataset(file_name='synthetic.fasta',
                             n_segments=36,
                             segment_length=250,
                             motif_length=20,
                             proportions=None,
                             seed=42)
    print(motif)


def main2():  # TODO clean up this
    backend = matplotlib.get_backend()
    interactive = matplotlib.is_interactive()
    print('Using backend', backend, ', interactive = ', interactive)

    parsed = SeqIO.parse('upstream250.txt', 'fasta')
    # parsed = SeqIO.parse('synthetic.fasta', 'fasta')
    segment_IDs = []
    segments = []
    for record in parsed:
        segment_IDs.append(record.id)
        segments.append(str(record.seq))

    """
    random.seed(42)
    for k in range(8, 42, 2):
        samples = sample_random_relative_entropy(dna=segments, k=k, n=10000, seed=None)
        mean = statistics.mean(samples)
        var = statistics.variance(samples)
        print(k, mean/k, var)
    """

    """
    best_scores = []
    for k in range(8, 42, 2):
        motifs, scores = gibbs_sampler(segments, k=k, n=20000, seed=42)
        print('\n k=', k)
        print(motifs)
        consensus = consensus_from_motifs(motifs)
        print('Consensus is', consensus)
        print('Best score is', min(scores), 'at iteration#', argmin(scores))
        best_scores.append(min(scores) / k)
        # plt.plot(scores)
        # plt.show(block=True)
    """
    k = 20
    motifs, scores = gibbs_sampler(segments, k=k, n=20000, seed=42)
    print('\n k=', k)
    print(motifs)
    consensus = consensus_from_motifs(motifs)
    print('Consensus is', consensus)
    print('Best score is', min(scores), 'at iteration#', argmin(scores))
    plt.plot(scores)
    plt.show(block=True)

    # plt.plot(range(8, 42, 2), best_scores)
    # plt.show(block=True)

    # median = median_string(segments, k=5)
    # print('Median string is', median)


if __name__ == '__main__':
    main2()


def shuffle_motifs(motifs, seed=None):
    """
    Given a collection of motifs, returns a collection of randomly generated, different motifs, with the same score.
    :param motifs: The given motifs, a sequence of strings.
    :param seed: The seed for initialization of the random number generator, an integer; if set to None, no initialization is performed.
    :return: The randomly generated motifs, with the same score as the given ones; a list of strings.
    """
    if seed is not None:
        random.seed = 42
    t = len(motifs)  # Number of motifs
    k = len(motifs[0])  # Length of k-mers
    motifs_T = list(map(list, zip(*motifs)))  # Transpose motifs
    shuffle(motifs_T)  # Shuffle
    motifs_T_rot = [row[pos:] + row[:pos] for row, pos in zip(motifs_T, choices(list(range(0, t)), k=k))]  # Rotate
    motifs = list(map(list, zip(*motifs_T_rot)))  # Transpose back
    return motifs
