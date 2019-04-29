import matplotlib.pyplot as plt
import matplotlib

plt.ion()
backend = matplotlib.get_backend()
interactive = matplotlib.is_interactive()
print('Using backend', backend, ', interactive = ', interactive)


def argsmin(list_like):
    min_so_far = float('inf')
    indices = []  # Unnecessary, here just for clarity
    for i, item in enumerate(list_like):
        if item == min_so_far:
            indices.append(i)
        elif item < min_so_far:
            min_so_far = item
            indices = [i]
    return indices


def argsmax(list_like):
    max_so_far = float('-inf')
    indices = []  # Unnecessary, here just for clarity
    for i, item in enumerate(list_like):
        if item == max_so_far:
            indices.append(i)
        elif item > max_so_far:
            max_so_far = item
            indices = [i]
    return indices


def compute_skew(genome):
    genome = str.lower(genome)
    skew = []
    gc_diff_count = 0
    for neuclotide in genome:
        if neuclotide == 'g':
            gc_diff_count += 1
        elif neuclotide == 'c':
            gc_diff_count -= 1
        skew.append(gc_diff_count)
    return skew


def min_skew(genome):
    skew = compute_skew(genome)
    res = argsmin(skew)
    res = [item + 1 for item in res]  # Count positions from 1 instead of 0, as requested by specs
    return res


def MinimumSkew(genome):
    return min_skew(genome)


# Input:  Two strings p and q
# Output: An integer value representing the Hamming Distance between p and q.
def HammingDistance(p, q):
    assert len(p) == len(q)
    dist = 0
    for char1, char2 in zip(p, q):
        if char1 != char2:
            dist += 1
    return dist


def approximate_pattern_matching(text, pattern, d):
    """
    Returns all starting positions where a given pattern appears as a substring of a text, with a given number of mismatches at most.
    :param text: The text, a string.
    :param pattern: The given pattern, a string.
    :param d: The maximum number of mismatches for every position in the text, an integer.
    :return: The list of amtching positions in text, zero-indexed, a list of integers in increasing order.
    """
    if len(pattern) > len(text):
        return []
    positions = []  # initializing list of positions
    sub_text = text[:len(pattern)]
    for i in range(0, len(text) - len(pattern) + 1):
        if i > 0:
            sub_text = sub_text[1:] + text[len(pattern) + i - 1]
        if HammingDistance(sub_text, pattern) <= d:
            positions.append(i)
    return positions


# Input:  Strings Pattern and Text, and an integer d
# Output: The number of times Pattern appears in Text with at most d mismatches
def approximate_pattern_count(Text, Pattern, d):
    """

    :param Text:
    :param Pattern:
    :param d:
    :return:
    """
    positions = approximate_pattern_matching(Text, Pattern, d)
    return len(positions)


def immediate_neighbors(kmer):
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


def frequent_words_with_mismatches(Text, k, d):
    """
    Returns the most frequent k-mers in a text with at most a given number of mismatches.
    :param Text: The text, a string.
    :param k: The length of every k-mer, an integer.
    :param d: The maximum number of mismatches allowed in each k-mer, an integer.
    :return: The most frequent k-mer appearing in text with at most d mismatches, a list of strings.
    """
    all_neighbors = set()
    for i in range(0, len(Text) - k + 1):
        current_kmer = Text[i:i + k]
        all_neighbors = all_neighbors.union(neighbors(current_kmer, d))

    all_neighbors = dict.fromkeys(all_neighbors, 0)
    max_so_far = float('-inf')
    for i in range(0, len(Text) - k + 1):
        current_kmer = Text[i:i + k]
        for kmer, count in all_neighbors.items():
            if HammingDistance(kmer, current_kmer) <= d:
                all_neighbors[kmer] = count + 1
                if count + 1 > max_so_far:
                    max_so_far = count + 1

    res = []
    for k, v in all_neighbors.items():
        if v == max_so_far:
            res.append(k)
    return res


def DNA_rev_complement(dna):
    """
    Returns the reverse complement of a DNA string.
    :param dna: The DNA string, a string composed of symbols from 'ACGT' (case sensitive).
    :return: The reverse complement, a string.
    """
    complements = {'A': 'T',
                   'T': 'A',
                   'G': 'C',
                   'C': 'G'}
    complement = [complements[nucleotide] for nucleotide in reversed(dna)]
    as_string = ''.join(complement)
    return as_string


# Place your FrequentWordsWithMismatchesAndReverseComplements() function here, along with any subroutines you need.
# Your function should return a list.
def FrequentWordsWithMismatchesAndReverseComplements(Text, k, d):
    all_neighbors = set()
    for i in range(0, len(Text) - k + 1):
        current_kmer = Text[i:i + k]
        all_neighbors = all_neighbors.union(neighbors(current_kmer, d))
        all_neighbors = all_neighbors.union(neighbors(DNA_rev_complement(current_kmer), d))

    all_neighbors = dict.fromkeys(all_neighbors, 0)
    max_so_far = float('-inf')
    for i in range(0, len(Text) - k + 1):
        current_kmer = Text[i:i + k]
        for kmer, count in all_neighbors.items():
            updated_count = count
            if HammingDistance(kmer, current_kmer) <= d:
                updated_count += 1
            kmer_c = DNA_rev_complement(kmer)
            if HammingDistance(kmer_c, current_kmer) <= d:
                updated_count += 1
            # assert kmer != kmer_c
            if updated_count > count:
                all_neighbors[kmer] = updated_count
                if updated_count > max_so_far:
                    max_so_far = updated_count

    res = []
    for k, v in all_neighbors.items():
        if v == max_so_far:
            res.append(k)
    return res


def plot_skew(running_skew):
    plt.plot(running_skew)
    plt.xlabel('position')
    plt.ylabel('running count')
    plt.title('G-C Skew')
    plt.grid(True)
    plt.show(block=True)


def main():  # TODO remove this!
    with open('skew_dataset.txt') as genome_file:
        genome = genome_file.read()

    running_skew = [0] + compute_skew(genome)
    mins = argsmin(running_skew)
    print("Minimum of skew at position(s)", mins)
    plot_skew(running_skew)


def compute_frequencies(text, k):
    frequencies = {}
    for i in range(0, len(text) - k + 1):
        kmer = text[i:i + k]
        current = frequencies.get(kmer, 0)
        frequencies[kmer] = current + 1
    return frequencies


def find_clumps(genome, k, l, t):
    """
    Returns the k-mers that form an (l,t)-clump within a given genome. A k-mer forms an (l,t)-clump in a string genome if there is an interval of genome of length l in which the k-mer appears at least t times.
    :param genome: The genome, a string.
    :param k: The size of the k-mers, an integer.
    :param l: The length of the cluster, an integer.
    :param t: The minimum number of times a k-mer is required to appear in the interval, an integer.
    :return: The k-mers forming an (l,t)-clump, a list of string.
    """

    assert t >= 1
    clumped_kmers = set()
    for i in range(0, len(genome) - l):
        text = genome[i:i + l]
        frequencies = compute_frequencies(text, k)
        clumped_kmers |= {kmer for kmer in frequencies.keys() if frequencies[kmer] >= t}
    return list(clumped_kmers)
