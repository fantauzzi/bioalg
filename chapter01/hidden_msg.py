import matplotlib.pyplot as plt
import matplotlib

plt.ion()
backend = matplotlib.get_backend()
interactive = matplotlib.is_interactive()
print('Using backend', backend, ', interactive = ', interactive)


def argsmin(list_like):
    """
    Given a sequence, returns the position in the sequence of its minimum item. If multiple items in the list have minimum value, the position of the first of them is returned.
    :param list_like: The sequence.
    :return: The 0-indexed position, an integer.
    """
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
    """
    Given a sequence, returns the position in the sequence of its maximum item. If multiple items in the list have maximum value, the position of the first of them is returned.
    :param list_like: The sequence.
    :return: The 0-indexed position, an integer.
    """
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
    """
    Returns the skew of a genome at every position in a genome. The skew at a given position is the difference between the counts of G and C from the beginning of the genome up to the given position included.
    :param genome: The genome, a string.
    :return: The skew, a list of integer numbers, as long as the genome.
    """
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


def hamming_distance(p, q):
    """
    Returns the Hamming distance betweem two strings of the same length.
    :param p: The first string.
    :param q: The second string.
    :return: The Hammind distance, an integer.
    """
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
        if hamming_distance(sub_text, pattern) <= d:
            positions.append(i)
    return positions


def approximate_pattern_count(text, pattern, d):
    """
    Returns the number of different places a given string appears as a substring of a text with at most a given number of mismatches.
    :param text: The text, a string.
    :param pattern: The given string.
    :param d: The maximum allowed number of mismatches at every matching location, an integer.
    :return: The count of approximate matches, an integer.
    """
    positions = approximate_pattern_matching(text, pattern, d)
    return len(positions)


def immediate_neighbors(kmer):
    """
    Returns the 1-neighbors of a given k-mer.
    :param kmer: The k-mer, a DNA string.
    :return: All the 1-neighbors of the k-mer, a list of strings.
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
    Returns the d-neighbors of a given k-mer.
    :param kmer: The k-mer, a DNA string.
    :param d: The value for d paramter, a non-negative integer.
    :return: All the d-neighbots of the k-mer, a list of strings.
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


def frequent_words_with_mismatches(text, k, d):
    """
    Returns the most frequent k-mers in a text with at most a given number of mismatches.
    :param text: The text, a string.
    :param k: The length of every k-mer, an integer.
    :param d: The maximum number of mismatches allowed in each k-mer, an integer.
    :return: The most frequent k-mer appearing in text with at most d mismatches, a list of strings.
    """
    all_neighbors = set()
    for i in range(0, len(text) - k + 1):
        current_kmer = text[i:i + k]
        all_neighbors = all_neighbors.union(neighbors(current_kmer, d))

    all_neighbors = dict.fromkeys(all_neighbors, 0)
    max_so_far = float('-inf')
    for i in range(0, len(text) - k + 1):
        current_kmer = text[i:i + k]
        for kmer, count in all_neighbors.items():
            if hamming_distance(kmer, current_kmer) <= d:
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


def frequent_words_with_complements(text, k, d):
    """
    Returns the k-mers that, along with their reverse complements, appear the most in a text, with at most a given number of mismatches.
    :param text: The text, a string.
    :param k: The length of every k-mer, an integer.
    :param d: The maximum number of mismatches allowed each k-mer, an integer.
    :return: The most frequent k-mers, with reverse complements and at most d mismatches, in the text, a list of strings.
    """
    all_neighbors = set()
    for i in range(0, len(text) - k + 1):
        current_kmer = text[i:i + k]
        all_neighbors = all_neighbors.union(neighbors(current_kmer, d))
        all_neighbors = all_neighbors.union(neighbors(DNA_rev_complement(current_kmer), d))

    all_neighbors = dict.fromkeys(all_neighbors, 0)
    max_so_far = float('-inf')
    for i in range(0, len(text) - k + 1):
        current_kmer = text[i:i + k]
        for kmer, count in all_neighbors.items():
            updated_count = count
            if hamming_distance(kmer, current_kmer) <= d:
                updated_count += 1
            kmer_c = DNA_rev_complement(kmer)
            if hamming_distance(kmer_c, current_kmer) <= d:
                updated_count += 1
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
    """
    Plots the chart of a running skew with matplotlib.
    :param running_skew: The running skew, a sequence of integer numbers, as returned by compute_skew()
    """
    plt.plot(running_skew)
    plt.xlabel('position')
    plt.ylabel('running count')
    plt.title('G-C Skew')
    plt.grid(True)
    plt.show(block=True)


def main():
    """
    Plots a sample C-G skew chart.
    """
    with open('skew_dataset.txt') as genome_file:
        genome = genome_file.read()

    running_skew = [0] + compute_skew(genome)
    mins = argsmin(running_skew)
    print("Minimum of skew at position(s)", mins)
    plot_skew(running_skew)


def compute_frequencies(text, k):
    """
    For every k-mer that appears in a text, returns the k-mer and the number of times it appears.
    :param text: The text, a string.
    :param k: The size of the k-mers, an integer.
    :return: A dictionary, associating every k-mer that appears at least once in the text with the count of its appearances (a string to integer association).
    """
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
