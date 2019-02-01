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


# Input:  Strings Pattern and Text along with an integer d
# Output: A list containing all starting positions where Pattern appears
# as a substring of Text with at most d mismatches
def ApproximatePatternMatching(Text, Pattern, d):
    if len(Pattern) > len(Text):
        return []
    positions = []  # initializing list of positions
    # your code here
    sub_text = Text[:len(Pattern)]
    for i in range(0, len(Text) - len(Pattern) + 1):
        if i > 0:
            sub_text = sub_text[1:] + Text[len(Pattern) + i - 1]
        if HammingDistance(sub_text, Pattern) <= d:
            positions.append(i)
    return positions


# Input:  Strings Pattern and Text, and an integer d
# Output: The number of times Pattern appears in Text with at most d mismatches
def ApproximatePatternCount(Text, Pattern, d):
    positions = ApproximatePatternMatching(Text, Pattern, d)
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


# Write your FrequentWordsWithMismatches() function here, along with any subroutines you need.
# Your function should return a list.
def FrequentWordsWithMismatches(Text, k, d):
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


def DNA_complement(dna):
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
        all_neighbors = all_neighbors.union(neighbors(DNA_complement(current_kmer), d))

    all_neighbors = dict.fromkeys(all_neighbors, 0)
    max_so_far = float('-inf')
    for i in range(0, len(Text) - k + 1):
        current_kmer = Text[i:i + k]
        for kmer, count in all_neighbors.items():
            updated_count = count
            if HammingDistance(kmer, current_kmer) <= d:
                updated_count += 1
            kmer_c = DNA_complement(kmer)
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


def main():
    with open('skew_dataset.txt') as genome_file:
        genome = genome_file.read()

    running_skew = [0] + compute_skew(genome)
    mins = argsmin(running_skew)
    print("Minimum of skew at position(s)",mins)
    plot_skew(running_skew)


if __name__ == '__main__':
    main()
