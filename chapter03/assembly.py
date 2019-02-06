def composition(k, text):
    """
    Produces the list of k-mers in a string, in lexicographic order.
    :param k: The k-mer size (number of nucleotides).
    :param text: The string.
    :return: A list of strings, the k-mers in lexicographic order.
    """

    res = []
    for i in range(0, len(text) - k + 1):
        res.append(text[i:i + k])
    return sorted(res)


def path_to_genome(path):
    """
    Reconstruct a string from its genome path.
    :param path:The genome path, a sequence of k-mers.
    :return:The string corresponding to the given genome path.
    """

    nucleotides = [item[-1] for item in path[1:]]
    genome = path[0]+''.join(nucleotides)
    return genome


def main():
    path = []
    try:
        while True:
            item = input()
            path.append(item)
    except EOFError:
        pass
    res = path_to_genome(path)
    print(res)


if __name__ == '__main__':
    main()
