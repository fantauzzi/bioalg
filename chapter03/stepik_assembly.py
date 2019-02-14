import assembly


def fetch_int():
    return int(input())


def fetch_int_pair():
    line = input()
    k, d = line.split(sep=' ')
    return int(k), int(d)


def fetch_lines():
    lines = []
    try:
        while True:
            item = input()
            lines.append(item)
    except EOFError:
        pass
    return lines


def main_reconstruct_string():
    k = int(input())
    text = []
    try:
        while True:
            item = input()
            text.append(item)
    except EOFError:
        pass
    assert k == len(text[0])

    res = assembly.reconstruct_string_from_kmers(text)
    print(res)


def main_eulerian_cycle():
    """
    k = int(input())
    text = input()
    adj = de_brujin_graph(k, text)
    print_graph(adj)
    """

    text = []
    try:
        while True:
            item = input()
            text.append(item)
    except EOFError:
        pass
    adj = assembly.parse_graph(text)
    cycle = assembly.eulerian_cycle(adj)
    assembly.print_cycle(cycle)


def main_eulerian_path():
    """
    k = int(input())
    text = input()
    adj = de_brujin_graph(k, text)
    print_graph(adj)
    """

    text = []
    try:
        while True:
            item = input()
            text.append(item)
    except EOFError:
        pass
    adj = assembly.parse_graph(text)
    path = assembly.eulerian_path(adj)
    assembly.print_cycle(path)


def main_k_universal():
    k = int(input())
    res = assembly.make_k_universal_string(k)
    print(res)


def main_paired_reads():
    line = input()
    k, d = line.split(sep=' ')
    k, d = int(k), int(d)

    text = []
    try:
        while True:
            item = input()
            text.append(item)
    except EOFError:
        pass

    gapped_kmers = [item.split(sep='|') for item in text]
    assert k == len(gapped_kmers[0][0])
    res = assembly.reconstruct_string_from_paired_reads(d, gapped_kmers)
    print(res)


def main_max_non_br_paths():
    text = []
    try:
        while True:
            item = input()
            text.append(item)
    except EOFError:
        pass
    adj = assembly.parse_graph(text)
    paths = assembly.max_no_branch_paths(adj)
    for path in paths:
        assembly.print_cycle(path)
        print()


def main_contigs_from_kmers():
    kmers = []
    try:
        while True:
            item = input()
            kmers.append(item)
    except EOFError:
        pass
    contigs = assembly.contigs_from_kmers(kmers)
    print(*contigs, sep=' ')


def main_path_to_genome():
    kmers = fetch_lines()
    genome = assembly.path_to_genome(kmers)
    print(genome)


def main_overlap_graph():
    kmers = fetch_lines()
    adj, count = assembly.overlap_graph(kmers)
    assembly.print_overlap_graph(adj, count)


def main_de_brujin_from_kmers():
    kmers = fetch_lines()
    adj = assembly.de_brujin_graph_from_kmers(kmers)
    assembly.print_graph(adj)

if __name__ == '__main__':
    pass

