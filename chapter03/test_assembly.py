import assembly


def test_composition():
    k = 5
    text = 'CAATCCAAC'
    res = assembly.composition(k, text)
    assert res == ['AATCC', 'ATCCA', 'CAATC', 'CCAAC', 'TCCAA']


def test_path_to_genome():
    path = ['ACCGA', 'CCGAA', 'CGAAG', 'GAAGC', 'AAGCT']
    genome = assembly.path_to_genome(path)
    assert genome == 'ACCGAAGCT'


def test_overlap_graph():
    kmers = ['ATGCG', 'GCATG', 'CATGC', 'AGGCA', 'GGCAT']
    adj, count = assembly.overlap_graph(kmers)
    # print()
    # assembly.print_overlap_graph(adj, count)
    expected = {'ATGCG': [],
                'GCATG': ['CATGC'],
                'CATGC': ['ATGCG'],
                'AGGCA': ['GGCAT'],
                'GGCAT': ['GCATG']}
    assert adj == expected

    kmers2 = ['AAT', 'ATG', 'ATG', 'ATG', 'CAT', 'CCA', 'GAT', 'GCC', 'GGA', 'GGG', 'GTT', 'TAA', 'TGC', 'TGG', 'TGT']
    adj, count = assembly.overlap_graph(kmers2)
    # print()
    # assembly.print_overlap_graph(adj, count)
    assert len(adj) == len(count) == 13


def test_de_brujin_graph():
    k = 4
    text = 'AAGATTCTCTAAGA'
    adj = assembly.de_brujin_graph(k, text)
    expected = {'AAG': ['AGA', 'AGA'],
                'AGA': ['GAT'],
                'ATT': ['TTC'],
                'CTA': ['TAA'],
                'CTC': ['TCT'],
                'GAT': ['ATT'],
                'TAA': ['AAG'],
                'TCT': ['CTC', 'CTA'],
                'TTC': ['TCT'],
                }
    assert adj == expected


def test_de_brujin_graph_from_kmers():
    kmers = ['GAGG',
             'CAGG',
             'GGGG',
             'GGGA',
             'CAGG',
             'AGGG',
             'GGAG']

    expected = {'AGG': ['GGG'],
                'CAG': ['AGG', 'AGG'],
                'GAG': ['AGG'],
                'GGA': ['GAG'],
                'GGG': ['GGG', 'GGA']}

    adj = assembly.de_brujin_graph_from_kmers(kmers)
    assert adj == expected


def test_eulerian_cycle():
    adj = {'0': ['1'],
           '1': ['2'],
           '2': ['0']}
    cycle = assembly.eulerian_cycle(adj)
    print()
    print(cycle)
    # print(assembly.cycle_from_adj(cycle, '1'))

    adj = {'0': ['0']}
    cycle = assembly.eulerian_cycle(adj)
    print()
    print(cycle)

    adj = {'0': ['1'],
           '1': ['0']}
    cycle = assembly.eulerian_cycle(adj)
    print()
    print(cycle)
    # print(assembly.cycle_from_adj(cycle, '1'))

    adj = {'0': ['1'],
           '1': ['1', '2', '3'],
           '2': ['3'],
           '3': ['0', '1']}
    cycle = assembly.eulerian_cycle(adj)
    print()
    print(cycle)
    # print(assembly.cycle_from_adj(cycle, '1'))

    adj = {'0': ['1', '2', '3', '4'],
           '1': ['0', '2'],
           '2': ['0', '3'],
           '3': ['0', '4'],
           '4': ['0', '1']}
    cycle = assembly.eulerian_cycle(adj)
    print()
    print(cycle)
    # print(assembly.cycle_from_adj(cycle, '1'))
