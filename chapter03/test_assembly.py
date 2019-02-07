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
    print()
    assembly.print_overlap_graph(adj, count)
    assert True

    kmers2 = ['AAT', 'ATG', 'ATG', 'ATG', 'CAT', 'CCA', 'GAT', 'GCC', 'GGA', 'GGG', 'GTT', 'TAA', 'TGC', 'TGG', 'TGT']
    adj, count = assembly.overlap_graph(kmers2)
    print()
    assembly.print_overlap_graph(adj, count)
    assert True
