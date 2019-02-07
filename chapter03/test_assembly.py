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
    assembly.print_overlap_graph(adj, count)
    assert True

