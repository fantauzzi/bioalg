import assembly
from pathlib import Path


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
    expected = {'ATGCG': [],
                'GCATG': ['CATGC'],
                'CATGC': ['ATGCG'],
                'AGGCA': ['GGCAT'],
                'GGCAT': ['GCATG']}
    assert adj == expected

    kmers2 = ['AAT', 'ATG', 'ATG', 'ATG', 'CAT', 'CCA', 'GAT', 'GCC', 'GGA', 'GGG', 'GTT', 'TAA', 'TGC', 'TGG', 'TGT']
    adj, count = assembly.overlap_graph(kmers2)
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
    assert cycle == ['0', '1', '2', '0']
    res = assembly.is_eulerian_cycle(adj, cycle)
    assert res

    adj = {'0': ['0']}
    cycle = assembly.eulerian_cycle(adj)
    res = assembly.is_eulerian_cycle(adj, cycle)
    assert cycle == ['0', '0']
    assert res

    adj = {'0': ['1'],
           '1': ['0']}
    cycle = assembly.eulerian_cycle(adj)
    res = assembly.is_eulerian_cycle(adj, cycle)
    assert res

    adj = {'0': ['1'],
           '1': ['1', '2', '3'],
           '2': ['3'],
           '3': ['0', '1']}
    cycle = assembly.eulerian_cycle(adj)
    assert cycle == ['0', '1', '3', '1', '1', '2', '3', '0']
    res = assembly.is_eulerian_cycle(adj, cycle)
    assert res

    adj = {'0': ['1', '2', '3', '4'],
           '1': ['0', '2'],
           '2': ['0', '3'],
           '3': ['0', '4'],
           '4': ['0', '1']}
    cycle = assembly.eulerian_cycle(adj)
    assert cycle == ['0', '4', '1', '2', '3', '4', '0', '3', '0', '2', '0', '1', '0']
    res = assembly.is_eulerian_cycle(adj, cycle)
    assert res

    with open(Path('test/dataset01.txt')) as dataset_file:
        dataset = dataset_file.readlines()
    adj = assembly.parse_graph(dataset)
    cycle = assembly.eulerian_cycle(adj)
    assert len(cycle) == 4001
    res = assembly.is_eulerian_cycle(adj, cycle)
    assert res


def test_euleria_path():
    adj = {'0': ['2'],
           '1': ['3'],
           '2': ['1'],
           '3': ['0', '4'],
           '6': ['3', '7'],
           '7': ['8'],
           '8': ['9'],
           '9': ['6']}

    path = assembly.eulerian_path(adj)
    assert path == ['6', '7', '8', '9', '6', '3', '0', '2', '1', '3', '4']

    adj = {'0': ['1']}
    path = assembly.eulerian_path(adj)
    assert path == ['0', '1']

    for i in range(1, 6):
        with open(Path('test/test' + str(i) + '.txt')) as input_file:
            test_case_text = input_file.readlines()
        adj = assembly.parse_graph(test_case_text)
        path = assembly.eulerian_path(adj)
        with open(Path('test/test' + str(i) + 'out.txt')) as input_file:
            expected_txt = input_file.readlines()[0].rstrip()
        expected = expected_txt.split('->')
        assert expected == path, 'i={}'.format(i)
