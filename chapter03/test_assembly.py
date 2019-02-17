import assembly
from pathlib import Path


def fetch_int_pair_from_file(the_file):
    line = the_file.readline()
    k, d = line.split(sep=' ')
    return int(k), int(d)


def fetch_lines_from_file(the_file):
    lines = the_file.readlines()
    stripped = [line.rstrip() for line in lines]
    return stripped


def test_composition():
    k = 5
    text = 'CAATCCAAC'
    composition = assembly.composition(k, text)
    assert composition == ['AATCC', 'ATCCA', 'CAATC', 'CCAAC', 'TCCAA']


def test_path_to_genome():
    path = ['ACCGA', 'CCGAA', 'CGAAG', 'GAAGC', 'AAGCT']
    genome = assembly.path_to_genome(path)
    assert genome == 'ACCGAAGCT'


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
    kmers = ['GAGG', 'CAGG', 'GGGG', 'GGGA', 'CAGG', 'AGGG', 'GGAG']

    expected = {'AGG': ['GGG'],
                'CAG': ['AGG', 'AGG'],
                'GAG': ['AGG'],
                'GGA': ['GAG'],
                'GGG': ['GGG', 'GGA']}

    adj = assembly.de_brujin_graph_from_kmers(kmers)
    assert adj == expected


def test_de_brujin_graph_from_paired_reads():
    d = 2
    reads = [['GAGA', 'TTGA'],
             ['TCGT', 'GATG'],
             ['CGTG', 'ATGT'],
             ['TGGT', 'TGAG'],
             ['GTGA', 'TGTT'],
             ['GTGG', 'GTGA'],
             ['TGAG', 'GTTG'],
             ['GGTC', 'GAGA'],
             ['GTCG', 'AGAT']]

    adj = assembly.de_brujin_graph_from_paired_reads(reads)
    assert adj == {('GAG', 'TTG'): [('AGA', 'TGA')], ('TCG', 'GAT'): [('CGT', 'ATG')], ('CGT', 'ATG'): [('GTG', 'TGT')],
                   ('TGG', 'TGA'): [('GGT', 'GAG')], ('GTG', 'TGT'): [('TGA', 'GTT')], ('GTG', 'GTG'): [('TGG', 'TGA')],
                   ('TGA', 'GTT'): [('GAG', 'TTG')], ('GGT', 'GAG'): [('GTC', 'AGA')], ('GTC', 'AGA'): [('TCG', 'GAT')]}


def test_reconstruct_string_from_paired_reads():
    d = 2
    reads = [['GAGA', 'TTGA'],
             ['TCGT', 'GATG'],
             ['CGTG', 'ATGT'],
             ['TGGT', 'TGAG'],
             ['GTGA', 'TGTT'],
             ['GTGG', 'GTGA'],
             ['TGAG', 'GTTG'],
             ['GGTC', 'GAGA'],
             ['GTCG', 'AGAT']]
    genome = assembly.reconstruct_string_from_paired_reads(d, reads)
    assert genome == 'GTGGTCGTGAGATGTTGA'

    with open(Path('test/reconstruct_string_from_paired_reads_test.txt')) as input_file:
        k, d = fetch_int_pair_from_file(input_file)
        lines = fetch_lines_from_file(input_file)
    with open(Path('test/reconstruct_string_from_paired_reads_test_out.txt')) as input_file:
        expected = input_file.readline().rstrip()
    reads2 = [line.split('|') for line in lines]
    assert k == len(reads2[0][0])
    genome2 = assembly.reconstruct_string_from_paired_reads(d, reads2)
    assert genome2 == expected


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


def test_is_k_universal():
    assert assembly.is_k_universal(4, '0000110010111101')

    with open(Path('test/kuniversal_test.txt')) as input_file:
        string = input_file.readline().rstrip()
    assert assembly.is_k_universal(14, string)


def text_make_key_universal():
    for k in range(1, 10):
        string = assembly.make_k_universal_string(k)
        assert assembly.is_k_universal(k, string)


def test_reconstruct_string():
    kmers = ['CTTA', 'ACCA', 'TACC', 'GGCT', 'GCTT', 'TTAC']
    assert assembly.reconstruct_string_from_kmers(kmers) == 'GGCTTACCA'

    with open(Path('test/reconstruct_string_test.txt')) as input_file:
        kmers = input_file.readlines()
    kmers = [kmer.rstrip() for kmer in kmers]
    with open(Path('test/reconstruct_string_test_out.txt')) as input_file:
        expected = input_file.readline().rstrip()
    assert assembly.reconstruct_string_from_kmers(kmers) == expected


def test_gapped_path_to_genome():
    gapped_kmers = [['GACC', 'GCGC'],
                    ['ACCG', 'CGCC'],
                    ['CCGA', 'GCCG'],
                    ['CGAG', 'CCGG'],
                    ['GAGC', 'CGGA']]

    gen = assembly.gapped_path_to_genome(2, gapped_kmers)
    assert gen == 'GACCGAGCGCCGGA'


def test_max_no_branch_paths():
    adj = {'1': ['2'],
           '2': ['3'],
           '3': ['4', '5'],
           '6': ['7'],
           '7': ['6']}

    paths = assembly.max_no_branch_paths(adj)
    assert (sorted(paths) == sorted([['1', '2', '3'], ['3', '4'], ['3', '5'], ['6', '7', '6']])) or (
                sorted(paths) == sorted([['1', '2', '3'], ['3', '4'], ['3', '5'], ['7', '6', '7']]))

    adj = {'TA': ['AA'],
           'AA': ['AT'],
           'AT': ['TG', 'TG', 'TG'],
           'TG': ['GC', 'GG', 'GT'],
           'GT': ['TT'],
           'GC': ['CC'],
           'CC': ['CA'],
           'CA': ['AT'],
           'GG': ['GG', 'GA'],
           'GA': ['AT']}
    paths = assembly.max_no_branch_paths(adj)
    assert sorted(paths) == sorted(
        [['TA', 'AA', 'AT'], ['AT', 'TG'], ['AT', 'TG'], ['AT', 'TG'], ['TG', 'GC', 'CC', 'CA', 'AT'], ['TG', 'GG'],
         ['TG', 'GT', 'TT'], ['GG', 'GG'], ['GG', 'GA', 'AT']])


def test_contigs_from_kmers():
    kmers = [
        'ATG',
        'ATG',
        'TGT',
        'TGG',
        'CAT',
        'GGA',
        'GAT',
        'AGA']
    paths = assembly.contigs_from_kmers(kmers)
    assert sorted(paths) == sorted(['ATG', 'ATG', 'TGT', 'TGGA', 'CAT', 'GAT', 'AGA'])


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
    expected2 = {'AAT': ['ATG'], 'ATG': ['TGC', 'TGG', 'TGT'], 'CAT': ['ATG'], 'CCA': ['CAT'], 'GAT': ['ATG'],
                 'GCC': ['CCA'], 'GGA': ['GAT'], 'GGG': ['GGA', 'GGG'], 'GTT': [], 'TAA': ['AAT'], 'TGC': ['GCC'],
                 'TGG': ['GGA', 'GGG'], 'TGT': ['GTT']}
    assert adj == expected2


def test_break_gapped_reads():
    gapped_reads = [['AACCGG', 'CCGGTT'], ['GACTTG', 'CTTTAT'], ['TATATA', 'GCGCGC']]
    res = assembly.break_gapped_reads(2, gapped_reads)
    print(res)
