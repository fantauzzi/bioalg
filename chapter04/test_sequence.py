import sequence


def test_translate_RNA():
    rna = 'AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA'
    expected = 'MAMAPRTEINSTRING'
    res = sequence.translate_RNA(rna)
    assert res == expected

    assert sequence.translate_RNA('GAA') == 'E'
    assert sequence.translate_RNA('UAG') == ''


def test_peptide_encoding():
    dna = 'ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA'
    ammino = 'MA'
    res = sequence.peptide_encoding(dna, ammino)
    assert sorted(res) == sorted(['ATGGCC', 'GGCCAT', 'ATGGCC'])
