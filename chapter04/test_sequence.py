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


def test_cyclic_spectrum():
    peptide = 'LEQN'
    res = sequence.peptide_spectrum(peptide)
    assert res == [0, 113, 114, 128, 129, 227, 242, 242, 257, 355, 356, 370, 371, 484]

    peptide2 = 'YYTMDKDVWRAP'
    res2 = sequence.peptide_spectrum(peptide2)
    assert res2 == [0, 71, 97, 99, 101, 115, 115, 128, 131, 156, 163, 163, 168, 186, 214, 227, 232, 243, 243, 246, 260,
                    264, 285, 324, 326, 331, 342, 342, 347, 358, 374, 395, 400, 413, 423, 427, 441, 457, 475, 487, 489,
                    494, 510, 510, 512, 524, 528, 556, 558, 588, 590, 595, 609, 627, 638, 643, 650, 655, 673, 673, 684,
                    689, 724, 726, 751, 753, 755, 770, 772, 774, 799, 801, 836, 841, 852, 852, 870, 875, 882, 887, 898,
                    916, 930, 935, 937, 967, 969, 997, 1001, 1013, 1015, 1015, 1031, 1036, 1038, 1050, 1068, 1084, 1098,
                    1102, 1112, 1125, 1130, 1151, 1167, 1178, 1183, 1183, 1194, 1199, 1201, 1240, 1261, 1265, 1279,
                    1282, 1282, 1293, 1298, 1311, 1339, 1357, 1362, 1362, 1369, 1394, 1397, 1410, 1410, 1424, 1426,
                    1428, 1454, 1525]


def test_linear_spectrum():
    peptide = 'NQEL'
    res = sequence.peptide_spectrum(peptide, cyclic=False)
    assert res == [0, 113, 114, 128, 129, 242, 242, 257, 370, 371, 484]

    res2 = sequence.peptide_spectrum('VNLMFIITTVLEWGFTKAFHDLNCMNKVGRFPSVYGGCRKTSHESG', cyclic=False)
    assert res2==[0,57,57,57,57,57,71,87,87,87,97,99,99,99,99,101,101,101,101,103,103,113,113,113,113,113,114,114,114,114,115,128,128,128,129,129,131,131,137,137,144,147,147,147,147,156,156,156,160,163,184,186,186,188,199,200,202,204,212,213,213,214,216,217,217,218,220,224,226,227,227,227,228,229,229,234,242,242,243,244,244,245,248,252,259,260,262,266,273,277,278,283,284,284,284,300,301,303,305,312,313,315,315,316,316,319,325,326,327,330,331,341,341,342,346,348,348,349,353,353,355,358,360,365,372,373,373,373,376,376,380,385,387,390,391,391,398,399,400,406,410,414,414,428,428,430,433,440,440,442,444,445,446,447,447,453,454,457,457,459,461,462,463,470,472,472,474,476,479,479,483,485,487,488,491,497,501,503,504,504,504,505,512,519,527,527,527,529,536,541,543,544,545,554,556,560,566,575,575,575,575,576,582,582,583,584,584,586,587,590,593,594,598,598,602,604,605,609,617,618,619,620,626,628,632,632,632,635,640,643,643,650,651,656,663,664,669,674,684,685,685,689,689,690,690,697,699,701,703,706,707,711,712,713,717,718,722,726,729,729,731,731,731,733,738,742,746,748,749,763,765,769,769,771,786,787,788,788,798,800,802,805,806,806,810,812,818,819,819,819,825,825,826,827,830,832,832,832,832,837,841,842,846,850,852,859,860,861,863,863,864,870,882,885,898,899,902,903,905,916,917,918,918,920,926,928,928,929,931,931,932,933,933,933,935,947,951,951,955,955,955,959,960,962,966,966,966,974,974,974,984,985,985,989,1012,1015,1016,1016,1019,1023,1029,1031,1031,1032,1032,1032,1033,1034,1038,1042,1042,1045,1046,1047,1048,1049,1054,1059,1061,1073,1079,1088,1089,1090,1094,1099,1102,1102,1103,1111,1115,1118,1119,1122,1122,1130,1130,1131,1132,1135,1145,1146,1147,1147,1147,1159,1159,1160,1160,1162,1162,1173,1173,1175,1176,1178,1179,1195,1201,1202,1204,1205,1216,1217,1218,1218,1233,1233,1233,1233,1244,1250,1250,1258,1259,1260,1261,1262,1267,1272,1272,1274,1274,1275,1277,1278,1278,1279,1282,1290,1301,1304,1304,1306,1307,1307,1315,1316,1329,1331,1331,1332,1335,1346,1346,1346,1351,1361,1364,1364,1373,1374,1380,1381,1388,1391,1392,1400,1401,1402,1403,1406,1406,1407,1408,1414,1414,1416,1419,1419,1421,1430,1437,1438,1438,1444,1445,1445,1448,1457,1459,1460,1461,1478,1485,1488,1493,1495,1495,1495,1495,1501,1507,1511,1517,1517,1520,1531,1534,1535,1538,1543,1545,1548,1548,1549,1550,1550,1552,1558,1558,1559,1560,1561,1575,1594,1598,1598,1606,1606,1606,1608,1609,1613,1616,1630,1632,1632,1632,1635,1635,1644,1648,1648,1651,1651,1657,1658,1661,1664,1664,1665,1666,1679,1692,1697,1704,1705,1705,1712,1714,1722,1722,1723,1729,1731,1737,1743,1745,1745,1745,1749,1753,1754,1758,1760,1760,1761,1762,1763,1765,1779,1779,1780,1791,1792,1792,1793,1816,1825,1836,1837,1844,1848,1848,1850,1857,1858,1858,1859,1859,1860,1860,1861,1861,1861,1864,1868,1880,1882,1884,1890,1891,1891,1893,1905,1906,1915,1917,1918,1921,1940,1944,1947,1948,1958,1962,1964,1967,1971,1972,1973,1974,1981,1983,1988,1992,1992,1996,1997,2004,2005,2005,2007,2008,2020,2021,2034,2043,2045,2063,2064,2065,2070,2075,2075,2077,2077,2078,2085,2093,2096,2097,2102,2104,2104,2105,2106,2109,2111,2118,2121,2132,2133,2133,2134,2135,2136,2144,2162,2184,2188,2189,2190,2192,2192,2206,2206,2207,2207,2210,2210,2224,2224,2232,2232,2233,2233,2233,2234,2246,2248,2249,2249,2249,2251,2263,2289,2291,2295,2297,2307,2319,2320,2320,2320,2321,2325,2333,2335,2335,2336,2346,2347,2348,2348,2361,2362,2363,2363,2364,2377,2380,2380,2390,2412,2421,2423,2423,2433,2434,2434,2435,2445,2448,2450,2451,2454,2462,2462,2466,2466,2476,2476,2477,2480,2491,2493,2508,2511,2511,2524,2534,2537,2546,2547,2549,2549,2561,2563,2564,2568,2568,2575,2579,2579,2579,2580,2590,2590,2592,2594,2597,2604,2609,2625,2647,2650,2660,2663,2671,2677,2678,2680,2680,2686,2689,2689,2693,2693,2696,2697,2707,2707,2708,2710,2711,2717,2728,2754,2760,2765,2767,2776,2776,2790,2792,2794,2807,2808,2808,2811,2815,2822,2824,2824,2826,2827,2833,2839,2864,2873,2875,2877,2883,2884,2891,2895,2902,2904,2907,2909,2914,2923,2938,2938,2939,2940,2952,2955,2959,2962,2976,2978,2995,2996,2996,3004,3012,3020,3020,3032,3033,3037,3038,3043,3049,3051,3053,3056,3066,3070,3077,3091,3095,3106,3108,3113,3117,3120,3133,3139,3143,3151,3152,3156,3161,3165,3165,3167,3177,3190,3196,3198,3199,3200,3204,3222,3240,3248,3253,3255,3262,3264,3264,3264,3280,3297,3298,3299,3303,3305,3312,3321,3327,3337,3349,3351,3353,3354,3356,3378,3386,3395,3406,3409,3410,3411,3411,3428,3440,3450,3457,3466,3466,3467,3477,3482,3496,3508,3512,3515,3523,3523,3525,3539,3541,3553,3553,3570,3580,3581,3595,3610,3613,3613,3622,3624,3628,3640,3640,3652,3652,3670,3683,3694,3709,3721,3726,3727,3727,3739,3741,3741,3744,3765,3781,3796,3801,3808,3808,3828,3830,3839,3842,3854,3857,3858,3864,3868,3894,3907,3914,3925,3929,3955,3961,3965,3967,3971,3971,3981,3986,3993,4028,4038,4042,4066,4068,4070,4074,4080,4085,4094,4114,4117,4127,4137,4155,4179,4181,4184,4188,4195,4215,4230,4238,4245,4282,4287,4292,4302,4308,4339,4344,4346,4358,4395,4421,4433,4439,4443,4452,4459,4472,4508,4546,4565,4568,4570,4571,4573,4655,4660,4672,4683,4699,4712,4759,4786,4797,4812,4843,4896,4899,4926,4956,5013,5025,5070,5112,5169]

def print_cyclopeptide_sequences(seqs):
    # ammino_mass = sequence.get_ammino_mass()
    # masses = [[ammino_mass[ammino] for ammino in seq] for seq in seqs]
    for item in seqs:
        print(*item, sep='-', end='')
        print(' ', sep='', end='')
    print()

def test_sequence_cyclopeptide():
    spectrum = [0, 97, 97, 99, 101, 103, 196, 198, 198, 200, 202, 295, 297, 299, 299, 301, 394, 396, 398, 400, 400, 497]
    res = sequence.sequence_cyclopeptide(spectrum)
    print()
    print_cyclopeptide_sequences(res)

    spectrum2 = [0, 113, 128, 186, 241, 299, 314, 427]
    res2 = sequence.sequence_cyclopeptide(spectrum2)
    print_cyclopeptide_sequences(res2)

    spectrum3 = [0,87,99,113,115,115,137,137,147,147,202,228,234,236,246,252,252,260,284,339,347,349,351,365,375,383,383,399,438,462,462,486,498,498,512,512,520,577,585,585,599,599,611,635,635,659,698,714,714,722,732,746,748,750,758,813,837,845,845,851,861,863,869,895,950,950,960,960,982,982,984,998,1010,1097]
    res3 = sequence.sequence_cyclopeptide(spectrum3)
    print_cyclopeptide_sequences(res3)




