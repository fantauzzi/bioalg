from pathlib import Path
import rearrange


def fetch_list_of_lists(file_name):
    """
    Reads from file and returns the expected output for a test case of test_greedy_sorting().
    :param file_name: Name of the file, with relative path, containing the expected output.
    :return: The expected output, a list of lists of integer numbers.
    """

    res = []
    with open(file_name) as input_file:
        line = input_file.readline()
        while line:
            items = [int(item) for item in line.split(' ')]
            res.append(items)
            line = input_file.readline()
    return res


def fetch_stepik_permutations(file_name):
    res = []
    with open(file_name) as input_file:
        line = input_file.readline()
        while line:
            permutations = line.split(')(')
            res.append([])
            for permutation in permutations:
                permutation_as_int = [int(item) for item in permutation.lstrip('(').rstrip(')\n').split(' ')]
                res[-1].append(permutation_as_int)
            line = input_file.readline()
    return res


def test_greedy_sorting():
    p = [-3, +4, +1, +5, -2]
    res = rearrange.greedy_sorting(p)
    assert res == [[-1, -4, +3, +5, -2],
                   [+1, -4, +3, +5, -2],
                   [+1, +2, -5, -3, +4],
                   [+1, +2, +3, +5, +4],
                   [+1, +2, +3, -4, -5],
                   [+1, +2, +3, +4, -5],
                   [+1, +2, +3, +4, +5]]

    p = [-340, -324, +233, -345, +21, +120, -366, -419, +105, -25, +356, -148, +279, +59, -110, +231, -375, +330, +199,
         +48, +122, -42, +236, -365, -26, +275, -60, -227, +259, -338, +183, -126, -223, -354, +93, -288, -221, -332,
         +219, +400, -350, +147, -286, +143, +161, -213, -234, -152, -17, -109, -116, -176, +131, -55, +323, +267, +285,
         +82, +287, -61, -207, -83, -1, -7, -188, +242, -78, -396, -392, +274, -258, -371, -309, +69, +45, +57, +314,
         +228, +373, -90, +159, +348, -239, -257, -136, -63, -144, +162, +104, +278, +406, +182, -170, -158, -169, +265,
         +399, +226, -201, +33, -79, +357, -118, +212, -343, +196, -154, +326, -89, +187, -291, -269, -240, +115, +405,
         -130, +195, +386, +310, +281, +420, +85, +172, +128, -16, +49, +404, -304, +376, +331, +168, +73, +282, +124,
         +302, -232, +121, +398, +355, -74, +146, -247, +165, -75, +362, +5, -395, -253, +41, -316, -205, -56, -24, -37,
         -262, +150, -250, -409, +283, -290, -64, +10, +378, -211, -43, -51, -65, -370, +235, -35, -377, -84, -353, -40,
         -34, -103, +125, +28, -422, +209, -335, +108, +294, +328, +6, +417, +38, -29, +36, -284, +206, +296, +300,
         -100, -312, -237, +166, -76, -384, -244, -208, +245, +397, -140, -401, +119, +261, +299, +292, +347, +194, +80,
         -62, +117, +9, +238, +214, +382, +321, +71, +307, +252, -218, +190, -308, -3, -98, +30, +361, +171, +341, -106,
         +295, +298, +352, +418, -13, -413, -387, -129, -351, -254, -319, -58, -318, -94, +315, +177, -87, -317, +391,
         +181, -180, +416, -132, +402, -163, +151, -408, +97, +142, -77, +46, +67, -27, -342, -364, -173, -210, +134,
         -388, -349, -359, +368, -276, +91, -339, +92, -264, +255, +193, -11, +385, +263, +230, -191, +241, -225, +268,
         -320, +272, -149, -81, -31, -139, -337, +164, -204, +19, -403, -123, -222, +289, +360, -306, -367, +248, +23,
         +39, +380, -260, +135, -224, -22, -107, -383, -246, -12, -369, +297, -344, -99, -305, -280, +44, -412, +113,
         -167, -156, +145, -336, -325, -101, +203, +53, +47, +50, +389, +293, +8, -137, +410, +346, -220, +186, -256,
         +249, +215, +153, +311, +141, -127, -72, -138, +411, -421, -390, -303, -358, -216, +192, -266, +333, +313,
         -277, -270, -32, -102, -189, -2, -68, -414, +329, -334, -271, -198, -114, -174, +393, -155, -327, -175, -273,
         +70, +95, +217, -301, -4, +52, -88, -18, +229, -54, -415, +184, +200, -363, -243, +202, -112, -251, -379, +178,
         -179, -157, +197, -160, +185, -381, -96, +372, -374, -15, +133, -111, -394, +20, -407, +14, +322, -86, +66]

    res = rearrange.greedy_sorting(p)
    expected = fetch_list_of_lists(Path('test/testcase01.txt'))
    assert res == expected


def test_count_breakpoints():
    p = [+10, +6, -8, -7, +17, -20, +18, +19, -5, -16, -11, -4, -3, -2, +13, +14, -1, +9, -12, +15]
    count = rearrange.count_breakpoints(p)
    p = [-16, -20, +11, +12, -14, -13, -15, -6, -8, -19, -18, -17, -10, +4, -5, -2, +7, -3, +1, -9]
    count = rearrange.count_breakpoints(p)

    p = [+3, +4, +5, -12, -8, -7, -6, +1, +2, +10, +9, -11, +13, +14]
    count = rearrange.count_breakpoints(p)
    assert count == 8

    with open('test/testcase02.txt') as input_file:
        line = input_file.readline()
        p = [int(item) for item in line.split(' ')]
    count = rearrange.count_breakpoints(p)
    assert count == 178


def test_chromosome_to_cycle():
    ch = (+1, -2, -3, +4)
    cycle = rearrange.chromosome_to_cycle(ch)
    assert cycle == [1, 2, 4, 3, 6, 5, 7, 8]

    ch = (+1, +2, -3, -4, +5, -6, +7, +8, -9, +10, -11, -12, +13, -14, +15, -16, +17, +18, +19, -20, +21, -22, -23, -24,
          +25, +26, +27, -28, +29, -30, -31, +32, +33, +34, -35, +36, +37, +38, -39, +40, -41, -42, +43, -44, -45, -46,
          -47, +48, +49, -50, +51, -52, -53, -54, +55, +56, +57, -58, -59, -60, +61, +62, +63, -64, +65, +66, +67, -68,
          -69)
    cycle = rearrange.chromosome_to_cycle(ch)
    assert cycle == [1, 2, 3, 4, 6, 5, 8, 7, 9, 10, 12, 11, 13, 14, 15, 16, 18, 17, 19, 20, 22, 21, 24, 23, 25, 26, 28,
                     27, 29, 30,
                     32, 31, 33, 34, 35, 36, 37, 38, 40, 39, 41, 42, 44, 43, 46, 45, 48, 47, 49, 50, 51, 52, 53, 54, 56,
                     55, 57, 58,
                     60, 59, 62, 61, 63, 64, 65, 66, 67, 68, 70, 69, 71, 72, 73, 74, 75, 76, 78, 77, 79, 80, 82, 81, 84,
                     83, 85, 86, 88,
                     87, 90, 89, 92, 91, 94, 93, 95, 96, 97, 98, 100, 99, 101, 102, 104, 103, 106, 105, 108, 107, 109,
                     110, 111, 112,
                     113, 114, 116, 115, 118, 117, 120, 119, 121, 122, 123, 124, 125, 126, 128, 127, 129, 130, 131, 132,
                     133, 134, 136,
                     135, 138, 137]

    ch = (-1, -2, -3, +4, -5, +6, -7, +8, +9, +10, +11, -12, -13, +14, +15, +16, +17, -18, +19, +20, -21, +22, +23, -24,
          -25, -26, +27, +28, +29, +30, -31, +32, -33, -34, -35, +36, -37, +38, +39, +40, -41, -42, +43, +44, -45, +46,
          +47, -48, -49, +50, -51, -52, +53, +54, +55, -56, +57, -58, -59, +60, +61, -62, +63, -64, +65, -66, +67, -68)
    cycle = rearrange.chromosome_to_cycle(ch)
    assert cycle == [2, 1, 4, 3, 6, 5, 7, 8, 10, 9, 11, 12, 14, 13, 15, 16, 17, 18, 19, 20, 21, 22, 24, 23, 26, 25, 27,
                     28, 29, 30, 31, 32, 33, 34, 36, 35, 37, 38, 39, 40, 42, 41, 43, 44, 45, 46, 48, 47, 50, 49, 52, 51,
                     53, 54, 55, 56, 57, 58, 59, 60, 62, 61, 63, 64, 66, 65, 68, 67, 70, 69, 71, 72, 74, 73, 75, 76, 77,
                     78, 79, 80, 82, 81, 84, 83, 85, 86, 87, 88, 90, 89, 91, 92, 93, 94, 96, 95, 98, 97, 99, 100, 102,
                     101, 104, 103, 105, 106, 107, 108, 109, 110, 112, 111, 113, 114, 116, 115, 118, 117, 119, 120, 121,
                     122, 124, 123, 125, 126, 128, 127, 129, 130, 132, 131, 133, 134, 136, 135]


def test_cycle_to_chromosome():
    cycle = [1, 2, 4, 3, 6, 5, 7, 8]
    ch = rearrange.cycle_to_chromosome(cycle)
    assert ch == [1, -2, -3, 4]

    cycle = (
        1, 2, 3, 4, 5, 6, 8, 7, 9, 10, 11, 12, 14, 13, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 26, 25, 27, 28, 29, 30,
        31,
        32, 33, 34, 36, 35, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 48, 47, 49, 50, 52, 51, 54, 53, 56, 55, 57, 58, 60,
        59,
        61, 62, 63, 64, 65, 66, 68, 67, 70, 69, 72, 71, 73, 74, 76, 75, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88,
        89,
        90, 92, 91, 93, 94, 96, 95, 97, 98, 100, 99, 102, 101, 103, 104, 106, 105, 108, 107, 109, 110, 111, 112, 113,
        114,
        115, 116, 118, 117, 120, 119, 122, 121, 123, 124, 125, 126, 128, 127, 130, 129, 132, 131, 133, 134, 136, 135,
        137,
        138, 140, 139)
    ch = rearrange.cycle_to_chromosome(cycle)
    assert ch == [+1, +2, +3, -4, +5, +6, -7, +8, +9, +10, +11, +12, -13, +14, +15, +16, +17, -18, +19, +20, +21, +22,
                  +23, -24, +25, -26, -27, -28, +29, -30, +31, +32, +33, -34, -35, -36, +37, -38, +39, +40, +41, +42,
                  +43, +44, +45, -46, +47, -48, +49, -50, -51, +52, -53, -54, +55, +56, +57, +58, -59, -60, -61, +62,
                  +63, -64, -65, -66, +67, -68, +69, -70]

    cycle = [2, 1, 3, 4, 6, 5, 7, 8, 9, 10, 11, 12, 14, 13, 15, 16, 18, 17, 20, 19, 21, 22, 24, 23, 25, 26, 27, 28, 30,
             29, 31, 32, 33, 34, 36, 35, 38, 37, 40, 39, 41, 42, 44, 43, 46, 45, 48, 47, 49, 50, 51, 52, 54, 53, 55, 56,
             57, 58, 59, 60, 62, 61, 64, 63, 65, 66, 67, 68, 69, 70, 71, 72, 74, 73, 75, 76, 77, 78, 80, 79, 82, 81, 83,
             84, 85, 86, 88, 87, 89, 90, 91, 92, 94, 93, 95, 96, 98, 97, 99, 100, 102, 101, 103, 104, 106, 105, 107,
             108, 110, 109, 112, 111, 114, 113, 116, 115, 117, 118, 119, 120, 122, 121, 123, 124, 126, 125, 127, 128,
             129, 130]
    ch = rearrange.cycle_to_chromosome(cycle)
    assert ch == [-1, 2, -3, 4, 5, 6, -7, 8, -9, -10, 11, -12, 13, 14, -15, 16, 17, -18, -19, -20, 21, -22, -23, -24,
                  25, 26, -27, 28, 29, 30, -31, -32, 33, 34, 35, 36, -37, 38, 39, -40, -41, 42, 43, -44, 45, 46, -47,
                  48, -49, 50, -51, 52, -53, 54, -55, -56, -57, -58, 59, 60, -61, 62, -63, 64, 65]


def test_colored_edges():
    genome = ((+1, -2, -3), (+4, +5, -6))
    edges = rearrange.colored_edges(genome)
    assert edges == [(2, 4), (3, 6), (5, 1), (8, 9), (10, 12), (11, 7)]

    genome = ((+1, -2, -3, +4, -5, -6, -7, -8, +9, -10, -11, -12, +13, +14, +15, -16, -17, +18, +19, +20, +21), (
        -22, -23, -24, -25, -26, -27, -28, +29, +30, +31, -32, -33, +34, +35, -36, +37, +38, +39, -40, +41, +42, +43,
        +44),
              (+45, +46, +47, +48, -49, +50, -51, +52, -53, +54, +55, -56, +57, -58, +59, -60, -61, -62, -63, -64, +65,
               -66), (
                  -67, +68, +69, +70, -71, -72, +73, +74, -75, +76, -77, -78, +79, +80, +81, +82, -83, +84, +85, +86,
                  -87,
                  -88, -89, +90), (
                  -91, -92, +93, +94, -95, -96, +97, +98, -99, -100, -101, -102, +103, +104, -105, +106, -107, -108,
                  +109,
                  +110, +111, -112, +113, +114, -115, +116, +117, -118, +119), (
                  -120, +121, +122, -123, -124, -125, +126, +127, +128, +129, -130, +131, -132, +133, -134, -135, -136,
                  -137, +138, -139, +140, +141, +142, +143, -144, -145, +146, -147, +148), (
                  -149, -150, +151, -152, -153, -154, -155, -156, +157, +158, +159, +160, -161, +162, +163, +164, +165,
                  -166, +167, -168, -169, +170, -171, +172))
    edges = rearrange.colored_edges(genome)
    assert edges == [(2, 4), (3, 6), (5, 7), (8, 10), (9, 12), (11, 14), (13, 16), (15, 17), (18, 20), (19, 22),
                     (21, 24), (23, 25), (26, 27), (28, 29), (30, 32), (31, 34), (33, 35), (36, 37), (38, 39), (40, 41),
                     (42, 1), (43, 46), (45, 48), (47, 50), (49, 52), (51, 54), (53, 56), (55, 57), (58, 59), (60, 61),
                     (62, 64), (63, 66), (65, 67), (68, 69), (70, 72), (71, 73), (74, 75), (76, 77), (78, 80), (79, 81),
                     (82, 83), (84, 85), (86, 87), (88, 44), (90, 91), (92, 93), (94, 95), (96, 98), (97, 99),
                     (100, 102), (101, 103), (104, 106), (105, 107), (108, 109), (110, 112), (111, 113), (114, 116),
                     (115, 117), (118, 120), (119, 122), (121, 124), (123, 126), (125, 128), (127, 129), (130, 132),
                     (131, 89), (133, 135), (136, 137), (138, 139), (140, 142), (141, 144), (143, 145), (146, 147),
                     (148, 150), (149, 151), (152, 154), (153, 156), (155, 157), (158, 159), (160, 161), (162, 163),
                     (164, 166), (165, 167), (168, 169), (170, 171), (172, 174), (173, 176), (175, 178), (177, 179),
                     (180, 134), (181, 184), (183, 185), (186, 187), (188, 190), (189, 192), (191, 193), (194, 195),
                     (196, 198), (197, 200), (199, 202), (201, 204), (203, 205), (206, 207), (208, 210), (209, 211),
                     (212, 214), (213, 216), (215, 217), (218, 219), (220, 221), (222, 224), (223, 225), (226, 227),
                     (228, 230), (229, 231), (232, 233), (234, 236), (235, 237), (238, 182), (239, 241), (242, 243),
                     (244, 246), (245, 248), (247, 250), (249, 251), (252, 253), (254, 255), (256, 257), (258, 260),
                     (259, 261), (262, 264), (263, 265), (266, 268), (267, 270), (269, 272), (271, 274), (273, 275),
                     (276, 278), (277, 279), (280, 281), (282, 283), (284, 285), (286, 288), (287, 290), (289, 291),
                     (292, 294), (293, 295), (296, 240), (297, 300), (299, 301), (302, 304), (303, 306), (305, 308),
                     (307, 310), (309, 312), (311, 313), (314, 315), (316, 317), (318, 319), (320, 322), (321, 323),
                     (324, 325), (326, 327), (328, 329), (330, 332), (331, 333), (334, 336), (335, 338), (337, 339),
                     (340, 342), (341, 343), (344, 298)]


def test_graph_to_genome():
    graph = [(2, 4), (3, 6), (5, 1), (7, 9), (10, 12), (11, 8)]
    genome = rearrange.graph_to_genome(graph)
    assert genome == [[1, -2, -3], [-4, 5, -6]]

    graph = [(1, 2), (3, 6), (5, 7), (8, 4)]
    genome = rearrange.graph_to_genome(graph)
    assert genome == [[-1], [-2, -3, 4]]

    graph = [(1, 4), (3, 6), (5, 7), (8, 9), (10, 11), (12, 14), (13, 15), (16, 18), (17, 19), (20, 22), (21, 24),
             (23, 26), (25, 27), (28, 29), (30, 32), (31, 33), (34, 36), (35, 38), (37, 39), (40, 42), (41, 44),
             (43, 2), (45, 48), (47, 50), (49, 51), (52, 54), (53, 56), (55, 57), (58, 60), (59, 62), (61, 63),
             (64, 66), (65, 68), (67, 69), (70, 72), (71, 74), (73, 76), (75, 78), (77, 79), (80, 81), (82, 84),
             (83, 86), (85, 87), (88, 90), (89, 91), (92, 94), (93, 46), (96, 98), (97, 99), (100, 102), (101, 103),
             (104, 106), (105, 107), (108, 109), (110, 111), (112, 114), (113, 115), (116, 117), (118, 119), (120, 122),
             (121, 123), (124, 125), (126, 127), (128, 130), (129, 131), (132, 134), (133, 135), (136, 138), (137, 139),
             (140, 141), (142, 144), (143, 145), (146, 147), (148, 149), (150, 152), (151, 95), (154, 155), (156, 158),
             (157, 159), (160, 161), (162, 163), (164, 166), (165, 168), (167, 169), (170, 171), (172, 174), (173, 175),
             (176, 178), (177, 180), (179, 182), (181, 184), (183, 186), (185, 187), (188, 190), (189, 191), (192, 194),
             (193, 196), (195, 153), (198, 200), (199, 201), (202, 204), (203, 205), (206, 207), (208, 210), (209, 211),
             (212, 214), (213, 215), (216, 217), (218, 219), (220, 222), (221, 223), (224, 226), (225, 227), (228, 230),
             (229, 232), (231, 234), (233, 235), (236, 238), (237, 239), (240, 241), (242, 244), (243, 245), (246, 197),
             (248, 250), (249, 252), (251, 253), (254, 255), (256, 258), (257, 260), (259, 262), (261, 263), (264, 265),
             (266, 267), (268, 270), (269, 272), (271, 273), (274, 276), (275, 278), (277, 280), (279, 281), (282, 283),
             (284, 286), (285, 247), (288, 290), (289, 292), (291, 293), (294, 295), (296, 297), (298, 299), (300, 302),
             (301, 303), (304, 306), (305, 307), (308, 310), (309, 312), (311, 314), (313, 316), (315, 318), (317, 319),
             (320, 322), (321, 323), (324, 325), (326, 327), (328, 330), (329, 332), (331, 333), (334, 336), (335, 338),
             (337, 339), (340, 287)]
    genome = rearrange.graph_to_genome(graph)
    assert genome == [
        [-1, -2, -3, +4, +5, +6, -7, +8, -9, +10, -11, -12, -13, +14, +15, -16, +17, -18, -19, +20, -21, -22],
        [-23, -24, -25, +26, -27, -28, +29, -30, -31, +32, -33, -34, +35, -36, -37, -38, -39, +40, +41, -42, -43, +44,
         -45, +46, -47],
        [+48, -49, +50, -51, +52, -53, +54, +55, +56, -57, +58, +59, +60, -61, +62, +63, +64, -65, +66, -67, +68, -69,
         +70, +71, -72, +73, +74, +75, -76],
        [+77, +78, -79, +80, +81, +82, -83, -84, +85, +86, -87, +88, -89, -90, -91, -92, -93, +94, -95, +96, -97, -98],
        [+99, -100, +101, -102, +103, +104, -105, +106, -107, +108, +109, +110, -111, +112, -113, +114, -115, -116,
         -117, +118, -119, +120, +121, -122, +123],
        [+124, -125, -126, +127, +128, -129, -130, -131, +132, +133, +134, -135, -136, +137, -138, -139, -140, +141,
         +142, -143],
        [+144, -145, -146, +147, +148, +149, +150, -151, +152, -153, +154, -155, -156, -157, -158, -159, +160, -161,
         +162, +163, +164, -165, -166, +167, -168, -169, +170]]


def test_two_break_on_genome_graph():
    graph = [(2, 4), (3, 8), (7, 5), (6, 1)]
    res = rearrange.two_break_on_genome_graph(graph, 1, 6, 3, 8)
    assert res == [(2, 4), (3, 1), (7, 5), (6, 8)]

    res = rearrange.two_break_on_genome_graph(graph, 6, 1, 8, 3)
    assert sorted(res) == sorted([(2, 4), (3, 1), (7, 5), (6, 8)])

    res = rearrange.two_break_on_genome_graph(graph, 3, 8, 1, 6)
    assert sorted(res) == sorted([(2, 4), (3, 1), (7, 5), (6, 8)])

    res = rearrange.two_break_on_genome_graph(graph, 8, 3, 6, 1)
    assert res == [(2, 4), (3, 1), (7, 5), (6, 8)]

    graph = [(2, 4), (3, 6), (5, 7), (8, 1)]
    res = rearrange.two_break_on_genome_graph(graph, 1, 8, 3, 6)
    assert res == [(2, 4), (3, 1), (5, 7), (8, 6)]

    graph = [(2, 3), (4, 8), (7, 6), (5, 9), (10, 1)]
    res = rearrange.two_break_on_genome_graph(graph, 1, 10, 5, 9)
    assert res == [(2, 3), (4, 8), (7, 6), (5, 1), (10, 9)]

    res2 = rearrange.two_break_on_genome_graph(graph, 9, 5, 10, 1)
    assert res2 == res

    res = rearrange.two_break_on_genome_graph(graph, 1, 10, 7, 6)
    assert res == [(2, 3), (4, 8), (7, 1), (5, 9), (10, 6)]

    graph = [(2, 4), (3, 5), (6, 8), (7, 10), (9, 12), (11, 13), (14, 15), (16, 17), (18, 19), (20, 21), (22, 24),
             (23, 25), (26, 28), (27, 30), (29, 32), (31, 33), (34, 35), (36, 37), (38, 40), (39, 42), (41, 44),
             (43, 45), (46, 47), (48, 49), (50, 51), (52, 53), (54, 56), (55, 58), (57, 59), (60, 61), (62, 64),
             (63, 66), (65, 68), (67, 69), (70, 72), (71, 73), (74, 75), (76, 77), (78, 80), (79, 81), (82, 83),
             (84, 86), (85, 88), (87, 90), (89, 92), (91, 94), (93, 96), (95, 98), (97, 99), (100, 102), (101, 104),
             (103, 106), (105, 107), (108, 109), (110, 111), (112, 113), (114, 116), (115, 117), (118, 120), (119, 121),
             (122, 124), (123, 126), (125, 128), (127, 129), (130, 132), (131, 134), (133, 136), (135, 1)]

    res = rearrange.two_break_on_genome_graph(graph, 87, 90, 74, 75)
    assert res == [(2, 4), (3, 5), (6, 8), (7, 10), (9, 12), (11, 13), (14, 15), (16, 17), (18, 19), (20, 21), (22, 24),
                   (23, 25), (26, 28), (27, 30), (29, 32), (31, 33), (34, 35), (36, 37), (38, 40), (39, 42), (41, 44),
                   (43, 45), (46, 47), (48, 49), (50, 51), (52, 53), (54, 56), (55, 58), (57, 59), (60, 61), (62, 64),
                   (63, 66), (65, 68), (67, 69), (70, 72), (71, 73), (74, 87), (76, 77), (78, 80), (79, 81), (82, 83),
                   (84, 86), (85, 88), (75, 90), (89, 92), (91, 94), (93, 96), (95, 98), (97, 99), (100, 102),
                   (101, 104), (103, 106), (105, 107), (108, 109), (110, 111), (112, 113), (114, 116), (115, 117),
                   (118, 120), (119, 121), (122, 124), (123, 126), (125, 128), (127, 129), (130, 132), (131, 134),
                   (133, 136), (135, 1)]


def test_fix_graph_cycles():
    graph = [(2, 4), (3, 8), (7, 5), (6, 1)]
    broken = rearrange.two_break_on_genome_graph(graph, 1, 6, 3, 8)
    fixed, n_cycles = rearrange.fix_graph_cycles(broken)
    assert fixed == [(1, 3), (4, 2), (5, 7), (8, 6)]
    assert n_cycles == 2

    graph = [(2, 4), (3, 5), (6, 8), (7, 10), (9, 12), (11, 13), (14, 15), (16, 17), (18, 19), (20, 21), (22, 24),
             (23, 25), (26, 28), (27, 30), (29, 32), (31, 33), (34, 35), (36, 37), (38, 40), (39, 42), (41, 44),
             (43, 45), (46, 47), (48, 49), (50, 51), (52, 53), (54, 56), (55, 58), (57, 59), (60, 61), (62, 64),
             (63, 66), (65, 68), (67, 69), (70, 72), (71, 73), (74, 87), (76, 77), (78, 80), (79, 81), (82, 83),
             (84, 86), (85, 88), (75, 90), (89, 92), (91, 94), (93, 96), (95, 98), (97, 99), (100, 102),
             (101, 104), (103, 106), (105, 107), (108, 109), (110, 111), (112, 113), (114, 116), (115, 117),
             (118, 120), (119, 121), (122, 124), (123, 126), (125, 128), (127, 129), (130, 132), (131, 134),
             (133, 136), (135, 1)]

    fixed, n_cycles = rearrange.fix_graph_cycles(graph)
    assert n_cycles == 1


def test_two_break_on_genome():
    genome = [[1, -2, -4, 3], [5, -6, -8, 7]]
    res = rearrange.two_break_on_genome(genome, 1, 6, 3, 8)
    assert res == [[-1, 2], [-3, 4], [-5, -7, 8, 6]]

    # genome = [[+1, +2, +3, +4],[+5, +6],[+7,+8,+9]]
    # res = rearrange.two_break_on_genome(genome, )

    genome = [
        [-1, +2, -3, -4, -5, -6, -7, +8, +9, +10, +11, +12, +13, +14, +15, -16, -17, -18, +19, -20, +21, +22, -23, -24,
         +25, -26, -27, +28, -29, -30, +31, +32, +33, -34, -35, +36, +37, +38, +39, -40, +41, +42, -43, +44, +45, +46,
         -47, +48, +49, +50, -51, +52, +53, +54, -55, -56, -57, +58, -59, -60, -61, +62, -63, +64, +65, -66, -67]]
    res = rearrange.two_break_on_genome(genome, 116, 118, 73, 72)
    assert res == [
        [-1, 2, -3, -4, -5, -6, -7, 8, 9, 10, 11, 12, 13, 14, 15, -16, -17, -18, 19, -20, 21, 22, -23, -24, 25, -26,
         -27, 28, -29, -30, 31, 32, 33, -34, -35, 36, -59, -60, -61, 62, -63, 64, 65, -66, -67],
        [-37, -58, 57, 56, 55, -54, -53, -52, 51, -50, -49, -48, 47, -46, -45, -44, 43, -42, -41, 40, -39, -38]]

    genome = [(-1, -2, -3, 4, 5, -6, -7, -8, 9, -10, 11, 12, 13, 14, -15, 16, -17, -18, 19, 20, -21, 22, -23, 24, -25,
               -26, -27, 28, 29, -30, 31, 32, -33, -34, 35, -36, 37, 38, 39, 40, -41, 42, -43, 44, -45, -46, -47, -48,
               49, -50, 51, -52, -53, 54, 55, -56, 57, 58, -59, 60, -61)]
    res = rearrange.two_break_on_genome(genome, 4, 1, 108, 109)
    assert res == [[-1, 55, -56, 57, 58, -59, 60, -61],
                   [-2, -3, 4, 5, -6, -7, -8, 9, -10, 11, 12, 13, 14, -15, 16, -17, -18, 19, 20, -21, 22, -23, 24, -25,
                    -26, -27, 28, 29, -30, 31, 32, -33, -34, 35, -36, 37, 38, 39, 40, -41, 42, -43, 44, -45, -46, -47,
                    -48, 49, -50, 51, -52, -53, 54]]


def test_two_break_dist():
    ps = [(1, 2, 3, 4, 5, 6)]
    qs = [(2, -4), (1, -3, -6, -5)]
    dist = rearrange.two_break_distance(ps, qs)
    assert dist == 3

    permutations = fetch_stepik_permutations(Path('test/testcase03.txt'))
    assert len(permutations) == 2
    ps = permutations[0]
    qs = permutations[1]
    dist = rearrange.two_break_distance(ps, qs)
    assert dist == 8079

    permutations = fetch_stepik_permutations(Path('test/testcase04.txt'))
    assert len(permutations) == 2
    ps = permutations[0]
    qs = permutations[1]
    dist = rearrange.two_break_distance(ps, qs)
    assert dist == 8354


def test_two_break_sorting():
    ps = [[1, -2, -3, 4]]
    qs = [[1, 2, -4, -3]]
    res = rearrange.two_break_sorting(ps, qs)
    assert res == [[[1, -2, -3, 4]], [[-1], [-2, -3, 4]], [[-1, -3, 4, -2]], [[-1, 3, 4, -2]]]

    ps = [[+9, -8, +12, +7, +1, -14, +13, +3, -5, -11, +6, -2, +10, -4]]
    qs = [[-11, +8, -10, -2, +3, +4, +13, +6, +12, +9, +5, +7, -14, -1]]
    res = rearrange.two_break_sorting(ps, qs)
    assert res == [[[9, -8, 12, 7, 1, -14, 13, 3, -5, -11, 6, -2, 10, -4]],
                   [[-1], [-2, 10, -4, 9, -8, 12, 7, -14, 13, 3, -5, -11, 6]],
                   [[-1, 13, 3, -5, -11, 6, -2, 10, -4, 9, -8, 12, 7, -14]],
                   [[-1, -10, 2, -6, 11, 5, -3, -13, -4, 9, -8, 12, 7, -14]],
                   [[-1, -9, 4, 13, 3, -5, -11, 6, -2, 10, -8, 12, 7, -14]],
                   [[-1, -11, 6, -2, 10, -8, 12, 7, -14], [-3, -13, -4, 9, 5]],
                   [[-1, -11, 6, 2, 10, -8, 12, 7, -14], [-3, -13, -4, 9, 5]],
                   [[-1, -11, 6, 12, 7, -14], [-2, 8, -10], [-3, -13, -4, 9, 5]],
                   [[-1, -11, 8, -10, -2, 6, 12, 7, -14], [-3, -13, -4, 9, 5]],
                   [[-1, -11, 8, -10, -2, 3, -5, -9, 4, 13, 6, 12, 7, -14]],
                   [[-1, -11, 8, -10, -2, 3, -12, -6, -13, -4, 9, 5, 7, -14]],
                   [[-1, -11, 8, -10, -2, 3, 4, 13, 6, 12, 9, 5, 7, -14]]]

    ps = [(-9, -1, -2, -6, +8, -12, -11, +4, -7, +10, -5, +3)]
    qs = [(+2, +9, -4, +8, +12, +1, +6, -3, +7, +11, +5, -10)]
    res = rearrange.two_break_sorting(ps, qs)
    assert res == [[(-9, -1, -2, -6, 8, -12, -11, 4, -7, 10, -5, 3)],
                   [[-1], [-2, -6, 8, -12, -11, 4, -7, 10, -5, 3, -9]],
                   [[-1, 8, -12, -11, 4, -7, 10, -5, 3, -9, -2, -6]], [[-1, 11, 12, -8, 4, -7, 10, -5, 3, -9, -2, -6]],
                   [[-1, -4, 8, -12, -11, -7, 10, -5, 3, -9, -2, -6]], [[-1, -3, 5, -10, 7, 11, 12, -8, 4, -9, -2, -6]],
                   [[-1, 2, 9, -4, 8, -12, -11, -7, 10, -5, 3, -6]], [[-1, 7, 11, 12, -8, 4, -9, -2, 10, -5, 3, -6]],
                   [[-1, 5, -10, 2, 9, -4, 8, -12, -11, -7, 3, -6]], [[-1, 12, -8, 4, -9, -2, 10, -5, -11, -7, 3, -6]],
                   [[-1, -12, -8, 4, -9, -2, 10, -5, -11, -7, 3, -6]]]


def test_find_shared_kmers():
    s1 = 'TCAGTTGGCCTACAT'
    s2 = 'CCTACATGAGGTCTG'
    pos = rearrange.find_shared_kmers(3, s1, s2)
    assert sorted(pos) == sorted([(8, 0), (9, 1), (10, 2), (11, 3), (12, 4), (12, 5), (0, 6), (8, 8), (1, 12)])

    s1 = 'AAACTCATC'
    s2 = 'TTTCAAATC'
    pos = rearrange.find_shared_kmers(3, s1, s2)
    assert sorted(pos) == sorted([(0, 4), (0, 0), (4, 2), (6, 6)])

    with open(Path('test/testcase05.txt')) as input_file:
        s1 = input_file.readline().rstrip('\n')
        s2 = input_file.readline().rstrip('\n')
    pos = rearrange.find_shared_kmers(17, s1, s2)
    assert sorted(pos) == sorted(
        [(120069, 237), (120070, 238), (195106, 238), (120071, 239), (195107, 239), (52667, 401), (162091, 401),
         (162090, 402), (3522, 424), (3521, 425), (71719, 425), (14521, 425), (152918, 686), (32618, 687),
         (149215, 687), (152919, 687), (32619, 688), (4462, 914), (4463, 915), (198123, 931), (159694, 932),
         (198124, 932), (198125, 933), (84508, 1046), (76663, 1046), (27099, 1046), (27100, 1047), (168848, 1162),
         (36531, 1484), (89268, 1484), (36530, 1485), (89267, 1485), (15761, 1485), (32619, 2034), (32618, 2035),
         (149215, 2035), (152919, 2035), (115465, 2462), (136933, 2661), (144396, 3189), (144395, 3190), (6310, 3190),
         (28253, 3219), (28254, 3220), (28255, 3221), (58146, 3328), (133968, 3329), (58145, 3329), (133967, 3330),
         (58144, 3330), (52970, 3330), (68429, 3525), (68428, 3526), (68427, 3527), (188490, 3527), (201061, 3527),
         (52667, 3607), (162091, 3607), (162090, 3608), (167039, 3685), (74888, 3686), (19381, 3686), (167040, 3686),
         (49567, 3789), (49566, 3790), (49565, 3791), (135498, 3791), (6504, 3792), (49564, 3792), (135497, 3792),
         (84510, 4405), (76665, 4405), (84509, 4406), (76664, 4406), (84508, 4407), (76663, 4407), (27099, 4407),
         (179922, 4694), (94600, 4695), (179923, 4695), (94601, 4696), (179924, 4696), (176648, 4696), (144395, 5255),
         (6310, 5255), (6309, 5256), (34545, 5494), (84318, 5590), (42954, 5689), (146189, 5689), (169271, 5689),
         (14501, 5944), (19075, 5944), (120395, 6227), (14268, 6227), (51480, 6422), (23261, 6422), (180945, 6423),
         (51479, 6423), (23260, 6423), (177700, 6485), (67641, 6485), (96780, 6485), (67640, 6486), (96779, 6486),
         (146187, 7073), (42953, 7074), (146188, 7074), (42954, 7075), (146189, 7075), (169271, 7075), (11409, 7153),
         (87154, 7153), (22956, 7263), (22957, 7264), (144498, 7264), (58629, 7264), (22958, 7265), (144499, 7265),
         (15760, 7531), (36530, 7532), (89267, 7532), (15761, 7532), (94600, 8107), (179923, 8107), (94601, 8108),
         (179924, 8108), (176648, 8108), (185037, 8134), (185038, 8135), (84508, 8314), (76663, 8314), (27099, 8314),
         (84507, 8315), (76662, 8315), (27098, 8315), (76661, 8316), (20477, 8529), (171653, 8529), (11409, 9010),
         (87154, 9010), (146187, 9381), (42953, 9382), (146188, 9382), (42954, 9383), (146189, 9383), (169271, 9383),
         (58630, 9540), (22957, 9541), (144498, 9541), (58629, 9541), (58628, 9542), (40141, 9953), (174984, 9953),
         (50951, 9953), (95836, 10214), (81907, 10760), (115249, 10760), (25795, 10760), (94600, 11523),
         (179923, 11523), (94601, 11524), (179924, 11524), (176648, 11524), (99487, 11559), (101465, 11559),
         (101466, 11560), (52969, 11700), (133967, 11701), (58144, 11701), (52970, 11701), (133968, 11702),
         (58145, 11702), (58146, 11703), (58147, 11704), (58148, 11705), (96781, 11968), (177700, 11969),
         (67641, 11969), (96780, 11969), (167041, 12019), (74888, 12020), (19381, 12020), (167040, 12020),
         (19380, 12021), (19379, 12022), (104474, 12275), (4465, 12275), (104475, 12276), (4466, 12276), (4467, 12277),
         (4468, 12278), (4469, 12279), (3521, 12334), (71719, 12334), (14521, 12334), (71720, 12335), (22959, 12493),
         (22958, 12494), (144499, 12494), (22957, 12495), (144498, 12495), (58629, 12495), (150822, 12610),
         (150823, 12611), (32618, 12642), (149215, 12642), (152919, 12642), (137746, 13232), (137745, 13233),
         (52609, 13233), (137744, 13234), (52608, 13234), (159694, 13445), (198124, 13445), (159695, 13446),
         (137745, 13519), (52609, 13519), (51478, 13588), (23259, 13588), (180945, 13589), (51479, 13589),
         (23260, 13589), (51480, 13590), (23261, 13590), (68427, 14064), (188490, 14064), (201061, 14064),
         (92272, 14208), (92271, 14209), (92270, 14210), (92269, 14211), (166516, 14284), (36529, 14460),
         (36530, 14461), (89267, 14461), (15761, 14461), (15762, 14462), (167039, 14527), (74888, 14528),
         (19381, 14528), (167040, 14528), (74889, 14529), (180945, 14676), (51479, 14676), (23260, 14676),
         (94219, 14879), (94220, 14880), (49283, 14880), (174983, 15320), (40141, 15321), (174984, 15321),
         (50951, 15321), (174985, 15322), (50952, 15322), (174986, 15323), (94220, 15845), (49283, 15845),
         (94219, 15846), (94218, 15847), (81906, 15969), (115248, 15969), (81907, 15970), (115249, 15970),
         (25795, 15970), (40141, 16026), (174984, 16026), (50951, 16026), (120070, 16076), (195106, 16076),
         (120069, 16077), (120068, 16078), (65305, 16235), (7748, 16236), (65306, 16236), (7749, 16237), (7750, 16238),
         (164125, 16354), (164124, 16355), (1652, 16355), (188491, 16369), (201062, 16369), (171656, 16948),
         (171655, 16949), (171654, 16950), (20477, 16951), (171653, 16951), (20476, 16952), (3522, 17131),
         (3521, 17132), (71719, 17132), (14521, 17132), (6504, 17707), (49564, 17707), (135497, 17707), (164124, 17976),
         (1652, 17976), (164123, 17977), (1651, 17977), (164122, 17978), (6504, 18851), (49564, 18851), (135497, 18851),
         (49563, 18852), (135496, 18852), (52969, 19423), (133967, 19424), (58144, 19424), (52970, 19424),
         (81907, 19460), (115249, 19460), (25795, 19460), (120396, 19696), (120395, 19697), (14268, 19697),
         (7747, 20094), (7748, 20095), (65306, 20095), (65307, 20096), (99486, 20374), (99487, 20375), (101465, 20375),
         (101466, 20376)]
    )


def xtest_graph_from_permutations():
    p = (1, -2, -3, 4)
    graph = rearrange.breakpoint_from_permutations([p])
    assert graph == {2: 4, 3: 6, 5: 7, 8: 1}

    q = (1, 3, 2, -4)
    graph = rearrange.breakpoint_from_permutations([q])
    assert graph == {2: 5, 6: 3, 4: 8, 7: 1}

    p = [(1, -2), (-4, 3)]
    graph = rearrange.breakpoint_from_permutations(p)
    assert graph == {2: 4, 3: 1, 7: 5, 6: 8}

    p = [(2, -4), (1, -3, -6, -5)]
    graph = rearrange.breakpoint_from_permutations(p)
    assert graph == {4: 8, 7: 3, 2: 6, 5: 12, 11: 10, 9: 1}

    p = [(1, 2, 3, 4, 5, 6)]
    graph = rearrange.breakpoint_from_permutations(p)
    assert graph == {2: 3, 4: 5, 6: 7, 8: 9, 10: 11, 12: 1}


def xtest_breakpoint_graph():
    q = [(2, -4), (1, -3, -6, -5)]
    p = [(1, 2, 3, 4, 5, 6)]
    adj = rearrange.breakpoint_graph(p, q, 'red', 'blue')
    assert adj == {2: [(3, 'red'), (6, 'blue')],
                   3: [(2, 'red'), (7, 'blue')],
                   4: [(5, 'red'), (8, 'blue')],
                   5: [(4, 'red'), (12, 'blue')],
                   6: [(7, 'red'), (2, 'blue')],
                   7: [(6, 'red'), (3, 'blue')],
                   8: [(9, 'red'), (4, 'blue')],
                   9: [(8, 'red'), (1, 'blue')],
                   10: [(11, 'red'), (11, 'blue')],
                   11: [(10, 'red'), (10, 'blue')],
                   12: [(1, 'red'), (5, 'blue')],
                   1: [(12, 'red'), (9, 'blue')]}

    p = (1, -2, -3, 4)
    q = (1, 3, 2, -4)
    adj = rearrange.breakpoint_graph([p], [q], color_ps='red', color_qs='blue')
    assert adj == {2: [(4, 'red'), (5, 'blue')],
                   4: [(2, 'red'), (8, 'blue')],
                   3: [(6, 'red'), (6, 'blue')],
                   6: [(3, 'red'), (3, 'blue')],
                   5: [(7, 'red'), (2, 'blue')],
                   7: [(5, 'red'), (1, 'blue')],
                   8: [(1, 'red'), (4, 'blue')],
                   1: [(8, 'red'), (7, 'blue')]}


def xtest_permutations_from_breakpoint():
    p = [1, -2, -3, 4]
    adj = rearrange.breakpoint_from_permutations([p])
    p2 = rearrange.permutations_from_breakpoint(adj)
    assert p == p2

    q = [1, 3, 2, -4]
    adj = rearrange.breakpoint_from_permutations([q])
    q2 = rearrange.permutations_from_breakpoint(adj)
    assert q == q2

    p = [[1, -2], [-4, 3]]
    adj = rearrange.breakpoint_from_permutations(p)
    p2 = rearrange.permutations_from_breakpoint(adj)
    assert p2 == [[1, -2], [3, -4]]


def xtest_two_break():
    p = [[1, -2, -4, 3]]
    adj = rearrange.breakpoint_from_permutations(p)
    rearrange.two_break(adj, 1, 6, 3, 8)
    assert adj == {2: 4, 3: 1, 7: 5, 6: 8}


def xtest_two_break_on_breakpoints():
    ps = [[1, -2, -3, 4]]
    qs = [[1, 2, -4, -3]]
    adj = rearrange.breakpoint_graph(ps, qs)
    rearrange.two_break_on_breakpoints(adj, 8, 1, 4, 2, 'red')
    assert adj == {2: [(3, 'blue'), (1, 'red')],
                   4: [(8, 'blue'), (8, 'red')],
                   3: [(6, 'red'), (2, 'blue')],
                   6: [(3, 'red'), (7, 'blue')],
                   5: [(7, 'red'), (1, 'blue')],
                   7: [(5, 'red'), (6, 'blue')],
                   8: [(4, 'blue'), (4, 'red')],
                   1: [(5, 'blue'), (2, 'red')]}

    rearrange.two_break_on_breakpoints(adj, 2, 1, 3, 6, 'red')
    assert adj == {2: [(3, 'blue'), (3, 'red')],
                   4: [(8, 'blue'), (8, 'red')],
                   3: [(2, 'blue'), (2, 'red')],
                   6: [(7, 'blue'), (1, 'red')],
                   5: [(7, 'red'), (1, 'blue')],
                   7: [(5, 'red'), (6, 'blue')],
                   8: [(4, 'blue'), (4, 'red')],
                   1: [(5, 'blue'), (6, 'red')]}

    rearrange.two_break_on_breakpoints(adj, 6, 1, 7, 5, 'red')
    assert adj == {2: [(3, 'blue'), (3, 'red')],
                   4: [(8, 'blue'), (8, 'red')],
                   3: [(2, 'blue'), (2, 'red')],
                   6: [(7, 'blue'), (7, 'red')],
                   5: [(1, 'blue'), (1, 'red')],
                   7: [(6, 'blue'), (6, 'red')],
                   8: [(4, 'blue'), (4, 'red')],
                   1: [(5, 'blue'), (5, 'red')]}

    ps = [[1, -2, -3, 4]]
    qs = [[1, 2, -4, -3]]
    res = rearrange.shortest_rearrangement(ps, qs)
    pass
