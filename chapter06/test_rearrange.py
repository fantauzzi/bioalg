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
    p = [+3, +4, +5, -12, -8, -7, -6, +1, +2, +10, +9, -11, +13, +14]
    count = rearrange.count_breakpoints(p)
    assert count == 8

    with open('test/testcase02.txt') as input_file:
        line = input_file.readline()
        p = [int(item) for item in line.split(' ')]
    count = rearrange.count_breakpoints(p)
    assert count == 178


def test_graph_from_permutations():
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


def test_breakpoint_graph():
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


def test_two_break_dist():
    ps = [(1, 2, 3, 4, 5, 6)]
    qs = [(2, -4), (1, -3, -6, -5)]
    dist = rearrange.two_break_dist(ps, qs)
    assert dist == 3

    permutations = fetch_stepik_permutations(Path('test/testcase03.txt'))
    assert len(permutations) == 2
    ps = permutations[0]
    qs = permutations[1]
    dist = rearrange.two_break_dist(ps, qs)
    assert dist == 8079

    permutations = fetch_stepik_permutations(Path('test/testcase04.txt'))
    assert len(permutations) == 2
    ps = permutations[0]
    qs = permutations[1]
    dist = rearrange.two_break_dist(ps, qs)
    assert dist == 8354


def test_permutations_from_breakpoint():
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


def test_two_break():
    p = [[1, -2, -4, 3]]
    adj = rearrange.breakpoint_from_permutations(p)
    rearrange.two_break(adj, 1, 6, 3, 8)
    assert adj == {2: 4, 3: 1, 7: 5, 6: 8}


def test_two_break_on_breakpoints():
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


def test_shortest_rearrangement():
    ps = [[1, -2, -3, 4]]
    qs = [[1, 2, -4, -3]]
    res = rearrange.shortest_rearrangement(ps, qs)
    pass
