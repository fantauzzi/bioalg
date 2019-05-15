from pathlib import Path
import numpy as np
import clustering
from stepik_clustering import fetch_farthest_first_centers_input, pretty_print_matrix, fetch_sq_error_dist_input, \
    fetch_hierarchical_cluster_input


def testfarthest_first_centers():
    data_points = [[0.0, 0.0],
                   [5.0, 5.0],
                   [0.0, 5.0],
                   [1.0, 1.0],
                   [2.0, 2.0],
                   [3.0, 3.0],
                   [1.0, 2.0]]
    centers = clustering.farthest_first_centers(3, data_points)
    assert centers == [[0., 0.], [5., 5.], [0., 5.]]

    k, data_points = fetch_farthest_first_centers_input(Path('test/testcase01.txt'))
    centers = clustering.farthest_first_centers(k, data_points)
    assert centers == [[0.8, 12.0, 17.5, 0.9, 7.2], [0.3, 16.4, 8.9, 34.6, 24.6], [32.3, 1.9, 5.1, 16.2, 8.8],
                       [23.1, 31.1, 3.6, 0.8, 0.3]]

    k, data_points = fetch_farthest_first_centers_input(Path('test/testcase02.txt'))
    centers = clustering.farthest_first_centers(k, data_points)
    assert centers == [[5.4, 21.0, 12.9], [36.5, 4.7, 6.9], [21.9, 2.4, 28.5], [12.9, 0.3, 0.4]]


def test_sq_error_distortion():
    data_points, centers, k = fetch_sq_error_dist_input(Path('test/testcase12.txt'))
    dist = clustering.sq_error_distortion(data_points, centers)
    print()
    print(dist)

    data_points, centers, k = fetch_sq_error_dist_input(Path('test/testcase03.txt'))
    dist = clustering.sq_error_distortion(data_points, centers)
    assert np.isclose(dist, 18.245559999999994)

    data_points, centers, k = fetch_sq_error_dist_input(Path('test/testcase04.txt'))
    dist = clustering.sq_error_distortion(data_points, centers)
    assert np.isclose(dist, 36.762834331337366)

    data_points, centers, k = fetch_sq_error_dist_input(Path('test/testcase05.txt'))
    dist = clustering.sq_error_distortion(data_points, centers)
    assert np.isclose(dist, 28.69774363476738)


def test_lloyd_cluster():
    k, data_points = fetch_farthest_first_centers_input(Path('test/testcase06.txt'))
    centers = clustering.lloyd_cluster(data_points=data_points, centers=data_points[:k])
    assert np.all(np.isclose(centers, [[1.800, 2.867], [1.060, 1.140]], rtol=0, atol=1e-3))

    k, data_points = fetch_farthest_first_centers_input(Path('test/testcase07.txt'))
    centers = clustering.lloyd_cluster(data_points=data_points, centers=data_points[:k])
    expected = [[7.561403508771929, 6.1672514619883065, 16.568421052631578, 6.078362573099413, 7.0964912280701755],
                [18.23220338983051, 6.14745762711864, 5.4677966101694935, 6.577966101694915, 6.0533898305084755],
                [7.03664596273292, 17.298757763975154, 6.927329192546586, 5.495031055900623, 7.027950310559],
                [7.7124999999999995, 7.232812500000005, 6.916406250000001, 18.717187499999998, 6.838281249999997],
                [6.041666666666666, 6.278846153846154, 5.708333333333335, 7.014102564102566, 17.408333333333328],
                [5.15814814814815, 4.55925925925926, 5.112592592592593, 5.143703703703705, 4.718888888888891]]
    assert np.all(np.isclose(centers, expected, rtol=0, atol=1e-3))

    k, data_points = fetch_farthest_first_centers_input(Path('test/testcase08.txt'))
    centers = clustering.lloyd_cluster(data_points=data_points, centers=data_points[:k])
    expected = [[4.307142857142859, 10.86111111111111, 16.488888888888887],
                [5.606140350877197, 11.50438596491228, 4.507894736842106],
                [17.24605809128631, 7.417842323651452, 6.301244813278007],
                [8.157777777777778, 20.985555555555557, 7.939999999999995],
                [4.922340425531912, 3.6191489361702125, 4.70345744680851],
                [8.710439560439553, 3.365384615384615, 14.840659340659354]]
    assert np.all(np.isclose(centers, expected, rtol=0, atol=1e-3))


def test_hierarchical_clusters():
    d = [[0, 20, 9, 11],
         [20, 0, 17, 11],
         [9, 17, 0, 8],
         [11, 11, 8, 0]]
    c1=(0, 2)
    c2=(1, 3)
    n = len(d)
    d = {(row,): {(col,): d[row][col] for col in range(0, n)} for row in range(0, n)}
    the_dist = clustering.dist(d, c1, c2)
    print()
    print(the_dist)

    d = [[0.00, 0.74, 0.85, 0.54, 0.83, 0.92, 0.89],
         [0.74, 0.00, 1.59, 1.35, 1.20, 1.48, 1.55],
         [0.85, 1.59, 0.00, 0.63, 1.13, 0.69, 0.73],
         [0.54, 1.35, 0.63, 0.00, 0.66, 0.43, 0.88],
         [0.83, 1.20, 1.13, 0.66, 0.00, 0.72, 0.55],
         [0.92, 1.48, 0.69, 0.43, 0.72, 0.00, 0.80],
         [0.89, 1.55, 0.73, 0.88, 0.55, 0.80, 0.00]]

    steps = clustering.hierarchical_clusters(d)
    assert steps == [(5, 3), (4, 6), (2, 5, 3), (0, 1), (2, 5, 3, 4, 6), (0, 1, 2, 5, 3, 4, 6)]

    d = fetch_hierarchical_cluster_input(Path('test/testcase09.txt'))
    steps = clustering.hierarchical_clusters(d)
    assert steps == [(0, 17), (8, 2), (3, 5), (10, 4), (14, 3, 5), (1, 0, 17), (15, 19), (9, 13), (6, 18), (8, 2, 12),
                     (11, 1, 0, 17), (11, 1, 0, 17, 7), (16, 15, 19), (6, 18, 16, 15, 19), (14, 3, 5, 8, 2, 12),
                     (9, 13, 10, 4), (11, 1, 0, 17, 7, 14, 3, 5, 8, 2, 12),
                     (6, 18, 16, 15, 19, 11, 1, 0, 17, 7, 14, 3, 5, 8, 2, 12),
                     (6, 18, 16, 15, 19, 11, 1, 0, 17, 7, 14, 3, 5, 8, 2, 12, 9, 13, 10, 4)]

    d = fetch_hierarchical_cluster_input(Path('test/testcase10.txt'))
    steps = clustering.hierarchical_clusters(d)
    assert steps == [(15, 3), (5, 16), (13, 7), (0, 19), (2, 17), (1, 5, 16), (15, 3, 13, 7), (4, 6), (11, 2, 17),
                     (8, 12), (4, 6, 10), (1, 5, 16, 9), (15, 3, 13, 7, 0, 19), (15, 3, 13, 7, 0, 19, 18),
                     (8, 12, 4, 6, 10), (11, 2, 17, 14), (15, 3, 13, 7, 0, 19, 18, 11, 2, 17, 14),
                     (1, 5, 16, 9, 15, 3, 13, 7, 0, 19, 18, 11, 2, 17, 14),
                     (1, 5, 16, 9, 15, 3, 13, 7, 0, 19, 18, 11, 2, 17, 14, 8, 12, 4, 6, 10)]


def test_max_distance():
    # data = [(2, 8), (2, 5), (6, 9), (7, 5), (5, 2) ]
    data = [(2, 6), (4, 9), (5, 7), (6, 5), (8, 3)]
    # centers = [(3, 5), (5, 4)]
    centers = [(4, 5), (7, 4)]
    dist, center = clustering.max_distance(data, centers)
    print()
    print(dist, center)


def test_hidden_matrix():
    data = [(2, 8), (2, 5), (6, 9), (7, 5), (5, 2)]
    centers = [(3, 5), (5, 4)]
    h_m = clustering.hidden_matrix(data, centers)
    print()
    print(*h_m, sep='\n')
