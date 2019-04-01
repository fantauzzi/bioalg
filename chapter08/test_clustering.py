from pathlib import Path
import numpy as np
import clustering
from stepik_clustering import fetch_farthest_first_centers_input, pretty_print_matrix, fetch_sq_error_dist_input


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
