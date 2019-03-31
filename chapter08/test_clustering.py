from pathlib import Path
import clustering
from stepik_clustering import fetch_farthest_first_centers_input, pretty_print_matrix


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
