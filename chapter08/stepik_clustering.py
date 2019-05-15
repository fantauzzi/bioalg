"""
Convenience functions to input/output data in the format expected by Stepik coding challenges
"""


def parse_matrix(input_file):
    lines = input_file.readlines()
    points = []
    for line in lines:
        point = [float(component) for component in line.rstrip('\n').split(' ')]
        points.append(point)
    return points


def fetch_farthest_first_centers_input(file_name):
    with open(file_name) as input_file:
        k, m = input_file.readline().rstrip('\n').split(' ')
        k, m = int(k), int(m)
        points = parse_matrix(input_file)
        assert len(points[0]) == m

    return k, points


def fetch_sq_error_dist_input(file_name):
    with open(file_name) as input_file:
        k, m = input_file.readline().rstrip('\n').split(' ')
        k, m = int(k), int(m)
        lines = input_file.readlines()
        points = []
        centers = []
        for count, line in enumerate(lines):
            if count == k:
                continue
            point = [float(component) for component in line.rstrip('\n').split(' ')]
            assert len(point) == m
            if count < k:
                centers.append(point)
            elif count > k:
                points.append(point)
        return points, centers, k


def fetch_hierarchical_cluster_input(file_name):
    with open(file_name) as input_file:
        n = input_file.readline().rstrip('\n')
        n = int(n)
        matrix = parse_matrix(input_file)
        assert len(matrix) == n
        return matrix


def pretty_print_matrix(m):
    for row in m:
        print(*row, sep=' ')


