def fetch_farthest_first_centers_input(file_name):
    with open(file_name) as input_file:
        k, m = input_file.readline().rstrip('\n').split(' ')
        k, m = int(k), int(m)
        lines = input_file.readlines()
        points = []
        for line in lines:
            point = [float(component) for component in line.rstrip('\n').split(' ')]
            assert len(point) == m
            points.append(point)

    return k, points


def pretty_print_matrix(m):
    for row in m:
        print(*row, sep=' ')
