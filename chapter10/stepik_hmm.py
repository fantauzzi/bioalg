from hmm import HMM


def parse_and_split(input_file):
    line = input_file.readline().rstrip('\n')
    items = line.split()
    return items


def fetch_hmm(file_name):
    def parse_matrix(input_file):
        columns = parse_and_split(input_file)
        matrix = {}
        while True:
            items = parse_and_split(input_file)
            if len(items) <= 1:  # Found the end of the matrix or of the file
                break
            row = items[0]
            matrix[row] = {}
            for col, item in zip(columns, items[1:]):
                matrix[row][col] = float(item)

        return matrix

    with open(file_name) as input_file:
        string = input_file.readline().rstrip('\n')
        input_file.readline()
        alphabet = parse_and_split(input_file)
        input_file.readline()
        states = parse_and_split(input_file)
        input_file.readline()
        transision_matrix = parse_matrix(input_file)
        emission_matrix = parse_matrix(input_file)
        return string, HMM(alphabet=alphabet, states=states, transition=transision_matrix, emission=emission_matrix)


def fetch_alignment(file_name):
    with open(file_name) as input_file:
        params = parse_and_split(input_file)
        theta, sigma = (float(params[0]), float(params[1])) if len(params) == 2 else (float(params[0]), .0)
        input_file.readline()
        alphabet = parse_and_split(input_file)
        input_file.readline()
        alignment = input_file.readlines()
        alignment = [line.rstrip('\n') for line in alignment]
        return theta, sigma, alphabet, alignment


def ugly_print_matrix(matrix, row_labels, col_labels):  # It is ugly, but it is the format stepik/Rosalind expects
    def make_print_label(label):
        n = len(label)
        print_label = 'S' if label[0] == 'S' and n == 2 else 'E' if label[0] == 'E' and n == 2 else label if n == 1 else \
        label[0] + str(label[1])
        return print_label

    print_labels = [make_print_label(label) for label in col_labels if label[0]!='S']
    if ('S', None) in col_labels:
        print('S   ', end='')
    else:
        print('    ', end='')
    print(*print_labels, sep='  ')
    for row_label in row_labels:
        row_print_label = make_print_label(row_label)
        print(row_print_label, '  ', sep='', end='')
        if len(row_print_label) == 1:
            print(' ', end='')
        values = [matrix[row_label][col_label] for col_label in col_labels]
        f = ['{:0.3f}   ' if value != 0 else '{:0.0f}  ' for value in values]
        f = ''.join(f)
        f = f[:len(f) - 2]
        print(f.format(*values))


def ugly_print_matrices(transition, emission, states, alphabet):
    ugly_print_matrix(transition, states, states)
    print('--------')
    ugly_print_matrix(emission, states, alphabet)
