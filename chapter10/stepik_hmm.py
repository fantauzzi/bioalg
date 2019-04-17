from hmm import HMM


def fetch_hmm(file_name):
    def parse_and_split(input_file):
        line = input_file.readline().rstrip('\n')
        items = line.split()
        return items

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
