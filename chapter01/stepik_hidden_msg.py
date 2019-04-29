def fetch_clumps_finding_input(file_name):
    with open(file_name) as input_file:
        genome = input_file.readline().rstrip('\n')
        line = input_file.readline().rstrip('\n')
    a, b, c = line.split(' ')
    a, b, c = int(a), int(b), int(c)
    return genome, a, b, c
