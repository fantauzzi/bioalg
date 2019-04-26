def fetch_ints(file_name):
    with open(file_name) as input_file:
        line = input_file.readline().rstrip('\n')
    ints = line.split(' ')
    ints = [int(item) for item in ints]
    return ints


def pretty_print_adj(graph):
    if not graph.edges():
        return
    lines = []
    for edge in sorted(graph.edges()):
        node1, node2 = edge
        ammino = graph[node1][node2]['ammino']
        lines.append(str(node1)+'->'+str(node2)+':'+ammino)
    print(*lines, sep='\n')

