def pretty_print_trie(trie):
    for node, adjs in sorted(trie.items()):
        for symbol, node2 in adjs.items():
            print(node, '->', node2, ":", symbol, sep='')


def pretty_print_edge_labels(labels):
    for label in labels:
        if not label:
            continue
        if label[-1] == 'x':
            label = label[:len(label) - 1] + '$'
        print(label)
