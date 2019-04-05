def pretty_print_trie(trie):
    for node, adjs in sorted(trie.items()):
        for symbol, node2 in adjs.items():
            print(node, '->', node2, ":", symbol, sep='')
