"""
NOTE
====

Tries in the code below are represented with their adjacency lists, where each node is associated with (is adjacent to)
its children via labelled edges. The adjacency lists are a Python dictionary, with key->value associations, where
'key' is an integer representing the node, and 'value' is another dictionary, with the adjacency information for that
node. Node 0 is always the root. The adjacency of a given node is a dictionary with key_adj->value_adj, where 'key_adj'
is the label for an edge, and 'value_adj' is the node on the other end of the egde. Note that adjacencies only report
children, not the parent node. If a node has no children, then its adjacencies are an empty dictionary {}.
"""


def trie_from_strings(strings):
    """
    Returns the adjacency lists for a trie corresponding to given strings.
    :param strings: The strings, a sequence of strings.
    :return: The trie, represented as its adjacency lists, a sequence of dictionaries.
    """
    trie = {0: {}}  # A trie containing only the root node, numbered as 0
    next_node = 1  # Keep track of the number for the next node to be inserted in the trie
    for string in strings:
        current_node = 0  # The root
        for symbol in string:
            # Is there an edge from current_node labelled with symbol?
            adj_node = trie[current_node].get(symbol)
            # If so, then follow it, and its other end-point becomes the current node
            if adj_node is not None:
                current_node = adj_node
            else:  # Otherwise, add such an edge, and a new node at its other end
                trie[current_node][symbol] = next_node
                trie[next_node] = {}
                current_node = next_node
                next_node += 1
    return trie
