from collections import namedtuple

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

''' The edge leaving a given node, in a modified suffix trie, it provides: the node at the other end of the edge,
and a weight. '''
Edge = namedtuple('Edge', ['node', 'weight'])


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


def prefix_trie_matching(text, trie):
    """
    Given a text and a trie, returns a string encoded by the trie that is a prefix of the text, if such a string exists,
    otherwise returns None.
    :param text: The text, a string.
    :param trie: The trie, a dictionary with its adjacency lists.
    :return: The string from trie which is a prefix for text, if it exists, otherwise None.
    """
    # Generate symbols from text, left to right, one at a time
    text_symbols = (symbol for symbol in text)
    node = 0  # The root of the trie, start looking for a match from there
    match = []  # Will collect the list of matching symbols here
    while trie[node]:  # As long as the currently visited node is not a leaf...
        try:
            symbol = next(text_symbols)  # Get the next symbol from the text
        except StopIteration:  # Handle the case where the text is shorter than any string encoded by 'trie'
            return None
        # Is there an outgoing edge from 'node' labelled with 'symbol'?
        next_node = trie[node].get(symbol)
        ''' If yes, append the symbol to the matching string in 'match', and visit the node at the other end of the edge
        at the next iteration '''
        if next_node is not None:
            match.append(symbol)
            node = next_node
        else:  # If not, then there is no match in trie for the given text
            return None
    return ''.join(match)


def trie_matching(text, trie):
    """
    Given a text and a trie, returns all starting positions in text where a string encoded by the trie appears as a
    substring, and the substring.
    :param text: The text, a string.
    :param trie: The trie, a dictionary of its adajcency lists.
    :return: A list with one element for every match, each element a pair with the matching string as the first item,
    and its position within the text as the second item. Positions are numbered starting from 0. If no match is found,
    the returned list is empty.
    """
    res = []
    for pos in range(0, len(text)):
        match = prefix_trie_matching(text[pos:], trie)
        if match is not None:
            res.append((match, pos))
    return res


def suffix_trie_from_text(text):
    """
    Returns the adjacency lists for the modified suffix trie of a text.
    :param text: The text, a string. It is not allowed to contain symbol '$.
    :return: The modified suffix trie, a dictionary with its adjacency lists.
    """
    assert '$' not in text
    text = text + '$'
    trie = {0: {}}  # Begin with a trie with the root node only; the root is numbered as 0
    next_node = 1  # Keep track of the next node to be inserted in the trie; 0 is the root, already inserted
    # Leaves will be labelled with the starting position in text of the string spelled by the path from root to the leaf
    leaf_labels = {}
    for i in range(0, len(text)):  # Process every suffix of text, starting from the longest
        current_node = 0  # Start from the root
        for j, symbol in enumerate(text[i:]):
            # If there is an edge labelled with 'symbol' leaving the current node, then update its weight and traverse it
            adj = trie[current_node].get(symbol)
            if adj is not None:
                current_node = adj.node
            else:  # If not, add to the trie a new edge (and node) for 'symbol', then traverse it
                trie[current_node][symbol] = Edge(node=next_node, weight=j + i)
                trie[next_node] = {}
                current_node = next_node
                next_node += 1
        # If the current node is a leaf, then set its label
        if not trie[current_node]:
            leaf_labels[current_node] = i
    return trie, leaf_labels


def suffix_tree_from_text(text):
    trie, leaf_labels = suffix_trie_from_text(text)
    # Collect all nodes that begin a branching path, i.e. with a fan-out of at least 2
    branching_nodes = [node for node in trie if len(trie[node]) > 1]
    position = {}
    length = {}
    for branching_node in branching_nodes:
        # Follow each branching path
        for label, edge in trie[branching_node].items():
            current_node = edge.node
            ''' Keep following the path as long as the current node is neither a leaf nor the start of another 
            branching path; remove the nodes within the path from the trie.'''
            path_length = 1
            while len(trie[current_node]) == 1:
                next_node = next(iter(trie[current_node].values())).node
                del trie[current_node]
                current_node = next_node
                path_length += 1
            ''' Add to the trie and edge from the branching node at the beginning of the path, to the last node of the
            path ('current_node') '''
            trie[branching_node][label] = Edge(node=current_node, weight=None)
            position[(branching_node, current_node)] = edge.weight
            length[(branching_node, current_node)] = path_length
    return trie, position, length
