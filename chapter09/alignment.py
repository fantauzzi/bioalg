from collections import namedtuple, deque

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


"""
NOTE
====
In the rest of the code, a trie or tree is a group of SuffixNode instances, linked to each other via their attributes
'parent' and 'symbol_to_child'.
"""


class SuffixNode:
    """
    The node for a text suffix trie or tree. It may have multiple children; the edge connecting to a child may have
    an associated symbol, weight, starting position and length. It may also have a label, an integer.
    The instance has the following attributes:
        data: the data contained by the node.
        parent: the parent of the node; another SuffixNode, or None, if the node has no parent.
            symbol_to_child: the edges leaving from the node, each edge has an associated symbol and enters another
            node, a dictionary symbol -> SuffixNode.
        weights: a dictionary that gives the weight of outgoing edges with given symbol, symbol -> weight.
        label: a label for the node, an integer or None.
        position: a dictionary that gives the starting position in the text of the prefix beginning at the current node
            and with a given symbol, symbol -> starting position (a non-negative integer).
        length:a dictionary giving the length of in the text of the prefix beginning at the current node and with a
            given symbol, symbol -> length (a non-negative integer)
    """

    def __init__(self, data, parent=None):
        """
        :param data: The data to be contained in the node.
        :param parent:The parent of the node, another SuffixNode, or None if the node has no parent.
        """
        self.data = data
        self.parent = parent
        self.symbol_to_child = {}
        self.weights = {}
        self.label = None
        self.position = {}
        self.length = {}
        self.color = None


def visit_trie_level_order(root):
    """
    Returns a generator that enumerates nodes in a suffix trie or tree in level-order (breadth-first).
    :param root: The root of the trie/tree, a SuffixNode.
    :return: The generator.
    """
    to_be_visited = deque([root])
    while to_be_visited:
        current_node = to_be_visited.popleft()
        for child_node in current_node.symbol_to_child.values():
            to_be_visited.append(child_node)
        yield current_node


def suffix_trie_from_text(text):
    """
    Returns the the modified suffix trie for a given text.
    :param text: The text, a string. It is not allowed to contain symbol '$'.
    :return: The root of the modified suffix trie, a SuffixNode.
    """
    assert '$' not in text
    text = text + '$'
    root = SuffixNode(data=0)  # Begin with a trie with the root node only; the root is numbered as 0
    next_node = 1  # Keep track of the next node number to be inserted in the trie; 0 is the root, already inserted
    for i in range(0, len(text)):  # Process every suffix of text, starting from the longest
        current_node = root  # Start from the root
        for j, symbol in enumerate(text[i:]):
            # If there is an edge labelled with 'symbol' leaving the current node, then update its weight and traverse it
            adj = current_node.symbol_to_child.get(symbol)
            if adj is not None:
                current_node = adj
            else:  # If not, add to the trie a new edge (and node) for 'symbol', then traverse it
                newly_inserted_node = SuffixNode(data=next_node, parent=current_node)
                current_node.symbol_to_child[symbol] = newly_inserted_node
                current_node.weights[symbol] = j + i
                current_node = newly_inserted_node
                next_node += 1
        # If the current node is a leaf, then set its label
        if not current_node.symbol_to_child:
            current_node.label = i
    return root


def suffix_tree_from_text(text):
    """
    Returns the suffix tree for a given text.
    :param text: The text, a string. It is not allowed to contain symbol '$'.
    :return: The root of the suffix tree, a SuffixNode.
    """
    tree_root = suffix_trie_from_text(text)
    # Collect all nodes that begin a branching path, i.e. with at least 2 children
    branching_nodes = [node for node in visit_trie_level_order(tree_root) if len(node.symbol_to_child) > 1]
    for branching_node in branching_nodes:
        # Follow each branching path
        for symbol, current_node in branching_node.symbol_to_child.items():
            ''' Keep following the path as long as the current node is neither a leaf nor the start of another 
            branching path; remove from the tree the nodes within the path.'''
            path_length = 1
            while len(current_node.symbol_to_child) == 1:
                next_node = next(iter(current_node.symbol_to_child.values()))
                del current_node
                current_node = next_node
                path_length += 1
            ''' Add to the tree an edge from the branching node at the beginning of the path, to the last node of the
            path, 'current_node' '''
            branching_node.symbol_to_child[symbol] = current_node
            current_node.parent = branching_node
            branching_node.length[symbol] = path_length
            branching_node.position[symbol] = branching_node.weights[symbol]
            del branching_node.weights[symbol]
    return tree_root


def longest_repeat(text):
    """
    Returns the longest repeat in a text.
    :param text: The text, a string.
    :return: The longest repeat within the given text, a string.
    """

    def longest_path_to_internal(node, text):
        """
        Returns the string corresponding to the longest path in a suffix tree from a given node to its descendants
        that are internal (i.e. non-leaf) nodes.
        :param node:The given node, a SuffixNode.
        :param text:The text which originated the suffix tree.
        :return: The substring of 'text' corresponding to the longest path from the given node to its non-leaf
        descendants, a string. If multiple such substrings exist, only one of them is returned.
        """
        # If you don't find anything better, the returned longest path will be the empty string
        longest_path_symbols = ''
        # Find if any non-leaf child of 'node' leads to a longer path
        for symbol, child in node.symbol_to_child.items():
            if not child.symbol_to_child:  # If the child is a leaf, skip it
                continue
            symbols = text[
                      node.position[symbol]: node.position[symbol] + node.length[symbol]] + longest_path_to_internal(
                child, text)
            if len(symbols) > len(longest_path_symbols):
                longest_path_symbols = symbols
        return longest_path_symbols

    root = suffix_tree_from_text(text)
    # Find the longest path from the root to any branching node, and that corresponds to the longest repeat in 'text'.
    symbols = longest_path_to_internal(root, text)
    return symbols


def color_suffix_tree(root, no_of_blue_leaves):
    ripe_nodes = [node for node in visit_trie_level_order(root) if not node.symbol_to_child]
    ripe_nodes = deque(ripe_nodes)
    count = {'red': 0, 'blue': 0, 'purple': 0}
    while ripe_nodes:
        current_node = ripe_nodes.popleft()
        color = None
        if not current_node.symbol_to_child:
            color = 'blue' if current_node.label < no_of_blue_leaves else 'red'
        for child in current_node.symbol_to_child.values():
            assert child.color is not None
            if color == None:
                color = child.color
            elif color != child.color:
                color = 'purple'
                break
        current_node.color = color
        count[color] += 1
        if current_node.parent is None:
            continue
        for child_of_parent in current_node.parent.symbol_to_child.values():
            if child_of_parent.color is None:
                break
        else:
            ripe_nodes.append(current_node.parent)
    return count

# TODO change suffix_trie_from_text and suffix_tree_from_text to assume the input string already ends by '$'