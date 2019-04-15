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
    terminator = text[-1]
    assert terminator not in text[0:len(text) - 1]
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

    root = suffix_tree_from_text(text + '$')
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


def longest_shared_substring(string1, string2):
    def longest_path_to_purple(node, text):
        """
        Returns the string corresponding to the longest path in a suffix tree from a given node to its descendants
        that are colored as purple.
        :param node:The given node, a SuffixNode.
        :param text:The text which originated the suffix tree.
        :return: The substring of 'text' corresponding to the longest path from the given node to its purple
        descendants, a string. If multiple such substrings exist, only one of them is returned.
        """
        # If you don't find anything better, the returned longest path will be the empty string
        longest_path_symbols = ''
        # Find if any non-leaf child of 'node' leads to a longer path
        for symbol, child in node.symbol_to_child.items():
            if child.color != 'purple':  # If child is not purple, skip it.
                continue
            symbols = text[
                      node.position[symbol]: node.position[symbol] + node.length[symbol]] + longest_path_to_purple(
                child, text)
            if len(symbols) > len(longest_path_symbols):
                longest_path_symbols = symbols
        return longest_path_symbols

    assert '#' not in string1 + string2
    assert '$' not in string1 + string2
    text = string1 + '#' + string2 + '$'
    root = suffix_tree_from_text(text)
    color_suffix_tree(root, len(string1))
    substring = longest_path_to_purple(root, text)
    return substring


def shortes_substring_not_appearing(string1, string2):
    def longest_match_in_prefix_tree(root, text, string):
        matched = ''
        current_node = root
        mismatch = False
        ''' Stop traversing the tree if: 'current_node' is a leaf, the whole string to be matched has been matched, or
        a mismatch has been found '''
        while current_node.symbol_to_child and len(matched) < len(string) and not mismatch:
            # Find an edge leaving from 'current_node' whose symbol is the next symbol to be matched
            symbol = string[len(matched)]
            child = current_node.symbol_to_child.get(symbol)
            # If it doesn't exist, then there is a mismatch, stop traversing the tree
            if child is None:
                break
            matching = string[
                       len(matched):len(matched) + min((len(string) - len(matched), current_node.length[symbol]))]
            for i in range(0, len(matching)):
                try:
                    if matching[i] != text[current_node.position[symbol] + i]:
                        mismatch = True
                        matched = matched + matching[:i]
                        break
                except IndexError:
                    pass
            else:
                matched = matched + matching
            current_node = child
        return matched

    assert '$' not in string1 + string2
    string2 = string2 + '$'
    root = suffix_tree_from_text(string2)
    shortes_longest_match = None
    shortes_longest_match_pos = None
    for i in range(0, len(string1)):
        suffix = string1[i:]
        longest_match = longest_match_in_prefix_tree(root, string2, suffix)
        longest_match = longest_match.rstrip('$')
        if shortes_longest_match is None or (
                len(longest_match) < len(shortes_longest_match) and len(longest_match) + i < len(string1)):
            shortes_longest_match = longest_match
            shortes_longest_match_pos = i
    shortest_mismatch = string1[shortes_longest_match_pos: shortes_longest_match_pos + len(shortes_longest_match) + 1]
    return shortest_mismatch


def suffix_array_for_text(text):
    suffixes = {}
    for i in range(0, len(text)):
        suffix = text[i:]
        suffixes[suffix] = i

    res = [value for (_, value) in sorted(suffixes.items())]
    return res


def burrows_wheeler_transform(text):
    rotations = [text[i:] + text[:i] for i in range(0, len(text))]
    rotations.sort()
    bw_transformed = [item[-1] for item in rotations]
    return ''.join(bw_transformed)


def inverted_burrow_wheeler(text):
    assert text.count('$') == 1
    sorted_text = sorted(text)
    left_count = [0] * len(text)
    left_symbols_count = {}
    for pos, symbol in enumerate(sorted_text):
        count = left_symbols_count.get(symbol, -1)
        count += 1
        left_count[pos] = count
        left_symbols_count[symbol] = count

    right_symbols_count = {}
    right_count = {}
    for pos, symbol in enumerate(text):
        count = right_symbols_count.get(symbol, -1)
        count += 1
        right_symbols_count[symbol] = count
        right_count[(symbol, count)] = pos

    inverted = []
    current_pos = text.find('$')
    assert current_pos >= 0
    while len(inverted) < len(text):
        current_left_symbol = sorted_text[current_pos]
        inverted.append(current_left_symbol)
        current_pos = right_count[(current_left_symbol, left_count[current_pos])]

    return ''.join(inverted)


def count_matches(text, pattern):
    def last_to_first(transformed_text):
        assert transformed_text.count('$') == 1
        sorted_text = sorted(transformed_text)

        left_symbols_count = {}  # TODO refactor common code and choose between left/right and first/last nomenclature
        left_symbols = {}
        for pos, symbol in enumerate(sorted_text):
            count = left_symbols_count.get(symbol, -1)
            count += 1
            left_symbols_count[symbol] = count
            left_symbols[(symbol, count)] = pos

        l_to_f = [0] * len(transformed_text)
        right_symbols_count = {}
        for pos, symbol in enumerate(transformed_text):
            count = right_symbols_count.get(symbol, -1)
            count += 1
            right_symbols_count[symbol] = count
            l_to_f[pos] = left_symbols[(symbol, count)]
        return l_to_f

    if not isinstance(pattern, str):
        counts = [count_matches(text, one_pattern) for one_pattern in pattern]
        return counts

    transformed_text = text
    l_to_f = last_to_first(transformed_text)

    top, bottom = 0, len(text) - 1
    n = len(transformed_text)
    while top <= bottom:
        if pattern:
            symbol = pattern[-1]
            pattern = pattern[:len(pattern) - 1]
            try:
                top_index = transformed_text.index(symbol, top, bottom + 1)
            except ValueError:
                return 0
            bottom_index = n - transformed_text[::-1].index(symbol, n - 1 - bottom, n - top) - 1
            top = l_to_f[top_index]
            bottom = l_to_f[bottom_index]
        else:
            return bottom - top + 1
    assert False


def better_BW_matching(first_occurrence, last_column, pattern, count):
    top = 0
    bottom = len(last_column) - 1
    while top <= bottom:
        if pattern:
            symbol = pattern[-1]
            pattern = pattern[:len(pattern) - 1]
            if symbol in last_column[top: bottom + 1]:
                top = first_occurrence[symbol] + count[symbol][top]
                bottom = first_occurrence[symbol] + count[symbol][bottom + 1] - 1
            else:
                return 0
        else:
            return bottom - top + 1


def better_count_matches(transformed_text, patterns):
    first_column = sorted(transformed_text)

    first_occurrence = {}
    for pos, symbol in enumerate(first_column):
        if first_occurrence.get(symbol) is None:
            first_occurrence[symbol] = pos

    symbols = set(transformed_text)
    count = {symbol: [0] * (len(transformed_text) + 1) for symbol in symbols}
    for pos, pos_symbol in enumerate(transformed_text):
        for symbol in symbols:
            count[symbol][pos + 1] = count[symbol][pos] + 1 if symbol == pos_symbol else count[symbol][pos]

    counters = []
    for pattern in patterns:
        the_counter = better_BW_matching(first_occurrence, transformed_text, pattern, count)
        counters.append(the_counter)

    return counters


def matching_positions(first_occurrence, last_column, pattern, count):
    top = 0
    bottom = len(last_column) - 1
    while top <= bottom:
        if pattern:
            symbol = pattern[-1]
            pattern = pattern[:len(pattern) - 1]
            if symbol in last_column[top: bottom + 1]:
                top = first_occurrence[symbol] + count[symbol][top]
                bottom = first_occurrence[symbol] + count[symbol][bottom + 1] - 1
            else:
                return None, None
        else:
            return top, bottom


def approx_match(string1, string2, d):
    assert len(string1) == len(string2)
    mismatches = 0
    for symbol1, symbol2 in zip(string1, string2):
        if symbol1 != symbol2:
            mismatches += 1
            if mismatches > d:
                return False
    return True


def approx_matching_positions(text, suffix_array, first_occurrence, last_column, pattern, count, d=1):
    assert text[-1] == '$'
    # Break up 'pattern' into d+1 seeds (substrings)
    n = len(pattern)
    k = n // (d + 1)  # The length of each seed, except the last one
    seeds = [pattern[i * k:i * k + k] for i in range(0, d)]
    seeds.append(pattern[d * k:])
    seeds_pos = [i * k for i in range(0, d + 1)]  # The starting position of each seed withing 'pattern'
    assert sum([len(seed) for seed in seeds]) == len(pattern)

    all_matchin_pos = set()
    #  For every seed ...
    for seed, seed_pos in zip(seeds, seeds_pos):
        # ... try an exact match between the seed and the text
        top, bottom = matching_positions(first_occurrence, last_column, seed, count)
        if top is None:  # If there is no exact match, nothing to be done with the current seed
            continue
        # If there are matches of 'seed', determine all the corresponding positions in 'text' where 'pattern' should begin
        positions = [suffix_array[pos] - seed_pos for pos in range(top, bottom + 1)]
        ''' Only keep positions such that 'pattern' is entirely within 'text' and where 'pattern' has not been matched
        already  (if 'pattern' has at least one exact match in 'text', then different seeds will match to the same
        position) '''
        positions = [pos for pos in positions if pos >= 0 and pos + n < len(
            text) and pos not in all_matchin_pos]  # Keep in mind that 'text' ends with '$'
        ''' For every determined position, check if a prefix of 'text' starting at that position matches
        'pattern' with at most d differences '''
        matching_pos = {pos for pos in positions if approx_match(text[pos: pos + n], pattern, d)}
        all_matchin_pos |= matching_pos

    return list(all_matchin_pos)


def find_all(text, patterns, d=0):
    transformed_text = burrows_wheeler_transform(text)
    suffix_array = suffix_array_for_text(text)

    first_column = sorted(text)
    first_occurrence = {}
    for pos, symbol in enumerate(first_column):
        if first_occurrence.get(symbol) is None:
            first_occurrence[symbol] = pos

    symbols = set(transformed_text)
    count = {symbol: [0] * (len(transformed_text) + 1) for symbol in symbols}
    for pos, pos_symbol in enumerate(transformed_text):
        for symbol in symbols:
            count[symbol][pos + 1] = count[symbol][pos] + 1 if symbol == pos_symbol else count[symbol][pos]

    all_positions = []
    for pattern in patterns:
        if d == 0:
            top, bottom = matching_positions(first_occurrence, transformed_text, pattern, count)
            if top is not None:
                positions = [suffix_array[pos] for pos in range(top, bottom + 1)]
                all_positions.extend(positions)
        else:
            positions = approx_matching_positions(text, suffix_array, first_occurrence, transformed_text, pattern,
                                                  count, d)
            all_positions.extend(positions)
    return sorted(all_positions)
