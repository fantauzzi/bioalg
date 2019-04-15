from collections import deque, namedtuple
from alignment import visit_trie_level_order


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


def pretty_print_trie_adj_lists2(root):
    to_be_visited = deque([root])
    while to_be_visited:
        current_node = to_be_visited.popleft()
        print(current_node.data, ' -> ', end='')
        for label, child_node in current_node.symbol_to_child.items():
            print('(child={}, label={}, weight={}) '.format(child_node.data, label, current_node.weights[label]),
                  end='')
            to_be_visited.append(child_node)
        print('(parent={})'.format(current_node.parent.data if current_node.parent is not None else 'None'))


def pretty_print_trie_adj_lists(root):
    nodes = visit_trie_level_order(root)
    for current_node in nodes:
        print(current_node.data, ' -> ', end='')
        for label, child_node in current_node.symbol_to_child.items():
            print('(child={}, label={}, weight={}) '.format(child_node.data, label, current_node.weights[label]),
                  end='')
        print('(parent={})'.format(current_node.parent.data if current_node.parent is not None else 'None'))


def fetch_string(file_name):
    with open(file_name) as input_file:
        text = input_file.readline().rstrip('\n')
    return text


def fetch_sequence_of_int(file_name):
    with open(file_name) as input_file:
        line = input_file.readline().rstrip('\n')

    separator = ', ' if line.find(',') >=0 else ' '
    seq = [int(item) for item in line.split(separator)]
    return seq


def fetch_BW_matching_input(file_name):
    with open(file_name) as input_file:
        text = input_file.readline().rstrip('\n')
        text2 = input_file.readline().rstrip('\n')
    pattern = text2.split(' ')
    return text, pattern


ChildInfo = namedtuple('ChildInfo', ['data', 'symbol', 'weight', 'position', 'length'])
NodeInfo = namedtuple('NodeInfo', ['data', 'parent_data', 'label', 'children'])


def fetch_find_all_input(file_name):
    with open(file_name) as input_file:
        text = input_file.readline().rstrip('\n')
        patterns = input_file.readlines()

    patterns = [item.rstrip('\n') for item in patterns]
    return text, patterns



def serialise_suffix_tree(root):
    nodes = visit_trie_level_order(root)
    serialised = []
    for current_node in nodes:
        serialised_node = []
        for label, child_node in current_node.symbol_to_child.items():
            child_info = ChildInfo(data=child_node.data,
                                   symbol=label,
                                   weight=current_node.weights.get(label),
                                   position=current_node.position.get(label),
                                   length=current_node.length.get(label))
            serialised_node.append(child_info)
        node_info = NodeInfo(data=current_node.data,
                             parent_data=current_node.parent.data if current_node.parent is not None else None,
                             label=current_node.label,
                             children=sorted(serialised_node))
        serialised.append(node_info)
    return sorted(serialised)
