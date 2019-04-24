import matplotlib.pyplot as plt
import networkx.drawing.nx_pylab as nxp
import phylogeny
from phylogeny import add_node, BinTreeAdj
from pathlib import Path
import networkx as nx
import networkx.algorithms.isomorphism as iso


def pretty_print_matrix(matrix):
    """
    Prints the matrix with distances between nodes in a graph, as returned by dist_between_leaves().
    :param matrix: The matrix.
    """
    items = sorted(matrix)
    for row_item in items:
        line = [matrix[row_item][col_item] for col_item in items]
        print(*line, sep=' ')


def pretty_print_adjacencies(tree):
    for node in sorted(tree):
        for adj_node, adj_dist in tree[node]:
            print(node, '->', adj_node, ':', adj_dist, sep='')


def pretty_print_parsimony(tree, results, score):
    print(score)
    for node, adjs in tree.items():
        if adjs.left is not None:
            dist = phylogeny.hamming_distance(results[node], results[adjs.left])
            print(results[node], '->', results[adjs.left], ':', dist, sep='')
        if adjs.right is not None:
            dist = phylogeny.hamming_distance(results[node], results[adjs.right])
            print(results[node], '->', results[adjs.right], ':', dist, sep='')
        if adjs.parent is not None:
            dist = phylogeny.hamming_distance(results[node], results[adjs.parent])
            print(results[node], '->', results[adjs.parent], ':', dist, sep='')


def fetch_stepik_input(file_name):
    """
    Fetches the adjacency lists of a tree from a file, with input from the Stepik challenge "Distances Between Leaves
    Problem"
    :param file_name: The file name, with its relative path.
    :return: The adjacency lists, as expected by dist_between_leaves().
    """
    tree = {}
    with open(file_name) as input_file:
        input_file.readline().rstrip('\n')
        lines = input_file.readlines()
        for line in lines:
            line = line.rstrip('\n')
            left, right = line.split('->')
            node1 = int(left)
            right_node2, right_weight = right.split(':')
            node2 = int(right_node2)
            weight = int(right_weight)
            adjs = tree.get(node1, [])
            adjs.append((node2, weight))
            tree[node1] = adjs

    return tree


def parse_stepik_matrix(input_file):
    """
    Fetches and returns a matrix of distances between nodes in a graph from a given file.
    :param file_name: The file name, with its relative path. The content of the file must be formatted like the
    output of the stepik challenge "Distances Between Leaves Problem".
    :return: The distances matrix, a dictionary of dictionaries, same as returned by dist_between_leaves()
    """
    matrix = {}
    lines = input_file.readlines()

    for row_i, line in enumerate(lines):
        matrix[row_i] = {}
        line = line.rstrip('\n')
        for col_i, item in enumerate(line.split(' ')):
            matrix[row_i][col_i] = int(item)

    return matrix


def fetch_stepik_result(file_name):
    with open(file_name) as input_file:
        matrix = parse_stepik_matrix(input_file)
    return matrix


def fetch_stepik_limb_length_input(file_name):
    with open(file_name) as input_file:
        n = input_file.readline().rstrip('\n')
        n = int(n)
        j = input_file.readline().rstrip('\n')
        j = int(j)
        matrix = parse_stepik_matrix(input_file)
    return n, j, matrix


def fetch_stepik_additive_phylogeny_input(file_name):
    with open(file_name) as input_file:
        n = input_file.readline().rstrip('\n')
        n = int(n)
        matrix = parse_stepik_matrix(input_file)
    return n, matrix


def fetch_small_parsimony_input(file_name):
    tree = nx.Graph()
    with open(file_name) as input_file:
        n_leaf = input_file.readline().rstrip('\n')
        n_leaf = int(n_leaf)
        for i_leaf in range(0, n_leaf):
            line = input_file.readline().rstrip('\n')
            left, right = line.rstrip('\n').split('->')
            parent = int(left)
            tree.add_edge(parent, i_leaf)
            tree.nodes[i_leaf]['label'] = right
            tree.nodes[i_leaf]['parent'] = parent
        lines = input_file.readlines()
        for line in lines:
            left, right = line.rstrip('\n').split('->')
            node1, node2 = int(left), int(right)
            tree.add_edge(node1, node2)
            tree.nodes[node2]['parent'] = node1
        count = 0
        for node in tree.nodes():
            if tree.nodes[node].get('parent') is None:
                count += 1
                assert count == 1
                tree.nodes[node]['parent'] = None
    return tree


def fetch_small_parsimony_unrooted(file_name):
    tree = nx.Graph()
    with open(file_name) as input_file:
        n_leaf = input_file.readline().rstrip('\n')
        n_leaf = int(n_leaf)
        for i_leaf in range(0, n_leaf):
            input_file.readline()
            line = input_file.readline().rstrip('\n')
            left, right = line.rstrip('\n').split('->')
            parent = int(left)
            tree.add_edge(parent, i_leaf)
            tree.nodes[i_leaf]['label'] = right
        lines = input_file.readlines()
        for line in lines:
            left, right = line.rstrip('\n').split('->')
            node1, node2 = int(left), int(right)
            tree.add_edge(node1, node2)
    return tree


def pretty_print_small_parsimony(tree, score):
    print(score)
    for node1, node2 in tree.edges():
        label1 = tree.nodes[node1]['label']
        label2 = tree.nodes[node2]['label']
        dist = phylogeny.hamming_distance(label1, label2)
        print(label1, '->', label2, ':', dist, sep='')
        print(label2, '->', label1, ':', dist, sep='')


def fetch_nn_input(file_name):
    tree = nx.Graph()
    with open(file_name) as input_file:
        line = input_file.readline().rstrip('\n')
        a, b = line.split(' ')
        a, b = int(a), int(b)
        lines = input_file.readlines()
        for line in lines:
            node1, node2 = line.split('->')
            node1, node2 = int(node1), int(node2)
            tree.add_edge(node1, node2)

    return a, b, tree


def pretty_print_nn_trees(tree1, tree2):
    for tree in (tree1, tree2):
        for n1, n2 in tree.edges():
            print(n1, '->', n2, sep='')
            print(n2, '->', n1, sep='')
        print()
