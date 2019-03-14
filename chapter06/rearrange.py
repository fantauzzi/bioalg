from copy import deepcopy


def pretty_print(printme):
    """
    Prints the result of greedy_sorting() in a format accepted by the Stepik challenge.
    :param printme: The list of lists of numbers to be printed.
    """
    fmat = '{:+} ' * len(printme[0])
    fmat = fmat.rstrip(' ')
    for line in printme:
        print(fmat.format(*line))


def greedy_sorting(p):
    """
    Returns the sequence of permutations leading from a given permutation to the indentity permutation, implementing greedy sorting.
    :param p: The given permutation, a list of integer numbers optionally with sign; positions in the permutation must be numbered starting from 1 (not 0).
    :return: A list of permutations that, applying greedy sorting, lead from the given one to the identity one, a list of lists of integer numbers. The list does not contain the given permutation, but ends with the identity permutation.
    """
    steps = []
    for i in range(0, len(p)):
        if abs(p[i]) != i + 1:
            try:
                pos = p.index(i + 1, i + 1)
            except ValueError:
                pos = p.index(-(i + 1), i + 1)
            p = p[:i] + [-item for item in p[i:pos + 1][::-1]] + p[pos + 1:]
            steps.append(deepcopy(p))
            assert abs(p[i]) == i + 1
        if p[i] == - (i + 1):
            p[i] = i + 1
            steps.append(deepcopy(p))

    return steps


def count_breakpoints(p):
    """
    Returns the number of breakpoints in a given permutation.
    :param p: The permutation, a list of integer numbers.
    :return: The number of breakpoints, an integer.
    """
    p = [0]+p+[len(p)+1]
    count = sum([p[i+1]-p[i] != 1 for i in range(0, len(p)-1)])
    return count
