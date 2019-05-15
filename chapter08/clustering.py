import numpy as np
from scipy.spatial.distance import euclidean


def farthest_first_centers(k, data_points):
    """
    Returns cluster centers for given data points, based on farthest-first traversal clustering heuristic.
    :param k: The number of required cluster centers, a positive integer.
    :param data_points: The data points, a sequence of sequences of numbers.
    :return: The centers, a sequence of sequences of numbers.
    """
    points = np.array(data_points, dtype=float)
    centers = np.array([points[0]])  # Initialise the choice of centers with the first point in the list
    while len(centers) < k:
        # Select the point from data_points with the longest distance from the already selected center points
        max_dist = float('-inf')
        max_dist_point = None
        for point in points:
            ''' The distance from a given point to the selected center points is the minimum of the distance
            between the given point and any of the center points. '''
            dist = min([euclidean(point, center) for center in centers])
            if dist > max_dist:
                max_dist = dist
                max_dist_point = point
        centers = np.append(centers, np.expand_dims(max_dist_point, axis=0), axis=0)
    return centers.tolist()


def sq_error_distortion(data_points, centers):
    """
    Computes and returns the squared error distortion of given data points with respect to given centers.
    :param data_points: The data points, a sequence of sequences of numbers.
    :param centers:The centers, a sequence of sequences of numbers.
    :return: The squared error distortion, a number.
    """
    points = np.array(data_points, dtype=float)
    centers = np.array(centers)
    dist = sum([min([euclidean(point, center) for center in centers]) ** 2 for point in points]) / len(points)
    return dist


def lloyd_cluster(data_points, centers):
    """
    Clusters data points based on Lloyd algorithm, and returns the cluster centers.
    :param data_points: The data points, a sequence of sequences of numbers.
    :param centers: The center points from where to start iterating the clustering algorithms; they do not need to be
    included in data_points.
    :return: The obtained cluster centers, as many as there where starting points in centers, a sequence of sequence of
    numbers.
    """

    points = np.array(data_points, dtype=float)
    centers = np.array(centers)
    k = len(centers)  # Number of clusters

    prev_centers = None
    # As long as centers have moved significantly in the last iteration...
    while prev_centers is None or not np.all(np.isclose(prev_centers, centers, atol=0, rtol=1e-6)):
        # Associate each cluster (numbered based on its appearance in 'centers') with its data points
        clusters = {center_idx: [] for center_idx in range(0, k)}
        # Assign each data point to the closest center, and therefore to its cluster
        for point in points:
            closest_center_idx = np.argmin([euclidean(point, center) for center in centers])
            clusters[closest_center_idx].append(point)
        # Re-assign each cluster center to the center of gravity of data points belonging to the cluster
        prev_centers = np.ndarray.copy(centers)  # Keep a copy of the current centers, to see if they move after update
        for clusters_idx, cluster_points in clusters.items():
            cluster_points = np.array(cluster_points, dtype=float)
            center_of_gravity = np.mean(cluster_points, axis=0)
            centers[clusters_idx] = center_of_gravity

    return centers.tolist()


def dist(d, cluster1, cluster2):
    the_distance = sum([d[(row,)][(col,)] for row in cluster1 for col in cluster2]) / (
            len(cluster1) * len(cluster2))
    return the_distance


def hierarchical_clusters(d):
    def arg_min_dist(d, clusters):
        min_dist = float('inf')
        min_idx = None
        n = len(d)
        for row in clusters:
            for col in clusters:
                if d[row][col] < min_dist and row != col:
                    min_dist = d[row][col]
                    min_idx = row, col
        return min_idx

    n = len(d)
    d = {(row,): {(col,): d[row][col] for col in range(0, n)} for row in range(0, n)}

    clutering_steps = []  # Keep track of the succession of clustering operations, the return value
    clusters = set(d)  # Keep track of clusters, initialised with one cluster per item (leaf)
    while len(clusters) > 1:
        # Find the two closest clusters
        c_i, c_j = arg_min_dist(d, clusters)
        # Make a new cluster from their unionlumns corresponding to clusters c_i and c_j
        c_new = c_i + c_j
        # Update the set of clusters
        clusters = (clusters | {c_new}) - {c_i, c_j}
        # Update the list of performed steps (unions)
        clutering_steps.append(c_new)
        # Compute the new row (and column) to be inserted in d for cluster c_new
        new_row = {col: dist(d, c_new, col) for col in d.keys()}
        new_row[c_new] = .0
        # Insert it
        d[c_new] = new_row
        for cluster, row in d.items():
            row[c_new] = new_row[cluster]

    return clutering_steps


def max_distance(points, centers):
    max_dist = float('-inf')
    max_dist_point = None
    for point in points:
        dist = min([euclidean(point, center) for center in centers])
        if dist > max_dist:
            max_dist = dist
            max_dist_point = point
    return max_dist_point, max_dist
from math import exp

def hidden_matrix(data, centers):
    """h_m = [[1 / (euclidean(data[j], centers[i]) ** 2 )/ sum(
        [1 / (euclidean(data[j], centers[t])**2) for t in range(0, len(centers))]) for j in range(0, len(data))] for i in
           range(0, len(centers))]"""
    beta =1
    h_m = [[exp(-beta*euclidean(data[j], centers[i])) / sum(
        [exp(-beta*euclidean(data[j], centers[t])) for t in range(0, len(centers))]) for j in range(0, len(data))] for i
           in
           range(0, len(centers))]

    return h_m
