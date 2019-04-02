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
