import numpy as np
from scipy.spatial.distance import euclidean


def farthest_first_centers(k, data_points):
    points = np.array(data_points, dtype=float)
    centers = np.array([points[0]])
    while len(centers) < k:
        max_dist = float('-inf')
        max_dist_point = None
        for point in points:
            dist = min([euclidean(point, center) for center in centers])
            if dist > max_dist:
                max_dist = dist
                max_dist_point = point
        centers = np.append(centers, np.expand_dims(max_dist_point, axis=0), axis=0)
    return centers.tolist()


def sq_error_distortion(data_points, centers):
    points = np.array(data_points, dtype=float)
    centers = np.array(centers)
    dist = sum([min([euclidean(point, center) for center in centers]) ** 2 for point in points]) / len(points)
    return dist
