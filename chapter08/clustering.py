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


def lloyd_cluster(data_points, centers):
    points = np.array(data_points, dtype=float)
    centers = np.array(centers)
    k = len(centers)

    # As long as centers have moved in the previous iteration
    prev_centers = None
    while prev_centers is None or not np.all(np.isclose(prev_centers, centers, atol=0, rtol=1e-6)):
        clusters = {center_idx: [] for center_idx in range(0, k)}
        # Assign each data point to the closest center
        for point in points:
            closest_center_idx = np.argmin([euclidean(point, center) for center in centers])
            current_list = clusters[closest_center_idx]
            current_list.append(point)
            clusters[closest_center_idx] = current_list
            # clusters[closest_center_idx].append(point)
        # Re-assign each cluster center to the center of gravity for that cluster
        prev_centers = np.ndarray.copy(centers)
        for clusters_idx, cluster_points in clusters.items():
            cluster_points = np.array(cluster_points, dtype=float)
            center_of_gravity = np.mean(cluster_points, axis=0)
            centers[clusters_idx] = center_of_gravity

    return centers.tolist()
