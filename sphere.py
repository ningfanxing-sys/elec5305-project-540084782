import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial import SphericalVoronoi

# Lowest elevation used in both databases: SADIE II onsets are implausible below this (check_data.py)
EL_MIN = -65


def to_xyz(az, el):
    az, el = np.deg2rad(az), np.deg2rad(el)
    return np.column_stack((np.cos(el) * np.cos(az),
                            np.cos(el) * np.sin(az), np.sin(el)))


def angles(a, b):
    return np.rad2deg(np.arccos(np.clip(a @ b.T, -1, 1)))


def fibonacci_layout(M, el_min=-90):
    z_min = np.sin(np.deg2rad(el_min))
    z = z_min + (np.arange(M) + 0.5) / M * (1 - z_min)
    az = np.arange(M) * (180 * (3 - np.sqrt(5))) % 360
    return az, np.rad2deg(np.arcsin(z))


def random_layout(M, seed, el_min=-90):
    rng = np.random.default_rng(seed)
    az = rng.uniform(0, 360, M)
    z = rng.uniform(np.sin(np.deg2rad(el_min)), 1, M)
    return az, np.rad2deg(np.arcsin(z))


def nominal_layouts(M):
    return {"uniform": fibonacci_layout(M, EL_MIN), "random1": random_layout(M, 1, EL_MIN),
            "random2": random_layout(M, 2, EL_MIN), "gap": fibonacci_layout(M, el_min=-30)}


def map_to_grid(nominal_xyz, grid_xyz):
    distance = angles(nominal_xyz, grid_xyz)
    rows, indices = linear_sum_assignment(distance)
    return indices, distance[rows, indices]


def voronoi_weights(xyz):
    area = SphericalVoronoi(xyz).calculate_areas()
    return area * (4 * np.pi / area.sum())
