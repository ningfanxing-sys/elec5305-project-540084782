import numpy as np
from numpy.polynomial.legendre import leggauss

from sh import sh_fit, sh_matrix
from sphere import fibonacci_layout, map_to_grid, to_xyz, voronoi_weights


def test_orthonormality():
    N = 8
    z, weight = leggauss(N + 1)
    az, el = np.meshgrid(np.arange(2 * N + 1) * 360 / (2 * N + 1),
                         np.rad2deg(np.arcsin(z)))
    Y = sh_matrix(N, az.ravel(), el.ravel())
    w = np.repeat(weight, 2 * N + 1) * 2 * np.pi / (2 * N + 1)
    np.testing.assert_allclose(Y.T @ (w[:, None] * Y), np.eye((N + 1) ** 2),
                               atol=1e-10, rtol=0)


def test_round_trip():
    rng = np.random.default_rng(0)
    az, el = fibonacci_layout(300)
    Y = sh_matrix(6, az, el)
    coef = rng.normal(size=(49, 2, 12)) + 1j * rng.normal(size=(49, 2, 12))
    h = np.einsum("pc,cef->pef", Y, coef)
    fitted = sh_fit(Y, h, voronoi_weights(to_xyz(az, el)))
    reconstructed = np.einsum("pc,cef->pef", Y, fitted)
    assert np.max(np.abs(reconstructed - h)) < 1e-8


def test_map_to_grid():
    xyz = to_xyz(*fibonacci_layout(80))
    subset = np.array([31, 4, 66, 12, 0])
    indices, error = map_to_grid(xyz[subset], xyz)
    np.testing.assert_array_equal(indices, subset)
    assert np.max(error) < 2e-6
    assert len(np.unique(indices)) == len(indices)
