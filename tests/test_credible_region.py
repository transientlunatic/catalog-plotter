import numpy as np


def find_threshold(z, x_pts, y_pts, credible_fraction=0.90):
    """Find density threshold on a grid `z` that encloses `credible_fraction` of the mass.
    This mirrors the algorithm used in the notebook: compute area-per-pixel, total mass,
    then binary-search the threshold t so that mass(z >= t) ~= credible_fraction * total_mass.
    """
    dx = x_pts[1] - x_pts[0] if len(x_pts) > 1 else 1.0
    dy = y_pts[1] - y_pts[0] if len(y_pts) > 1 else 1.0
    area_per_pixel = dx * dy

    flat = np.ravel(z)
    vals = flat[np.isfinite(flat)]
    if vals.size == 0:
        return np.nan

    total_mass = np.sum(vals) * area_per_pixel
    if total_mass <= 0:
        # fallback to percentile
        return np.nanpercentile(vals, credible_fraction * 100)

    target = credible_fraction * total_mass

    def mass_above(t):
        return np.sum(vals[vals >= t]) * area_per_pixel

    lo = np.min(vals)
    hi = np.max(vals)

    # quick bounds check
    if mass_above(lo) < target:
        return np.nanpercentile(vals, credible_fraction * 100)

    left, right = lo, hi
    for _ in range(80):
        mid = 0.5 * (left + right)
        if mass_above(mid) > target:
            left = mid
        else:
            right = mid
    return 0.5 * (left + right)


def test_find_threshold_on_gaussian():
    # create a normalized-ish 2D Gaussian grid (not normalized to 1, but consistent)
    sigma = 1.0
    x_pts = np.linspace(-6, 6, 400)
    y_pts = np.linspace(-6, 6, 400)
    xx, yy = np.meshgrid(x_pts, y_pts)
    z = np.exp(-(xx ** 2 + yy ** 2) / (2 * sigma ** 2))

    # find threshold that should enclose ~90% of the mass on the grid
    thresh = find_threshold(z, x_pts, y_pts, credible_fraction=0.90)

    # compute fraction actually enclosed
    dx = x_pts[1] - x_pts[0]
    dy = y_pts[1] - y_pts[0]
    area = dx * dy
    flat = z.ravel()
    vals = flat[np.isfinite(flat)]
    total_mass = np.sum(vals) * area
    enclosed = np.sum(vals[vals >= thresh]) * area
    frac = enclosed / total_mass

    # allow small numerical tolerance
    assert abs(frac - 0.90) < 2e-3, f"enclosed fraction {frac:.6f} differs from 0.9"
