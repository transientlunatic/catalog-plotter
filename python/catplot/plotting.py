import os
import numpy as np
import yaml
import h5py
import matplotlib
from matplotlib.mlab import GaussianKDE
import matplotlib.patheffects as pe
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import h5py as h5
#from analysis import calculate_divergences

import matplotlib.ticker as ticker

from pyfonts import load_google_font, set_default_font

from .asimov import get_results_file
from . import plotting
from . import binning
from . import io
from .colour import colours
from scipy.stats import gaussian_kde
from typing import Optional, Sequence, Tuple, List

font = load_google_font("Source Code Pro")
set_default_font(font)

display_parameters = {
    "total_mass_source": { 
        "name": "Total mass (source frame)",
        "symbol": r"$M_{tot, src}$"
    },
    "mass_ratio": {
        "name": "Mass ratio",
        "symbol": r"$q$"
    }
}

ANALYSES = ["GWTC-1", "GWTC-2.1", "GWTC-3", "IAS-1", "IAS-2", "IAS-3", "IAS-4", "1-OGC", "2-OGC", "3-OGC", "4-OGC", "Beyond"]
analysis_hatching = {
    "IAS-1": "/",
    "Ares": "/",
    "IAS-2": "//",
    "IAS-4": "//",
    "IAS-3": "//",
    "1-OGC": "\\",
    "2-OGC": "\\\\",
    "4-OGC": "\\",
    "IAS-HM": "//",
}





class Bounded_2d_kde(gaussian_kde):
    r"""Two-dimensional Gaussian KDE that mirrors mass across optional bounds.

    This mirrors samples across the supplied rectangular boundaries so the KDE
    naturally accounts for soft boundaries by reflecting mass. Either lower or
    upper bounds may be ``None`` to indicate a single-sided boundary.
    """

    def __init__(self, pts, xlow=None, xhigh=None, ylow=None, yhigh=None, *args, **kwargs):
        assert pts.ndim == 2, 'Bounded_kde can only be two-dimensional'
        super(Bounded_2d_kde, self).__init__(pts, *args, **kwargs)
        self._xlow = xlow
        self._xhigh = xhigh
        self._ylow = ylow
        self._yhigh = yhigh

    @property
    def xlow(self):
        return self._xlow

    @property
    def xhigh(self):
        return self._xhigh

    @property
    def ylow(self):
        return self._ylow

    @property
    def yhigh(self):
        return self._yhigh

    def evaluate(self, points):
        points = np.atleast_2d(points)
        assert points.ndim == 2, 'points must be two-dimensional'
        x, y = points
        pdf = super(Bounded_2d_kde, self).evaluate(points)
        if self.xlow is not None:
            pdf += super(Bounded_2d_kde, self).evaluate([2*self.xlow - x, y])
        if self.xhigh is not None:
            pdf += super(Bounded_2d_kde, self).evaluate([2*self.xhigh - x, y])
        if self.ylow is not None:
            pdf += super(Bounded_2d_kde, self).evaluate([x, 2*self.ylow - y])
        if self.yhigh is not None:
            pdf += super(Bounded_2d_kde, self).evaluate([x, 2*self.yhigh - y])
        if self.xlow is not None:
            if self.ylow is not None:
                pdf += super(Bounded_2d_kde, self).evaluate([2*self.xlow - x, 2*self.ylow - y])
            if self.yhigh is not None:
                pdf += super(Bounded_2d_kde, self).evaluate([2*self.xlow - x, 2*self.yhigh - y])
        if self.xhigh is not None:
            if self.ylow is not None:
                pdf += super(Bounded_2d_kde, self).evaluate([2*self.xhigh - x, 2*self.ylow - y])
            if self.yhigh is not None:
                pdf += super(Bounded_2d_kde, self).evaluate([2*self.xhigh - x, 2*self.yhigh - y])
        return pdf

    def __call__(self, pts):
        pts = np.atleast_2d(pts)
        out_of_bounds = np.zeros(pts.shape[1], dtype='bool')
        if self.xlow is not None:
            out_of_bounds[pts[:, 0] < self.xlow] = True
        if self.xhigh is not None:
            out_of_bounds[pts[:, 0] > self.xhigh] = True
        if self.ylow is not None:
            out_of_bounds[pts[:, 1] < self.ylow] = True
        if self.yhigh is not None:
            out_of_bounds[pts[:, 1] > self.yhigh] = True
        results = self.evaluate(pts)
        results[out_of_bounds] = 0.
        return results


class Posterior2DPlotter:
    """Plotter for 2D posterior KDEs with credible-region contours and labels.

    The class computes a gridded kernel density estimate (KDE) from posterior
    samples, determines a density threshold that encloses a requested credible
    fraction, draws the corresponding contour, and attempts to position a small
    textual label that lies inside the contour.

    Attributes
    ----------
    x_pts : Optional[np.ndarray]
        Grid x-coordinates used for KDE evaluation (monotonic increasing).
    y_pts : Optional[np.ndarray]
        Grid y-coordinates used for KDE evaluation (monotonic increasing).
    z : Optional[np.ndarray]
        KDE values evaluated on the Cartesian product of ``x_pts`` and
        ``y_pts`` with shape (len(y_pts), len(x_pts)).
    delta_x : Optional[float]
        Width of the x-range (x_max - x_min) used when plotting.
    delta_y : Optional[float]
        Height of the y-range (y_max - y_min) used when plotting.

    Examples
    --------
    >>> plotter = Posterior2DPlotter()
    >>> ax = plotter.plot(metafile_p, ("mass_1", "mass_2"), label="Event A")
    """

    def __init__(self) -> None:
        self.x_pts: Optional[np.ndarray] = None
        self.y_pts: Optional[np.ndarray] = None
        self.z: Optional[np.ndarray] = None
        self.delta_x: Optional[float] = None
        self.delta_y: Optional[float] = None

    def compute_credible_threshold(self, z: np.ndarray, x_pts: Sequence[float], y_pts: Sequence[float], credible_fraction: float = 0.90, max_iter: int = 60) -> float:
        dx = x_pts[1] - x_pts[0] if len(x_pts) > 1 else 1.0
        dy = y_pts[1] - y_pts[0] if len(y_pts) > 1 else 1.0
        area_per_pixel = dx * dy
        flat = np.ravel(z)
        vals = flat[np.isfinite(flat)]
        if vals.size == 0:
            return float('nan')
        total_mass = float(np.sum(vals) * area_per_pixel)
        target = credible_fraction * total_mass
        if total_mass <= 0 or target <= 0:
            return float(np.nanpercentile(vals, credible_fraction * 100))
        def mass_above(t: float) -> float:
            return float(np.sum(vals[vals >= t]) * area_per_pixel)
    def bilinear_interp(self, xv: float, yv: float) -> float:
        """
        Perform bilinear interpolation on the KDE grid at the given (xv, yv) coordinates.

        Args:
            xv (float): The x-coordinate at which to interpolate.
            yv (float): The y-coordinate at which to interpolate.

        Returns:
            float: The interpolated value at (xv, yv), or NaN if out of bounds.
        """
        assert self.x_pts is not None and self.y_pts is not None and self.z is not None
        x_pts = self.x_pts
        y_pts = self.y_pts
        z = self.z
        if xv < x_pts[0] or xv > x_pts[-1] or yv < y_pts[0] or yv > y_pts[-1]:
            return float('nan')
        ix = np.searchsorted(x_pts, xv) - 1
        iy = np.searchsorted(y_pts, yv) - 1
        ix = int(np.clip(ix, 0, len(x_pts) - 2))
        iy = int(np.clip(iy, 0, len(y_pts) - 2))
        x1, x2 = x_pts[ix], x_pts[ix + 1]
        y1, y2 = y_pts[iy], y_pts[iy + 1]
        Q11 = z[iy, ix]
        Q21 = z[iy, ix + 1]
        Q12 = z[iy + 1, ix]
        Q22 = z[iy + 1, ix + 1]
        tx = (xv - x1) / (x2 - x1) if x2 != x1 else 0.0
        ty = (yv - y1) / (y2 - y1) if y2 != y1 else 0.0
        return float(Q11 * (1 - tx) * (1 - ty) + Q21 * tx * (1 - ty) + Q12 * (1 - tx) * ty + Q22 * tx * ty)
        ix = np.searchsorted(x_pts, xv) - 1
        iy = np.searchsorted(y_pts, yv) - 1
        ix = int(np.clip(ix, 0, len(x_pts) - 2))
        iy = int(np.clip(iy, 0, len(y_pts) - 2))
    def get_contour_segments(self, contour) -> List[np.ndarray]:
        """
        Extracts contour segments from a matplotlib contour object.

        Args:
            contour: A matplotlib.contour.QuadContourSet object.

        Returns:
            List[np.ndarray]: A list of numpy arrays, each representing a contour segment.
        """
        segs: List[np.ndarray] = []
        try:
            for lvl_segs in contour.allsegs:
                if isinstance(lvl_segs, list):
                    segs.extend(lvl_segs)
        except Exception:
            pass
        return segs

    def get_contour_segments(self, contour) -> List[np.ndarray]:
    def find_nearest_seg_point(self, segs: Sequence[np.ndarray], x0: float, y0: float) -> Tuple[Optional[np.ndarray], Optional[int]]:
        """
        Find the nearest point on contour segments to a given (x0, y0) coordinate.

        Args:
            segs (Sequence[np.ndarray]): List of contour segments, each as an array of points.
            x0 (float): X-coordinate of the reference point.
            y0 (float): Y-coordinate of the reference point.

        Returns:
            Tuple[Optional[np.ndarray], Optional[int]]: The segment containing the nearest point and the index of that point within the segment.
        """
        best_seg: Optional[np.ndarray] = None
        best_idx: Optional[int] = None
        best_d2 = np.inf
        for seg in segs:
            if len(seg) == 0:
                continue
            d2 = (seg[:, 0] - x0) ** 2 + (seg[:, 1] - y0) ** 2
            idx = int(np.argmin(d2))
            if d2[idx] < best_d2:
                best_d2 = d2[idx]
                best_seg = seg
                best_idx = idx
    def place_label_inside(self, txt, segs: Sequence[np.ndarray], thresh: float) -> bool:
        """Position a Matplotlib text object inside a credible-region contour.

        Parameters
        ----------
        txt : matplotlib.text.Text
            Text object to be positioned.
        segs : Sequence[np.ndarray]
            Contour segments as Nx2 arrays describing the contour vertices.
        thresh : float
            Density threshold defining the interior (points with density >= thresh
            are considered inside the credible region).

        Returns
        -------
        bool
            True if the label was moved (or its alpha adjusted) to lie within the
            credible region; False only if no placement succeeded.
        """
        assert self.delta_x is not None and self.delta_y is not None
        x0, y0 = txt.get_position()
            d2 = (seg[:, 0] - x0) ** 2 + (seg[:, 1] - y0) ** 2
            idx = int(np.argmin(d2))
            if d2[idx] < best_d2:
                best_d2 = d2[idx]
                best_seg = seg
                best_idx = idx
        return best_seg, best_idx

    def place_label_inside(self, txt, segs: Sequence[np.ndarray], thresh: float) -> bool:
        assert self.delta_x is not None and self.delta_y is not None
        x0, y0 = txt.get_position()
        seg, idx = self.find_nearest_seg_point(segs, x0, y0)
        offset_mag = 0.02 * max(self.delta_x, self.delta_y)
        def test_and_place(nx: float, ny: float) -> bool:
            zval = self.bilinear_interp(nx, ny)
            if np.isfinite(zval) and zval >= thresh:
                txt.set_position((nx, ny))
                return True
            return False
        if seg is None:
            ok = test_and_place(x0, y0 + 0.02 * self.delta_y)
            txt.set_alpha(0.6)
            return ok
        if idx <= 0:
            p1 = seg[0]
            p2 = seg[1] if seg.shape[0] > 1 else seg[0]
        elif idx >= seg.shape[0] - 1:
            p1 = seg[-2] if seg.shape[0] > 1 else seg[-1]
            p2 = seg[-1]
        else:
            p1 = seg[max(0, idx - 1)]
            p2 = seg[min(seg.shape[0] - 1, idx + 1)]
        tangent = p2 - p1
        normal = np.array([-tangent[1], tangent[0]], dtype=float)
        nrm = np.hypot(normal[0], normal[1])
        norm_unit = normal / nrm if nrm != 0 else np.array([0.0, 1.0])
        for direction in (1.0, -1.0):
            nx = x0 + norm_unit[0] * offset_mag * direction
            ny = y0 + norm_unit[1] * offset_mag * direction
            if test_and_place(nx, ny):
                txt.set_alpha(0.6)
                return True
        for scale in (0.5, 0.25, 0.1, 0.0):
            nx = x0 + norm_unit[0] * offset_mag * scale
            ny = y0 + norm_unit[1] * offset_mag * scale
            if test_and_place(nx, ny):
                txt.set_alpha(0.6)
                return True
        cont_point = seg[idx]
        centroid = seg.mean(axis=0)
        dir_to_centroid = centroid - cont_point
        
    def plot(self,
             metafile_p: str,
             parameters: Tuple[str, str],
             ax=None,
             label: Optional[str] = None,
             gridsize: Tuple[int, int] = (100, 100),
             downsample: int = 2000,
             transform=None,
             x_range: Optional[Tuple[float, float]] = None,
             color: str = "k",
             linewidth: float = 1.0,
             y_range: Optional[Tuple[float, float]] = None,
             credible_fraction: float = 0.90):
        """
        Plot the 2D posterior KDE and credible contour for the given parameters.

        Args:
            metafile_p (str): Path to the metafile containing posterior samples.
            parameters (Tuple[str, str]): Tuple of parameter names to plot on x and y axes.
            ax (matplotlib.axes.Axes, optional): Matplotlib axis to plot on. If None, a new axis is created.
            label (str, optional): Label to display inside the credible region contour.
            gridsize (Tuple[int, int], optional): Number of grid points for KDE evaluation (default (100, 100)).
            downsample (int, optional): Number of posterior samples to use (default 2000).
            transform (callable, optional): Function to transform the data before KDE.
            x_range (Tuple[float, float], optional): Range for the x-axis.
            color (str, optional): Color for the contour and label (default "k").
            linewidth (float, optional): Line width for the contour (default 1.0).
            y_range (Tuple[float, float], optional): Range for the y-axis.
            credible_fraction (float, optional): Credible fraction for the contour (default 0.90).

        Returns:
            matplotlib.axes.Axes: The axis with the plot (ax).
        """
        if ax is None:
            f, ax = plt.subplots(1, 1, figsize=(10, 5), dpi=300)
             label: Optional[str] = None,
             gridsize: Tuple[int, int] = (100, 100),
             downsample: int = 2000,
             transform=None,
             x_range: Optional[Tuple[float, float]] = None,
             color: str = "k",
             linewidth: float = 1.0,
             y_range: Optional[Tuple[float, float]] = None,
             credible_fraction: float = 0.90):
        if ax is None:
            f, ax = plt.subplots(1, 1, figsize=(10, 5), dpi=300)
        with h5.File(metafile_p) as metafile:
            analyses = list(metafile.keys())
            for reserved in ("history", "version"):
                if reserved in analyses:
                    analyses.remove(reserved)
            analysis = analyses[0]
            posterior = metafile[analysis]['posterior_samples']
            class Posterior2DPlotter:
                """Plotter for 2D posterior KDEs with credible-region contours and labels.

                Parameters
                ----------
                None
                    This class holds state for a single plotting operation; call :meth:`plot`
                    to compute and render the contour.

                Notes
                -----
                The implementation is intentionally minimal and designed to work with the
                project's existing HDF5 PESummary metafiles. It expects that the named
                parameters are present in the metafile's posterior samples.
                """

                def __init__(self) -> None:
                    self.x_pts: Optional[np.ndarray] = None
                    self.y_pts: Optional[np.ndarray] = None
                    self.z: Optional[np.ndarray] = None
                    self.delta_x: Optional[float] = None
                    self.delta_y: Optional[float] = None

                def compute_credible_threshold(self,
                                               z: np.ndarray,
                                               x_pts: Sequence[float],
                                               y_pts: Sequence[float],
                                               credible_fraction: float = 0.90,
                                               max_iter: int = 60) -> float:
                    """Compute a density threshold that encloses a credible fraction.

                    Parameters
                    ----------
                    z : np.ndarray
                        2D array of density values (shape len(y_pts) x len(x_pts)).
                    x_pts, y_pts : Sequence[float]
                        Grid coordinates for x and y (monotonic increasing).
                    credible_fraction : float, optional
                        Fraction of total mass to enclose (default 0.90).
                    max_iter : int, optional
                        Maximum iterations for the binary search (default 60).

                    Returns
                    -------
                    float
                        Density threshold value.
                    """
                    dx = x_pts[1] - x_pts[0] if len(x_pts) > 1 else 1.0
                    dy = y_pts[1] - y_pts[0] if len(y_pts) > 1 else 1.0
                    area_per_pixel = dx * dy
                    flat = np.ravel(z)
                    vals = flat[np.isfinite(flat)]
                    if vals.size == 0:
                        return float('nan')
                    total_mass = float(np.sum(vals) * area_per_pixel)
                    target = credible_fraction * total_mass
                    if total_mass <= 0 or target <= 0:
                        return float(np.nanpercentile(vals, credible_fraction * 100))
                    def mass_above(t: float) -> float:
                        return float(np.sum(vals[vals >= t]) * area_per_pixel)
                    lo = float(np.min(vals))
                    hi = float(np.max(vals))
                    if mass_above(lo) < target:
                        return float(np.nanpercentile(vals, credible_fraction * 100))
                    left, right = lo, hi
                    for _ in range(max_iter):
                        mid = 0.5 * (left + right)
                        if mass_above(mid) > target:
                            left = mid
                        else:
                            right = mid
                    return 0.5 * (left + right)

                def bilinear_interp(self, xv: float, yv: float) -> float:
                    """Bilinearly interpolate the stored grid ``z`` at point ``(xv, yv)``.

                    Returns ``nan`` for points outside the grid.
                    """
                    assert self.x_pts is not None and self.y_pts is not None and self.z is not None
                    x_pts = self.x_pts
                    y_pts = self.y_pts
                    z = self.z
                    if xv < x_pts[0] or xv > x_pts[-1] or yv < y_pts[0] or yv > y_pts[-1]:
                        return float('nan')
                    ix = np.searchsorted(x_pts, xv) - 1
                    iy = np.searchsorted(y_pts, yv) - 1
                    ix = int(np.clip(ix, 0, len(x_pts) - 2))
                    iy = int(np.clip(iy, 0, len(y_pts) - 2))
                    x1, x2 = x_pts[ix], x_pts[ix + 1]
                    y1, y2 = y_pts[iy], y_pts[iy + 1]
                    Q11 = z[iy, ix]
                    Q21 = z[iy, ix + 1]
                    Q12 = z[iy + 1, ix]
                    Q22 = z[iy + 1, ix + 1]
                    tx = (xv - x1) / (x2 - x1) if x2 != x1 else 0.0
                    ty = (yv - y1) / (y2 - y1) if y2 != y1 else 0.0
                    return float(Q11 * (1 - tx) * (1 - ty) + Q21 * tx * (1 - ty) + Q12 * (1 - tx) * ty + Q22 * tx * ty)

                def get_contour_segments(self, contour) -> List[np.ndarray]:
                    """Flatten and return the contour segments from a Matplotlib QuadContourSet.

                    Parameters
                    ----------
                    contour : matplotlib.contour.QuadContourSet
                        Contour object returned by :meth:`matplotlib.axes.Axes.contour`.

                    Returns
                    -------
                    list of ndarray
                        A list of Nx2 arrays containing contour vertices for each segment.
                    """
                    segs: List[np.ndarray] = []
                    try:
                        for lvl_segs in contour.allsegs:
                            if isinstance(lvl_segs, list):
                                segs.extend(lvl_segs)
                    except Exception:
                        pass
                    return segs

                def find_nearest_seg_point(self, segs: Sequence[np.ndarray], x0: float, y0: float) -> Tuple[Optional[np.ndarray], Optional[int]]:
                    """Return (segment, index) of nearest contour point to (x0, y0)."""
                    best_seg: Optional[np.ndarray] = None
                    best_idx: Optional[int] = None
                    best_d2 = np.inf
                    for seg in segs:
                        if len(seg) == 0:
                            continue
                        d2 = (seg[:, 0] - x0) ** 2 + (seg[:, 1] - y0) ** 2
                        idx = int(np.argmin(d2))
                        if d2[idx] < best_d2:
                            best_d2 = d2[idx]
                            best_seg = seg
                            best_idx = idx
                    return best_seg, best_idx

                def place_label_inside(self, txt, segs: Sequence[np.ndarray], thresh: float) -> bool:
                    """Attempt to move a Matplotlib Text ``txt`` so it lies inside the contour ``z >= thresh``.

                    Strategy
                    --------
                    - find nearest contour segment and normal vector
                    - try small offsets along both normal directions and progressively smaller steps
                    - fallback to nudging toward the segment centroid

                    Returns
                    -------
                    bool
                        True if label was moved (or alpha set); always sets a subtle alpha for readability.
                    """
                    assert self.delta_x is not None and self.delta_y is not None
                    x0, y0 = txt.get_position()
                    seg, idx = self.find_nearest_seg_point(segs, x0, y0)
                    offset_mag = 0.02 * max(self.delta_x, self.delta_y)
                    def test_and_place(nx: float, ny: float) -> bool:
                        zval = self.bilinear_interp(nx, ny)
                        if np.isfinite(zval) and zval >= thresh:
                            txt.set_position((nx, ny))
                            return True
                        return False
                    # No segment: try small upward move
                    if seg is None:
                        ok = test_and_place(x0, y0 + 0.02 * self.delta_y)
                        txt.set_alpha(0.6)
                        return ok
                    # Build local normal
                    if idx <= 0:
                        p1 = seg[0]
                        p2 = seg[1] if seg.shape[0] > 1 else seg[0]
                    elif idx >= seg.shape[0] - 1:
                        p1 = seg[-2] if seg.shape[0] > 1 else seg[-1]
                        p2 = seg[-1]
                    else:
                        p1 = seg[max(0, idx - 1)]
                        p2 = seg[min(seg.shape[0] - 1, idx + 1)]
                    tangent = p2 - p1
                    normal = np.array([-tangent[1], tangent[0]], dtype=float)
                    nrm = np.hypot(normal[0], normal[1])
                    norm_unit = normal / nrm if nrm != 0 else np.array([0.0, 1.0])
                    # Try both normal directions
                    for direction in (1.0, -1.0):
                        nx = x0 + norm_unit[0] * offset_mag * direction
                        ny = y0 + norm_unit[1] * offset_mag * direction
                        if test_and_place(nx, ny):
                            txt.set_alpha(0.6)
                            return True
                    # Smaller steps
                    for scale in (0.5, 0.25, 0.1, 0.0):
                        nx = x0 + norm_unit[0] * offset_mag * scale
                        ny = y0 + norm_unit[1] * offset_mag * scale
                        if test_and_place(nx, ny):
                            txt.set_alpha(0.6)
                            return True
                    # Fallback: nudge toward centroid of the nearest segment
                    cont_point = seg[idx]
                    centroid = seg.mean(axis=0)
                    dir_to_centroid = centroid - cont_point
                    nrmdir = np.hypot(dir_to_centroid[0], dir_to_centroid[1])
                    if nrmdir > 0:
                        unit = dir_to_centroid / nrmdir
                        posx = cont_point[0] + unit[0] * 0.01 * max(self.delta_x, self.delta_y)
                        posy = cont_point[1] + unit[1] * 0.01 * max(self.delta_x, self.delta_y)
                        txt.set_position((posx, posy))
                    else:
                        txt.set_position((cont_point[0], cont_point[1]))
                    txt.set_alpha(0.6)
                    return True

                def plot(self,
                         metafile_p: str,
                         parameters: Tuple[str, str],
                         ax=None,
                         label: Optional[str] = None,
                         gridsize: Tuple[int, int] = (100, 100),
                         downsample: int = 2000,
                         transform=None,
                         x_range: Optional[Tuple[float, float]] = None,
                         color: str = "k",
                         linewidth: float = 1.0,
                         y_range: Optional[Tuple[float, float]] = None,
                         credible_fraction: float = 0.90):
                    """Compute KDE and draw a credible contour with an interior label.

                    Parameters
                    ----------
                    metafile_p : str
                        Path to a PESummary HDF5 metafile containing posterior samples.
                    parameters : tuple of str
                        Pair of parameter names to plot (x_param, y_param).
                    ax : matplotlib.axes.Axes, optional
                        Axes to draw on. If ``None`` a new figure and axes are created.
                    label : str, optional
                        Text label to place inside the contour (for example an event name).
                    gridsize : tuple of int, optional
                        Number of grid points in x and y (default (100, 100)).
                    downsample : int, optional
                        Number of posterior samples to draw for KDE estimation (default 2000).
                    transform : callable, optional
                        Optional transform applied to samples before KDE evaluation.
                    x_range, y_range : tuple, optional
                        (min, max) range for x and y axes. If None, ranges are inferred from
                        the posterior samples.
                    color : str, optional
                        Contour color (default 'k').
                    linewidth : float, optional
                        Contour line width (default 1.0).
                    credible_fraction : float, optional
                        Fraction of posterior mass to enclose with the contour (default 0.90).

                    Returns
                    -------
                    matplotlib.axes.Axes
                        Axes object the contour was drawn on.
                    """
                    if ax is None:
                        f, ax = plt.subplots(1, 1, figsize=(10, 5), dpi=300)
                    with h5.File(metafile_p) as metafile:
                        analyses = list(metafile.keys())
                        for reserved in ("history", "version"):
                            if reserved in analyses:
                                analyses.remove(reserved)
                        analysis = analyses[0]
                        posterior = metafile[analysis]['posterior_samples']
                        x = posterior[parameters[0]]
                        y = posterior[parameters[1]]
                        try:
                            selection = np.random.choice(len(x), downsample, replace=False)
                        except ValueError:
                            selection = np.arange(len(x))
                        data = np.vstack([x[selection], y[selection]])
                        if x_range is None:
                            x_range = (float(np.nanmin(x)), float(np.nanmax(x)))
                        if y_range is None:
                            y_range = (float(np.nanmin(y)), float(np.nanmax(y)))
                        self.delta_x = float(x_range[1] - x_range[0])
                        self.delta_y = float(y_range[1] - y_range[0])
                        if transform is None:
                            transform = lambda x: x
                        self.x_pts = np.linspace(x_range[0], x_range[1], gridsize[0])
                        self.y_pts = np.linspace(y_range[0], y_range[1], gridsize[1])
                        xx, yy = np.meshgrid(self.x_pts, self.y_pts)
                        positions = np.vstack([xx.ravel(), yy.ravel()])
                        kernel = Bounded_2d_kde(transform(data))
                        self.z = kernel(transform(positions)).T
                        self.z = np.reshape(self.z, xx.shape)
                        thresh = self.compute_credible_threshold(self.z, self.x_pts, self.y_pts, credible_fraction=credible_fraction)
                        if np.isnan(thresh):
                            thresh = float(np.nanpercentile(self.z, credible_fraction * 100))
                        contour = ax.contour(self.x_pts, self.y_pts, self.z, levels=[thresh], colors=color, linewidths=linewidth)
                        if label:
                            fmt = {thresh: label}
                            labels = ax.clabel(contour, fmt=fmt, inline=False, fontsize=8, colors=color)
                            segs = self.get_contour_segments(contour)
                            for txt in labels:
                                try:
                                    self.place_label_inside(txt, segs, thresh)
                                except Exception:
                                    try:
                                        txt.set_alpha(0.6)
                                    except Exception:
                                        pass
                        ax.set_xlabel(parameters[0])
                        ax.set_ylabel(parameters[1])
                    return ax
                            path_effects=[pe.Stroke(linewidth=1.5, foreground='w'), pe.Normal()])
                except KeyError:
                    pass
                except np.linalg.LinAlgError:
                    pass
                except ValueError:
                    print(f"\tNo prior found for {property}")
            #ax.set(yticklabels=[], xticklabels=[])
            ax.yaxis.set_major_formatter(plt.NullFormatter())
    fig.set_layout_engine("compressed")
    return fig

def plane_plot(catfile):
    properties = ("mass_1_source", "mass_2_source")
    
    events = io.read_catfile(catfile)

    property_map = {
        "mass_1_source": "$m_1$",
        "mass_2_source": "$m_2$"
        }
    
    
    bins = [np.linspace(0.1, 10, 200), np.linspace(0.1, 3, 200)]

    files = [event['metafile'] for event in events]
    
    data = binning.combine_analyses(files, properties,
                                    bins)

    f, ax = plt.subplots(1,1, figsize=(5*1.6, 3), dpi=300)
    ax.imshow(data.T, origin="lower", cmap="Greys", extent=[bins[0][0], bins[0][-1],
                                                            bins[1][0], bins[1][-1]])
    i = 0
    for event in events:
        if event.get('highlight'):
            colour = colours[i]
            i+= 1
        else:
            colour = 'k'
        cs, labels = binning.contours(event['metafile'], event['name'], ("mass_1_source", "mass_2_source"), 0.1, bins, colour=colour)
        
    ax.set_xlabel(property_map[properties[0]])
    ax.set_ylabel(property_map[properties[1]])
