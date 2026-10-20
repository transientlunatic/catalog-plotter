"""
Small multiples of marginal posteriors, one row per event.

This module is deliberately free of network access and of the heavier imports
in :mod:`catplot.plotting`, so it can be used (and tested) on its own.

The workflow is in two stages:

1. :func:`summarise_events` reads each event's posterior (and prior) samples
   once, and reduces them to densities on a plotting grid, a
   Jensen--Shannon divergence for each parameter, and a network SNR.
2. :func:`plot_small_multiples` and :func:`plot_divergence_snr` draw
   figures from those summaries.
"""

from dataclasses import dataclass, field

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d

from .io import read_metafile

# Parameter definitions.  ``xrange`` is the plotted range, ``bounds`` are
# hard physical/prior bounds used to reflect the KDE (None = unbounded),
# ``log`` means the density is estimated and plotted in ln(x).
PARAMETERS = {
    "chirp_mass_source": {
        "symbol": r"$\mathcal{M}_\mathrm{src}\,/\,M_\odot$",
        "xrange": (2, 100),
        "bounds": (None, None),
        "log": True,
        "plain_ticks": True,
        "ticks": [3, 10, 30, 100],
    },
    "mass_ratio": {
        "symbol": r"$q$",
        "xrange": (0.05, 1),
        "bounds": (None, 1),
        "ticks": [0.5, 1],
    },
    "chi_eff": {
        "symbol": r"$\chi_\mathrm{eff}$",
        "xrange": (-1, 1),
        "bounds": (-1, 1),
        "ticks": [-1, 0, 1],
    },
    "chi_p": {
        "symbol": r"$\chi_\mathrm{p}$",
        "xrange": (0, 1),
        "bounds": (0, 1),
        "ticks": [0, 0.5, 1],
    },
    "luminosity_distance": {
        "symbol": r"$d_L\,/\,\mathrm{Mpc}$",
        "xrange": (100, 20000),
        "bounds": (None, None),
        "log": True,
        "ticks": [100, 1000, 10000],
    },
    "chirp_mass": {
        "symbol": r"$\mathcal{M}\,/\,M_\odot$",
        "xrange": (0, 150),
        "bounds": (0, None),
        "ticks": [0, 50, 100, 150],
    },
    "theta_jn": {
        "symbol": r"$\theta_{JN}$",
        "xrange": (0, np.pi),
        "bounds": (0, np.pi),
        "ticks": [0, np.pi / 2, np.pi],
    },
}

DEFAULT_PROPERTIES = [
    "chirp_mass_source",
    "mass_ratio",
    "chi_eff",
    "chi_p",
    "luminosity_distance",
]

SNR_KEYS = ("network_matched_filter_snr", "network_optimal_snr")

_PAPER_STYLE = {
    "font.family": "serif",
    "font.serif": ["STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.linewidth": 0.6,
}

_MAX_SAMPLES = 20000
JSD_COLOUR = "#c8561e"


def _transform(samples, spec):
    """Map samples into the coordinate in which the density is estimated."""
    return np.log(samples) if spec.get("log") else np.asarray(samples, dtype=float)


def bounded_kde(samples, grid, bounds=(None, None)):
    """
    Gaussian KDE of ``samples`` on ``grid``, reflected at any finite bounds.

    Reflection stops the estimate from leaking past (and being depleted
    near) hard edges such as q = 1 or chi_p = 0.  The bandwidth is Scott's
    rule computed from the original, un-reflected samples.

    The returned density is normalised so that it integrates to one over the
    allowed region.
    """
    samples = np.asarray(samples, dtype=float)
    samples = samples[np.isfinite(samples)]
    rng = np.random.default_rng(0)
    if samples.size > _MAX_SAMPLES:
        samples = rng.choice(samples, _MAX_SAMPLES, replace=False)
    n = samples.size
    bandwidth = samples.std() * n ** (-1 / 5)
    if bandwidth == 0 or n < 2:
        raise ValueError("Samples have no spread.")

    lo, hi = bounds
    copies = [samples]
    if lo is not None:
        copies.append(2 * lo - samples)
    if hi is not None:
        copies.append(2 * hi - samples)
    points = np.concatenate(copies)

    # Bin the (reflected) samples on a fine grid and smooth with a Gaussian
    # of the KDE bandwidth, which is equivalent to a KDE but much faster.
    grid = np.asarray(grid, dtype=float)
    low = min(grid.min(), points.min()) - 4 * bandwidth
    high = max(grid.max(), points.max()) + 4 * bandwidth
    nbins = int(min(max((high - low) / bandwidth * 8, 512), 20000))
    counts, edges = np.histogram(points, bins=nbins, range=(low, high))
    width = edges[1] - edges[0]
    smooth = gaussian_filter1d(counts.astype(float), bandwidth / width, mode="constant")
    centres = 0.5 * (edges[1:] + edges[:-1])
    density = np.interp(grid, centres, smooth) / (n * width)

    # Enforce the support
    if lo is not None:
        density[grid < lo] = 0
    if hi is not None:
        density[grid > hi] = 0
    return density


def jensen_shannon_divergence(p, q):
    """
    Jensen--Shannon divergence, in bits (so between 0 and 1), between two
    densities sampled on the same uniform grid.
    """
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    p = p / p.sum()
    q = q / q.sum()
    m = 0.5 * (p + q)

    def kl(a, b):
        mask = a > 0
        return np.sum(a[mask] * np.log2(a[mask] / b[mask]))

    return 0.5 * kl(p, m) + 0.5 * kl(q, m)


@dataclass
class EventSummary:
    """Plot-ready reduction of one event's posterior and prior samples."""

    name: str
    snr: float = np.nan
    grid: dict = field(default_factory=dict)
    posterior: dict = field(default_factory=dict)
    prior: dict = field(default_factory=dict)
    jsd: dict = field(default_factory=dict)


def _summarise_property(post, prior, spec, ngrid=400):
    """Return (display_grid, post_density, prior_density, jsd) for a property."""
    bounds = spec.get("bounds", (None, None))
    xlo, xhi = spec["xrange"]
    tpost = _transform(post, spec)
    tprior = _transform(prior, spec) if prior is not None else None

    if spec.get("log"):
        display = np.linspace(np.log(xlo), np.log(xhi), ngrid)
    else:
        display = np.linspace(xlo, xhi, ngrid)
    p = bounded_kde(tpost, display, bounds)
    q = bounded_kde(tprior, display, bounds) if tprior is not None else None

    jsd = np.nan
    if tprior is not None:
        # The divergence is evaluated over the full extent of the samples,
        # not just the plotted window, so it isn't altered by the axis range.
        both = np.concatenate([tpost, tprior])
        lo = both.min() if bounds[0] is None else bounds[0]
        hi = both.max() if bounds[1] is None else bounds[1]
        pad = 0.05 * (hi - lo)
        full = np.linspace(
            lo - (pad if bounds[0] is None else 0),
            hi + (pad if bounds[1] is None else 0),
            1024,
        )
        jsd = jensen_shannon_divergence(
            bounded_kde(tpost, full, bounds), bounded_kde(tprior, full, bounds)
        )

    axis = np.exp(display) if spec.get("log") else display
    return axis, p, q, jsd


def summarise_events(events, properties=DEFAULT_PROPERTIES, analysis=None, fs=None):
    """
    Read posterior and prior samples for each event and reduce them to
    densities, divergences and a network SNR.

    Parameters
    ----------
    events : list of dict
        Each needs ``name`` and ``metafile`` keys, as returned by
        :func:`catplot.asimov.get_approved_event_data`.
    properties : list of str
        Keys of :data:`PARAMETERS`.

    Returns
    -------
    list of EventSummary
        Events whose files can't be read are skipped with a warning.
    """
    summaries = []
    for event in events:
        summary = EventSummary(name=event["name"])
        try:
            with read_metafile(event["metafile"], analysis=analysis, fs=fs) as meta:
                posterior = meta["posterior_samples"]
                names = posterior.dtype.names
                for key in SNR_KEYS:
                    if key in names:
                        summary.snr = float(np.median(np.real(posterior[key])))
                        break
                priors = meta["priors"]["samples"] if "priors" in meta else {}
                for prop in properties:
                    spec = PARAMETERS[prop]
                    prior = np.array(priors[prop]) if prop in priors else None
                    try:
                        axis, p, q, jsd = _summarise_property(
                            np.array(posterior[prop]), prior, spec
                        )
                    except ValueError:
                        continue
                    summary.grid[prop] = axis
                    summary.posterior[prop] = p
                    summary.prior[prop] = q
                    summary.jsd[prop] = jsd
        except (FileNotFoundError, OSError, KeyError) as error:
            print(f"Skipping {event['name']}: {error}")
            continue
        summaries.append(summary)
    return summaries


def sort_summaries(summaries, sort_by="snr"):
    """Order summaries by ``snr`` (highest first), ``name``, or a property's ``jsd``."""
    if sort_by == "name":
        return sorted(summaries, key=lambda s: s.name)
    if sort_by == "snr":
        return sorted(summaries, key=lambda s: -np.nan_to_num(s.snr, nan=-np.inf))
    raise ValueError(f"Unknown sort key {sort_by!r}")


def split_pages(items, rows_per_page):
    """Split ``items`` into the fewest pages of at most ``rows_per_page``, as evenly as possible."""
    if not rows_per_page or rows_per_page >= len(items):
        return [list(items)]
    npages = int(np.ceil(len(items) / rows_per_page))
    return [[items[i] for i in chunk] for chunk in np.array_split(np.arange(len(items)), npages)]


def _tick_label(prop, tick):
    if PARAMETERS[prop].get("log") and not PARAMETERS[prop].get("plain_ticks"):
        return rf"$10^{{{int(round(np.log10(tick)))}}}$"
    if prop == "theta_jn":
        return {0: "0", np.pi / 2: r"$\pi/2$", np.pi: r"$\pi$"}[tick]
    return f"{tick:g}"


def plot_small_multiples(
    summaries,
    properties=DEFAULT_PROPERTIES,
    colour="#2b5c8a",
    jsd=False,
    jsd_max=1.0,
    show_snr=True,
    width=7.0,
    row_height=0.28,
    fontsize=6.5,
    style="paper",
):
    """
    Draw one row per event and one column per parameter.

    Posterior (filled) and prior (line) share one scale in each cell, with the
    taller of the two reaching the top, so a posterior that is more peaked
    than the prior stands above it.

    Parameters
    ----------
    summaries : list of EventSummary
        In the order that they are to appear, top to bottom.
    colour : str
        A single colour for every row, or ``"cycle"`` to colour each row
        differently.
    jsd : bool
        If true, draw a bar below each cell whose length is the
        Jensen--Shannon divergence of that posterior from its prior, with the
        full cell width corresponding to ``jsd_max`` bits.
    show_snr : bool
        Print the median network matched-filter SNR beside the event name.
    """
    from .colour import colours as palette

    nrows = len(summaries)
    ncols = len(properties)
    label_cols = 2 if show_snr else 1
    header = 0.35
    footer = 0.4
    height = nrows * row_height + header + footer

    with plt.rc_context(_PAPER_STYLE if style == "paper" else {}):
        fig = plt.figure(figsize=(width, height))
        ratios = [1.5] + ([0.45] if show_snr else []) + [2.0] * ncols
        gs = fig.add_gridspec(
            nrows=nrows,
            ncols=len(ratios),
            width_ratios=ratios,
            left=0.01,
            right=0.985,
            top=1 - header / height,
            bottom=footer / height,
            hspace=0.0,
            wspace=0.14,
        )

        for row, summary in enumerate(summaries):
            row_colour = palette[row % len(palette)] if colour == "cycle" else colour

            ax = fig.add_subplot(gs[row, 0])
            ax.axis("off")
            ax.text(0, 0.5, summary.name, va="center", ha="left", fontsize=fontsize,
                    family="monospace")
            if show_snr:
                ax = fig.add_subplot(gs[row, 1])
                ax.axis("off")
                label = "" if np.isnan(summary.snr) else f"{summary.snr:.1f}"
                ax.text(1, 0.5, label, va="center", ha="right", fontsize=fontsize,
                        color="0.35")

            for col, prop in enumerate(properties):
                ax = fig.add_subplot(gs[row, label_cols + col])
                spec = PARAMETERS[prop]
                log = spec.get("log", False)
                if log:
                    ax.set_xscale("log")
                ax.set_xlim(*spec["xrange"])
                ax.set_ylim(-0.2 if jsd else -0.04, 1.08)
                ax.set_yticks([])
                for side in ("left", "right", "top"):
                    ax.spines[side].set_visible(False)
                ax.spines["bottom"].set_color("0.7")
                ax.set_xticks(spec["ticks"])
                ax.grid(axis="x", color="0.9", linewidth=0.5)
                ax.set_axisbelow(True)
                ax.tick_params(axis="x", length=0, labelbottom=False)
                ax.minorticks_off()

                if prop in summary.posterior:
                    x = summary.grid[prop]
                    post = summary.posterior[prop]
                    prior = summary.prior[prop]
                    scale = max(post.max(), prior.max() if prior is not None else 0)
                    ax.fill_between(x, post / scale, color=row_colour, alpha=0.55,
                                    linewidth=0)
                    ax.plot(x, post / scale, color=row_colour, linewidth=0.6)
                    if prior is not None:
                        ax.plot(x, prior / scale, color="0.25", linewidth=0.6,
                                linestyle=(0, (3, 1.5)))
                    if jsd and np.isfinite(summary.jsd[prop]):
                        # xmin/xmax are axes fractions, so the log axis doesn't distort the bar.
                        frac = min(summary.jsd[prop] / jsd_max, 1.0)
                        ax.axhspan(-0.17, -0.09, xmin=0, xmax=1, color="0.93",
                                   linewidth=0)
                        ax.axhspan(-0.17, -0.09, xmin=0, xmax=frac, color=JSD_COLOUR,
                                   linewidth=0)

                if row == 0:
                    ax.set_title(spec["symbol"], fontsize=fontsize + 2.5, pad=3)
                if row == nrows - 1:
                    ax.tick_params(axis="x", length=2.5, width=0.6, labelbottom=True,
                                   labelsize=fontsize, pad=1.5)
                    ax.set_xticklabels([_tick_label(prop, t) for t in spec["ticks"]])

        if jsd:
            _jsd_key(fig, header / height, fontsize, jsd_max)
    return fig


def _jsd_key(fig, header_frac, fontsize, jsd_max):
    """A small scale in the top-left corner explaining the divergence bars."""
    key = fig.add_axes([0.012, 1 - header_frac * 0.62, 0.105, header_frac * 0.3])
    key.set_xlim(0, 1)
    key.set_ylim(0, 1)
    key.axis("off")
    key.axhspan(0.15, 0.55, xmin=0, xmax=1, color="0.93", linewidth=0, clip_on=False)
    key.axhspan(0.15, 0.55, xmin=0, xmax=0.5, color=JSD_COLOUR, linewidth=0,
                clip_on=False)
    key.text(0, -0.25, "0", fontsize=fontsize - 0.5, ha="left", va="top")
    key.text(1, -0.25, f"{jsd_max:g} bit", fontsize=fontsize - 0.5, ha="right", va="top")
    key.text(0, 1.0, r"$D_\mathrm{JS}$(posterior, prior)", fontsize=fontsize - 0.5,
             ha="left", va="bottom", color=JSD_COLOUR)


def plot_divergence_snr(
    summaries,
    properties=DEFAULT_PROPERTIES,
    colour="#2b5c8a",
    width=7.0,
    fontsize=7,
    label_top=2,
    style="paper",
):
    """
    Jensen--Shannon divergence between posterior and prior against network
    SNR, one panel per parameter.  The ``label_top`` events with the largest
    divergence in each panel are labelled.
    """
    ncols = len(properties)
    with plt.rc_context(_PAPER_STYLE if style == "paper" else {}):
        fig, axes = plt.subplots(1, ncols, figsize=(width, 1.9), sharey=True)
        axes = np.atleast_1d(axes)
        for ax, prop in zip(axes, properties):
            pts = [(s.snr, s.jsd[prop], s.name) for s in summaries
                   if prop in s.jsd and np.isfinite(s.snr) and np.isfinite(s.jsd[prop])]
            snr, div, names = (np.array(v) for v in zip(*pts))
            ax.scatter(snr, div, s=12, color=colour, alpha=0.75, linewidths=0)
            ax.set_xscale("log")
            ax.set_xticks([6, 10, 20, 30])
            ax.set_xticklabels(["6", "10", "20", "30"])
            ax.minorticks_off()
            ax.set_xlim(5.5, 38)
            for i in np.argsort(div)[::-1][:label_top]:
                right = snr[i] > 15
                ax.annotate(names[i].replace("GW", ""), (snr[i], div[i]), fontsize=5,
                            xytext=(-3 if right else 3, -6 if right else 2),
                            ha="right" if right else "left",
                            textcoords="offset points", color="0.3")
            ax.set_title(PARAMETERS[prop]["symbol"], fontsize=fontsize + 2)
            ax.set_xlabel(r"Network SNR $\rho$", fontsize=fontsize)
            ax.tick_params(labelsize=fontsize - 0.5, width=0.6)
            ax.grid(color="0.92", linewidth=0.5)
            ax.set_axisbelow(True)
            for side in ("right", "top"):
                ax.spines[side].set_visible(False)
        axes[0].set_ylabel(r"$D_\mathrm{JS}$ (posterior, prior) / bits", fontsize=fontsize)
        fig.tight_layout(w_pad=0.6)
    return fig
