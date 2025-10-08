import os
from collections.abc import Sequence
from typing import Optional, Tuple


import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from matplotlib.mlab import GaussianKDE

#from analysis import calculate_divergences
from pyfonts import load_google_font, set_default_font
from scipy.stats import gaussian_kde

from .asimov import get_results_file

from .io import read_metafile

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

class Bounded_2d_kde(gaussian_kde):
    r"""Two-dimensional Gaussian KDE that mirrors mass across optional bounds.

    This mirrors samples across the supplied rectangular boundaries so the KDE
    naturally accounts for soft boundaries by reflecting mass. Either lower or
    upper bounds may be ``None`` to indicate a single-sided boundary.
    """

    def __init__(self, pts, xlow=None, xhigh=None, ylow=None, yhigh=None, *args, **kwargs):
        assert pts.ndim == 2, 'Bounded_kde can only be two-dimensional'
        super().__init__(pts, *args, **kwargs)
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
        pdf = super().evaluate(points)
        if self.xlow is not None:
            pdf += super().evaluate([2*self.xlow - x, y])
        if self.xhigh is not None:
            pdf += super().evaluate([2*self.xhigh - x, y])
        if self.ylow is not None:
            pdf += super().evaluate([x, 2*self.ylow - y])
        if self.yhigh is not None:
            pdf += super().evaluate([x, 2*self.yhigh - y])
        if self.xlow is not None:
            if self.ylow is not None:
                pdf += super().evaluate([2*self.xlow - x, 2*self.ylow - y])
            if self.yhigh is not None:
                pdf += super().evaluate([2*self.xlow - x, 2*self.yhigh - y])
        if self.xhigh is not None:
            if self.ylow is not None:
                pdf += super().evaluate([2*self.xhigh - x, 2*self.ylow - y])
            if self.yhigh is not None:
                pdf += super().evaluate([2*self.xhigh - x, 2*self.yhigh - y])
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


class Posterior2D:

    @classmethod
    def plot(cls,
             metafile: str,
             parameters: Sequence[str],
             xrange: Optional[Tuple[float, float]] = None,
             yrange: Optional[Tuple[float, float]] = None,
             color: str = "k",
             confidence_fraction: float = 0.9,
             gridsize: Tuple[int, int] = (100, 100),
             transform=None,
             label=None,
             contour_type: str = 'line',
             ax=None):
        """
        Plot a 2D posterior from a metafile.

        Parameters
        ----------
        metafile : str
            Path to the HDF5 metafile.
        parameters : Sequence[str]
            List of two parameters to plot.
        ax : matplotlib.axes.Axes, optional
            Axes to plot on. If None, a new figure and axes will be created.
        xrange : Tuple[float, float], optional
            Range for the x-axis. If None, the range will be determined from the data.
        yrange : Tuple[float, float], optional
            Range for the y-axis. If None, the range will be determined from the data.
        color : str, optional
            Color for the plot. Default is "k" (black).
        confidence_fraction : float, optional
            Confidence fraction for the contour level. Default is 0.9.
        gridsize : Tuple[int, int], optional
            Resolution of the grid for KDE evaluation. Default is (100, 100).
        contour_type : str, optional
            Type of contour to plot. Either 'line' for contour lines or 'filled' for filled contours.
            Default is 'line'.
        Returns
        -------
        ax : matplotlib.axes.Axes
            The axes with the plot.
        """

        if ax is None:
            f, ax = plt.subplots(1, 1, figsize=(10, 5), dpi=300)

        if transform is None:
            transform = lambda x: x

        with read_metafile(metafile) as data:

            posterior_x = np.array(data['posterior_samples'][parameters[0]]) 
            posterior_y = np.array(data['posterior_samples'][parameters[1]])

            pts = np.vstack([posterior_x, posterior_y])  
            kde = Bounded_2d_kde(transform(pts))

            if xrange is None:
                xrange = (np.min(posterior_x), np.max(posterior_x))
            if yrange is None:
                yrange = (np.min(posterior_y), np.max(posterior_y))

            x = np.linspace(xrange[0], xrange[1], gridsize[0])
            y = np.linspace(yrange[0], yrange[1], gridsize[1])
            X, Y = np.meshgrid(x, y)
            Z = kde(transform(np.vstack([X.ravel(), Y.ravel()]))).reshape(X.shape)

            # Calculate the 90% HPD level
            # Sort the flattened Z values in descending order
            sorted_z = np.sort(Z.ravel())[::-1]
            # Calculate cumulative probability
            cumsum = np.cumsum(sorted_z)
            cumsum /= cumsum[-1]  # Normalize to 1
            # Find the density level corresponding to the desired credible fraction
            level_cred = sorted_z[np.searchsorted(cumsum, confidence_fraction)]

            if contour_type == 'filled':
                try:
                    contour = ax.contourf(X, Y, Z, levels=[level_cred, Z.max()], colors=[color], alpha=0.3)
                except ValueError:
                    pass
            elif contour_type == 'line':
                ax.contour(X, Y, Z, levels=[level_cred], colors=color, labels=[label] if label else None)
            else:
                raise ValueError("contour_type must be either 'line' or 'filled'")
            ax.set_xlabel(display_parameters.get(parameters[0], {}).get('symbol', parameters[0]))
            ax.set_ylabel(display_parameters.get(parameters[1], {}).get('symbol', parameters[1]))

            # Shade out unphysical region where mass_2 > mass_1
            if (parameters[0] in ['mass_1', 'mass_1_source'] and 
                parameters[1] in ['mass_2', 'mass_2_source']):
                # Create mask where Y > X (mass_2 > mass_1)
                equal = np.linspace(xrange[0], xrange[1], 200)
                ax.fill_between(equal, equal, yrange[1], color='gray', alpha=0.3, 
                               zorder=10,)
                
            ax.set_xlim(xrange)
            ax.set_ylim(yrange)
                
                
        return ax



def comparison_plots_metafile(project, event_data, reviewed=False):
    """
    Compare results from a metafile.
    """


    properties = ['chirp_mass_source', 'mass_ratio', r'chi_eff', 'chi_p', 'luminosity_distance']
    #ias_property = {"chirp_mass_source": "mchirp", "mass_ratio": "lnq", "chi_eff": None, "chi_p": None, "luminosity_distance": "d_luminosity"}
    labels = [r'$\mathcal{M}/{\rm M}_{\odot}$', r'$q$', r'$\chi_{\rm eff}$', r'$\chi_{\rm p}$', r'$D_{\rm L}/{\rm Mpc}$']
    axis_types = ["l", "n", "n", "n", "l"]
    ranges = [[5,50], [0, 1], [-1, 1], [0, 1], [10, 20000]]


    colours = [
        "#3D49E3",
        "#74E33D",
        "#E36F3D",
        "#8E6756",
        "#88C6D4",
        "#C88EC7",
        "#48070E",
        "#FC851D",
        "#4DCA70",
        "#5BB75A",
        "#0D2421",
        "#A91031",
        "#3062D3",
        "#D7E83B",
        "#741D0A",
        "#BE8FCD"
    ]
    row = 0
    rows = []
    for event in event_data:
        for analysis in event['productions']:
            if not isinstance(analysis, dict) or not analysis:
                continue
            items = list(analysis.items())
            if not items:
                continue
            analysis_name, analysis_val = items[0]
            if "bilby" not in str(analysis_val.get('pipeline', '')).lower():
                continue
            rows.append(
                {"event": event['name'],
                 "analysis": analysis_name,
                 "metafile": get_results_file(project, event['name'], analysis_name),
                 "name": analysis_name}
            )
    
    fig = plt.figure(layout=None, figsize=(1.618*5, 0.33*len(rows)), dpi=300)
    gs = fig.add_gridspec(nrows=len(rows),
                          ncols=1+len(properties))

    events_so_far = []
    axes = {}            
    for row, event in enumerate(rows):
        events_so_far += [event['event']]
        events_so_far = list(set(events_so_far))


        analysis = event['analysis']
        ax = fig.add_subplot(gs[row, 0])
        ax.set(yticklabels=[])
        if not row-1 == len(rows):
            ax.tick_params(left=False, bottom=False, top=False, labelsize=6)
            ax.set(xticklabels=[])
        else:
            ax.spines['bottom'].set_visible(True)
            ax.tick_params(left=False, bottom=True, labelsize=5, labelrotation=90)
            ax.tick_params(which='major', width=1.00, length=5)
            ax.tick_params(which='minor', width=0.75, length=2.5)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.tick_params(left=False, bottom=False, labelsize=6)
        ax.grid(False)
        ax.set(yticklabels=[])
        ax.tick_params(left=False, bottom=False, labelsize=6)
        ax.spines[['bottom', 'left', 'right', 'top']].set_visible(False)
        ax.grid(False)

        ax.text(0,0.5, s=f"{event['event']}", fontname="serif", verticalalignment="center", size=7)
        for col, property in enumerate(properties):
            if property == "chirp_mass_source":
                axes[property] = ax = fig.add_subplot(gs[row, 1:2])
            else:
                axes[property] = ax = fig.add_subplot(gs[row, col+1])

            #ax.tick_params(top=False, right=False, left=False, bottom=False, labelsize=6)
            ax.grid(False)
            ax.spines[['bottom', 'left', 'right', 'top']].set_visible(False)

            print(row, event['event'], len(rows))

            if row+1 == len(rows):
                ax.tick_params(which='major', width=1.00, length=5, labelsize=6, right=False, left=False, labelrotation=90)
                ax.tick_params(which='minor', width=0.75, length=2.5, right=False, left=False)
            else:
                ax.tick_params(which='major', width=1.00, length=5, labelcolor="#00000000", bottom=False, left=False, right=False, top=False)
                ax.tick_params(which='minor', width=0.75, length=2.5, bottom=False, right=False, left=False)
            ax.set(yticklabels=[])


        if os.path.exists(event['metafile']):
            print(f"Reading {event['metafile']}")
            with h5.File(event['metafile']) as metafile:
                posterior_data = metafile
                analyses = list(metafile.keys())
                if "history" in analyses:
                    analyses.remove("history")
                if "version" in analyses:
                    analyses.remove("version")

                for col, property in enumerate(properties):
                    ax = axes[property]
                    if row == 0: 
                        ax.set_title(labels[col], fontsize=10)

                    if axis_types[col] == "l":
                        ax.set_xscale("log")
                        axis = np.logspace(np.log10(ranges[col][0]), np.log10(ranges[col][1]), 400)
                        ax.xaxis.set_major_locator(ticker.LogLocator(base=10, numticks=15))
                    else:
                        axis = np.linspace(ranges[col][0], ranges[col][1], 400)
                        if property in {"chi_p", "chi_eff", "mass_ratio"}:
                            ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
                            ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))

                    ax.set_xlim([axis[0],axis[-1]])
                    
                    posterior = posterior_data[event['name']]['posterior_samples']
                    color = colours[(len(events_so_far)-1)%len(colours)] # if analysis == "Beyond" else "k"
                    #hatch = analysis_hatching[analysis] if analysis != "Beyond" else None
                    above = 1 # if analysis == "Beyond" else -1
                    if property in [col for col, _ in posterior.dtype.descr]:
                        kde = GaussianKDE(posterior[property])(axis)
                        kde = kde / np.max(np.abs(kde))
                        ax.fill_between(axis, (above)*kde, color=color, alpha=0.3)#, hatch=hatch)
                    else:
                        print(f"{property} not in the samples dict for {analysis}/{event}")
                    ax.set_ylim([-1.2,1.2])

                    print("\tPlotting priors")
                    priors = posterior_data[event['name']]['priors']['samples']

                    kde = GaussianKDE(priors[property])(axis)
                    kde = kde / np.max(np.abs(kde))
                    ax.plot(axis, kde, color="black", linewidth=1, linestyle='dashed', path_effects=[pe.Stroke(linewidth=1.5, foreground='w'), pe.Normal()])
    
                ax.set_ylim([-1.2,1.2])

        
    return fig