import os
from collections.abc import Sequence
from typing import Optional, Tuple

from aiohttp import ClientPayloadError

import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from matplotlib.mlab import GaussianKDE
from matplotlib.ticker import MultipleLocator, AutoMinorLocator, LogLocator


#from analysis import calculate_divergences
from pyfonts import load_google_font, set_default_font
from scipy.stats import gaussian_kde

from .asimov import get_results_file

from .io import read_metafile
from .colour import color_cycle, color_iter

plt.rc('axes', prop_cycle=color_cycle)

font = load_google_font("Source Code Pro")
set_default_font(font)

display_parameters = {
    "theta_jn": {
        "name": "theta_JN",
        "symbol": r"$\theta_{JN}$",
        "xrange": (0, np.pi),
    },
    "chirp_mass": { 
        "name": "Chirp mass",
        "symbol": r"$\mathcal{M}$",
        "xrange": (0.8, 100),
    },
    "mass_1": { 
        "name": "Mass 1",
        "symbol": r"$m_1$ / M_\odot",
        "xrange": (0.8, 100),
    },
    "mass_2": { 
        "name": "Mass 2",
        "symbol": r"$m_2$ / M_\odot",
        "xrange": (0.8, 100),
    },
    "mass_1_source": { 
        "name": "Mass 1 (source frame)",
        "symbol": r"$m_{1, src} / M_\odot$",
        "xrange": (0.8, 100),
    },
    "mass_2_source": { 
        "name": "Mass 2 (source frame)",
        "symbol": r"$m_{2, src} / M_\odot$",
        "xrange": (0.8, 100),
    },
    "total_mass_source": { 
        "name": "Total mass (source frame)",
        "symbol": r"$M_{tot, src} / M_\odot$",
        "xrange": (0.8, 200),
    },
    "mass_ratio": {
        "name": "Mass ratio",
        "symbol": r"$q$",
        "xrange": (0.1, 1),
    },
    "chirp_mass_source": { 
        "name": "Chirp mass (source frame)",
        "symbol": r"$\mathcal{M}_{src} / M_\odot$",
        "xrange": (0.8, 100),
    },
    "chi_eff": { 
        "name": "Effective spin",
        "symbol": r"$\chi_{eff}$",
        "xrange": (-1, 1),
    },
    "chi_p": { 
        "name": "Precessing spin",
        "symbol": r"$\chi_{p}$",
        "xrange": (0, 1),
    },
    "luminosity_distance": { 
        "name": "Luminosity distance",
        "symbol": r"$d_L$ / Mpc",
        "xrange": (10, 20000),
        "scale": "log",
    },
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

    @staticmethod
    def _plot_mass_component_lines(ax, parameters, xrange, yrange, m_value=3.0, color='k', linestyle='--', zorder=20):
        """
        If parameters correspond to (total_mass*, mass_ratio) (in any order),
        plot dashed lines corresponding to m1 = m_value and m2 = m_value.

        The relations used are:
          M_tot = m1 + m2
          q = m2 / m1

        For m1 = m0 -> M = m0 G* (1 + q)
        For m2 = m0 -> M = m0 * (1 + 1/q)
        """
        total_mass_keys = {'total_mass', 'total_mass_source'}
        if 'mass_ratio' not in parameters or not any(tk in parameters for tk in total_mass_keys):
            return

        q_index = parameters.index('mass_ratio')
        q_is_x = (q_index == 0)

        q_range = xrange if q_is_x else yrange
        # M_range = yrange if q_is_x else xrange  # kept for clarity if needed later

        eps = 1e-6
        q_lo = max(q_range[0], eps)
        q_hi = q_range[1]
        if q_hi <= q_lo:
            return

        q_vals = np.linspace(q_lo, q_hi, 300)
        # m1 = m_value -> M = m_value * (1 + q)
        M_m1 = m_value * (1.0 + q_vals)
        # m2 = m_value -> M = m_value * (1 + 1/q)  (only for q > 0)
        q_nonzero = q_vals[q_vals > eps]
        if q_nonzero.size == 0:
            return
        M_m2 = m_value * (1.0 + 1.0 / q_nonzero)

        if q_is_x:
            ax.plot(q_vals, M_m1, linestyle=linestyle, color=color, linewidth=1, zorder=zorder)
            ax.plot(q_nonzero, M_m2, linestyle=linestyle, color=color, linewidth=1, zorder=zorder)
        else:
            ax.plot(M_m1, q_vals, linestyle=linestyle, color=color, linewidth=1, zorder=zorder)
            ax.plot(M_m2, q_nonzero, linestyle=linestyle, color=color, linewidth=1, zorder=zorder)


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
                
        cls._plot_mass_component_lines(ax, parameters, xrange, yrange, m_value=3.0, color='k', linestyle='--')
                
        ax.set_xlim(xrange)
        ax.set_ylim(yrange)

        return ax

class SmallMultiples: 

    @staticmethod
    def clear_axis(ax): 
        """ Clear the axis of ticks and spines. """
        ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
        for spine in ax.spines.values(): 
            spine.set_visible(False)
        ax.tick_params(which='major', length=0, width=0, labelbottom=False)
        ax.tick_params(which='minor', length=0, width=0, labelsize=8)
        ax.grid(False)

            

    @classmethod
    def _add_name_to_row(cls, fig, gs, row, event):
        ax = fig.add_subplot(gs[row, 0])
        ax.text(0, 0.5, 
                s=f"{event['name']}", 
                #fontname="serif", 
                verticalalignment="center", 
                size=7)
        cls.clear_axis(ax)
        return ax

    @staticmethod
    def _plot_1d_kde_estimate(ax, axis, data, colour="k", style="filled"):
        kde = GaussianKDE(data)(axis)
        kde = kde / np.max(np.abs(kde))
        if style == "filled":
            ax.fill_between(axis, kde, color=colour, alpha=0.3)
        elif style == "line":
            ax.plot(axis, kde, color=colour, linewidth=1, linestyle='--')
        return ax

    @classmethod
    def plot(cls, 
             events,
             properties = ['chirp_mass_source', 'mass_ratio', 'chi_eff', 'chi_p', 'luminosity_distance'],
             priors = True,
             **kwargs,
             ): 
        """ 
        Plot small multiples for a list of events and properties. 

        Parameters
        ----------
        events : list
            List of event data from catfile or elsewhere.
        properties : list
            List of properties to plot for each event.  
        priors : bool
            Whether to include prior distributions in the plots.

        Returns
        -------
        fig : matplotlib.figure.Figure
            The resulting figure object.
        """
        
        rows = len(events)
        columns = len(properties) + 1

        fig = plt.figure(layout=None, figsize=(9, 0.33*rows), dpi=300)
        gs = fig.add_gridspec(nrows=rows, ncols=columns)

        for row, event in enumerate(events): 
            event_color = next(color_iter)
            # Add the event name to the first column
            name_ax = cls._add_name_to_row(fig, gs, row, event)

            metafile = event.get('metafile', None)


            try:
                with read_metafile(event['metafile'], **kwargs) as meta: 
                    for col, prop in enumerate(properties): 
                        ax = fig.add_subplot(gs[row, col+1])
                        
                        if display_parameters[prop].get('scale', 'linear') == 'log':
                            ax.set_xscale('log')
                            axis = np.logspace(np.log10(display_parameters[prop]['xrange'][0]), np.log10(display_parameters[prop]['xrange'][1]), 400)
                        else:
                            axis = np.linspace(*np.array(display_parameters[prop]['xrange']), 400)
                        try:
                            ax = cls._plot_1d_kde_estimate(ax, axis, meta['posterior_samples'][prop], colour=event.get('color', event_color))
                        except np.linalg.LinAlgError:
                            pass

                        if priors:
                            try:
                                ax = cls._plot_1d_kde_estimate(ax, axis, meta['priors']['samples'][prop], colour='k', style="line")
                            except np.linalg.LinAlgError:
                                pass
                            except ValueError:
                                pass
                            except KeyError:
                                pass    

                        cls.clear_axis(ax)
            except(FileNotFoundError, ClientPayloadError) as e:
                print(e)
                pass


        top_row = fig.axes[1:columns]
        for col, prop in enumerate(properties):
            top_row[col].set_title(f"{display_parameters[prop]['symbol']}", fontsize=10)

        bottom_row = fig.axes[(-columns+1):]
        for col, prop in enumerate(properties):
            bottom_row[col].tick_params(bottom=True, labelbottom=True)
            if display_parameters[prop].get('scale', 'linear') == 'log':
                bottom_row[col].xaxis.set_major_locator(LogLocator(base=10.0, numticks=10))
                bottom_row[col].xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10)*0.1, numticks=10))
                bottom_row[col].tick_params(which='major', length=6, width=1, labelbottom=True)
                bottom_row[col].tick_params(which='minor', length=3, width=1, labelsize=8)
            else:
                bottom_row[col].xaxis.set_minor_locator(AutoMinorLocator())
                bottom_row[col].tick_params(which='major', length=6, width=1, labelbottom=True)
                bottom_row[col].tick_params(which='minor', length=3, width=1, labelsize=8)
                bottom_row[col].tick_params(which='minor', labelbottom=False)
        return fig
    


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