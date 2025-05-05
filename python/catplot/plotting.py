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

from .asimov import get_results_file
from . import plotting
from . import binning
from . import io
from .colour import colours


plt.rcParams["pdf.use14corefonts"] = True
plt.rcParams['font.family'] = 'monospace'
plt.rcParams["font.monospace"] = ["Sauce Code Pro"]
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


def comparison_plots_metafile(project, event_data, reviewed=False):
    """
    Compare results from a metafile.
    """

    fs = 10  # fontsize

    properties = ['chirp_mass_source', 'mass_ratio', r'chi_eff', 'chi_p', 'luminosity_distance']
    #ias_property = {"chirp_mass_source": "mchirp", "mass_ratio": "lnq", "chi_eff": None, "chi_p": None, "luminosity_distance": "d_luminosity"}
    labels = [r'$\mathcal{M}/{\rm M}_{\odot}$', r'$q$', r'$\chi_{\rm eff}$', r'$\chi_{\rm p}$', r'$D_{\rm L}/{\rm Mpc}$']
    units = [r"$$", "", r"", r"", r"", r"$$"]
    axis_types = ["l", "n", "n", "n", "l"]
    ranges = [[5,50], [0, 1], [-1, 1], [0, 1], [10, 20000]]

    Nevents = len(event_data)

# #4D4E63
# #54634D
# #A382A1
# #914D2E
    colours = """#3D49E3
#74E33D
#E36F3D
#8E6756
#88C6D4
#C88EC7
#48070E
#FC851D
#4DCA70
#5BB75A
#0D2421
#A91031
#3062D3
#D7E83B
#741D0A
#BE8FCD""".split()

    row = 0
    rows = []
    for event in event_data:
        
        for analysis in event['productions']:
            analysis_name, analysis = list(analysis.items())[0]
            if not "bilby" in analysis['pipeline'].lower(): continue
            rows.append(
                {"event": event['name'],
                 "analysis": analysis_name,
                 "metafile": get_results_file(project, event['name'], analysis_name),
                 "name": analysis_name}
            )
    
    fig = plt.figure(layout=None, figsize=(1.618*5, 0.33*len(rows)), dpi=300)
    gs = fig.add_gridspec(nrows=len(rows),
                          ncols=1+len(properties),
                          # left=0.05, right=0.75,
                          hspace=0.0, wspace=0.12)

    xticks = {
        "chirp_mass_source": np.linspace(10, 50, 5)
    }

    events_so_far = []
    axes = {}            
    for row, event in enumerate(rows):
        events_so_far += [event['event']]
        events_so_far = list(set(events_so_far))


        analysis = event['analysis']
        axes['name'] = ax = fig.add_subplot(gs[row, 0])
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
        beyond_only=True
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


        skips = []
        if os.path.exists(event['metafile']):
            print(f"Reading {event['metafile']}")
            beyond_only = False
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
                    alpha = 0.7 # if analysis == "Beyond" else 0.3
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
        #ax.clabel(cs, cs.levels, fmt=lambda x: f"{event['name']}", fontsize=5)
        
    ax.set_xlabel(property_map[properties[0]])
    ax.set_ylabel(property_map[properties[1]])
