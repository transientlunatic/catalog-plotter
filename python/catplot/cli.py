import os
import sys
import collections
import yaml
import glob

from . import plotting
from . import binning
from . import io
from .colour import colours

import pickle
from pesummary.io import read

import click

# @click.version_option(asimov.__version__)
@click.group()
@click.pass_context
def catplot(ctx):
    """
    Plot gravitational wave catalogues.
    """

    pass


@click.argument("project", required=True)
@catplot.command()
def family(project):
    """
    Create a single-page plot of all of the events in the catalogue
    with each analysis.
    """
    click.echo("Creating a family photo plot for the catalogue.")
    with open(os.path.join(project, ".asimov", "ledger.yml")) as datafile:
        event_data = yaml.safe_load(datafile)
        event_data = event_data['events']
    click.echo(f"Found {len(event_data)} events in the asimov project.")
    figure = plotting.comparison_plots_metafile(project, event_data)
    figure.savefig("family.png")

@click.option("--catfile")
@catplot.command()
def population(catfile):
    click.echo("Plotting a population plane plot")

    properties = ("mass_1_source", "mass_2_source")
    
    events = io.read_catfile(catfile)
    
    import numpy as np
    import matplotlib.pyplot as plt


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
    f.savefig("pop.pdf")
    
if __name__ == "__main__":
    catplot()
    
    # print("Loading events")
    # with open("../data/events.yaml", "r") as datafile:
    #     event_data = yaml.safe_load(datafile)
        
    # events = {}
    # for event, data in event_data.items():
    #     events[event] = data
            
    # fig = plotting.comparison_plots_metafile(events)
    # fig.savefig("/data/www.astro/daniel/beneath/comparsion-plots.pdf",  bbox_inches="tight")
