import os
import sys
import collections
import yaml
import glob

from . import plotting
from . import binning
from . import io
from .colour import colours
from . import latex

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

    f = plotting.plane_plot(catfile)
    f.savefig("pop.pdf")

@click.option("--catfile")
@catplot.command()
def macros(catfile):
    click.echo("Creating LaTeX macros for the catalogue.")
    cat_data = io.read_catfile(catfile)
    latex.properties_latex_macros(cat_data)

    click.echo(cat_data)

@catplot.command()
def json(catfile):
    click.echo("Creating json summary for the catalogue.")
    
    
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
