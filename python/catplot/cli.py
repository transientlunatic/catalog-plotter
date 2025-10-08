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
from . import data

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

@click.option("--project", default=None)
@click.option("--catfile", default=None)
@catplot.command()
def family(project, catfile):
    """
    Create a single-page plot of all of the events in the catalogue
    with each analysis.
    """
    click.echo("Creating a family photo plot for the catalogue.")

    if project is not None:
        with open(os.path.join(project, ".asimov", "ledger.yml")) as datafile:
            event_data = yaml.safe_load(datafile)
            event_data = event_data['events']
        click.echo(f"Found {len(event_data)} events in the asimov project.")
        cat_data = None

    elif catfile is not None:
        cat_data = io.read_catfile(catfile)
        project = None
        event_data = None
        
    figure = plotting.comparison_plots_metafile(project=project, event_data=event_data, catfile=cat_data)
    figure.savefig("family.png")

@click.option("--catfile")
@catplot.command()
def population(catfile):
    click.echo("Plotting a population plane plot")

    f = plotting.plane_plot(catfile)
    f.savefig("pop.pdf")

@click.option("--catfile")
@click.option("--output", "-o", default=None)
@catplot.command()
def macros(catfile, output):
    click.echo("Creating LaTeX macros for the catalogue.")
    cat_data = io.read_catfile(catfile)
    macros, table = latex.properties_latex_macros(cat_data)

    if output is None:
        output = os.getcwd()
    
    with open(os.path.join(output, "table.tex"), "w") as table_file:
        table_file.write(table)

    with open(os.path.join(output, "macros.tex"), "w") as macro_file:
        macro_file.write(macros)

@click.option("--catfile")
@click.option("--output", "-o", default=None)        
@catplot.command()
def json(catfile, output):
    click.echo("Creating json summary for the catalogue.")
    if catfile:
        cat_data = io.read_catfile(catfile)

    summary = data.CatalogSummary(cat_data)
    summary.save_json(output)
    
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
