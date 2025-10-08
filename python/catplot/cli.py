import os

import click
import yaml
import matplotlib.pyplot as plt

from . import data, io, latex, plotting
from .plotting import Posterior2D, ms2q

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

@click.option("--parameters", "-p", default=("mass_1_source", "mass_2_source"), multiple=True, help="Parameters to plot")
@click.option("--catfile", help="Path to the catalogue file")
@catplot.command(help="Create a population plot of all of the events in the catalogue.")
def population(catfile, parameters):
    """
    Create a population plot of all of the events in the catalogue.

    Parameters
    ----------
    catfile : str
        Path to the catalogue file.
    parameters : tuple of str
        Parameters to plot.
    """
    click.echo("Plotting a population plane plot")

    catfile = catplot.io.read_catfile(catfile)

    if (parameters[0] in ['mass_1', 'mass_1_source'] and 
                parameters[1] in ['mass_2', 'mass_2_source']):
        transform = ms2q
    else:
        transform = None

    f, ax = plt.subplots(1,1, figsize=(10,5), dpi=300)
    for event in catfile['events']:
        if event.get('highlight', False):
            color = event.get("color", "red")
            linewidth=2
        else:
            color="k"    
            linewidth=1
        ax = Posterior2D.plot(
            event['metafile'], 
            parameters, 
            xrange=[1,15],
            yrange=[0.1, 5],
            ax=ax, 
            contour_type='filled',
            color=color,
            confidence_fraction=0.9,
            transform=transform,
            label=event['name'],
            )

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
