from asyncio import events
import os

import click
import yaml
import matplotlib.pyplot as plt

from . import data, io, latex, plotting
from .plotting import Posterior2D, SmallMultiples
from .transforms import ms2q

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
        event_data = cat_data['events']
        
    # figure = plotting.comparison_plots_metafile(project=project, event_data=event_data, catfile=cat_data)
    f = SmallMultiples.plot(event_data);
    f.savefig("family.png")

def _load_events(project, catfile, approximant, pipeline, exclude=()):
    """Events with approved analyses from an asimov ledger, or from a catfile."""
    from . import asimov as asimov_utils

    if project is not None:
        criteria = [[{"pipeline": pipeline}, {"waveform.approximant": approximant}]]
        events = asimov_utils.get_approved_event_data(
            project=project,
            ledger_path=os.path.join(project, ".asimov", "ledger.yml"),
            analyses_criteria=criteria,
        )
        click.echo(f"Found {len(events)} events with approved {approximant} analyses.")
    elif catfile is not None:
        events = io.read_catfile(catfile)["events"]
    else:
        raise click.UsageError("Provide either --project or --catfile.")
    if exclude:
        events = [event for event in events if event["name"] not in exclude]
        click.echo(f"Excluded {', '.join(exclude)}.")
    return events


@catplot.command()
@click.option("--project", default=None, help="Asimov project; uses the approved analysis for each event.")
@click.option("--catfile", default=None, help="Catalogue file (alternative to --project).")
@click.option("--approximant", default="IMRPhenomXPHM", show_default=True)
@click.option("--pipeline", default="bilby", show_default=True)
@click.option("--exclude", multiple=True, help="Event name to leave out (repeatable).")
@click.option("--output", "-o", default="smallmultiples.pdf", show_default=True,
              help="Output file. With several pages, -1, -2... is added to the name.")
@click.option("--sort-by", type=click.Choice(["snr", "name"]), default="snr", show_default=True)
@click.option("--rows-per-page", type=int, default=0, help="Split over pages of at most this many rows.")
@click.option("--jsd/--no-jsd", default=False, help="Add a Jensen-Shannon divergence bar to each cell.")
@click.option("--colour", default="#2b5c8a", show_default=True, help="Fill colour, or 'cycle' for one per row.")
@click.option("--width", type=float, default=7.0, show_default=True, help="Figure width in inches.")
def smallmultiples(project, catfile, approximant, pipeline, exclude, output, sort_by,
                   rows_per_page, jsd, colour, width):
    """
    Small multiples of the marginal posteriors, one row per event.
    """
    from . import multiples

    events = _load_events(project, catfile, approximant, pipeline, exclude)
    summaries = multiples.sort_summaries(multiples.summarise_events(events), sort_by)
    pages = multiples.split_pages(summaries, rows_per_page)
    root, ext = os.path.splitext(output)
    for number, page in enumerate(pages, start=1):
        fig = multiples.plot_small_multiples(page, jsd=jsd, colour=colour, width=width)
        path = output if len(pages) == 1 else f"{root}-{number}{ext}"
        fig.savefig(path, dpi=300)
        plt.close(fig)
        click.echo(f"Saved {len(page)} events to {path}")


@catplot.command()
@click.option("--project", default=None)
@click.option("--catfile", default=None)
@click.option("--approximant", default="IMRPhenomXPHM", show_default=True)
@click.option("--pipeline", default="bilby", show_default=True)
@click.option("--exclude", multiple=True, help="Event name to leave out (repeatable).")
@click.option("--output", "-o", default="divergence.pdf", show_default=True)
@click.option("--width", type=float, default=7.0, show_default=True)
def divergence(project, catfile, approximant, pipeline, exclude, output, width):
    """
    Jensen-Shannon divergence of each posterior from its prior, against network SNR.
    """
    from . import multiples

    events = _load_events(project, catfile, approximant, pipeline, exclude)
    summaries = multiples.summarise_events(events)
    fig = multiples.plot_divergence_snr(summaries, width=width)
    fig.savefig(output, dpi=300)
    click.echo(f"Saved to {output}")


@click.argument("blueprint")
@click.option("--project", required=True, help="Path to the asimov project directory")
@click.option("--output", "-o", default="family.png")
@catplot.command()
def projectanalysis(blueprint, project, output):
    """Generate a small multiples plot from a ProjectAnalysis blueprint.

    BLUEPRINT is the path to a YAML file with kind: ProjectAnalysis.
    """
    from . import asimov as asimov_utils

    with open(blueprint) as f:
        bp = yaml.safe_load(f)

    subjects = bp.get("subjects", None)
    analyses_criteria = bp.get("analyses", None)
    properties = bp.get(
        "properties",
        ["chirp_mass_source", "mass_ratio", "chi_eff", "chi_p", "luminosity_distance"],
    )
    priors = bp.get("priors", True)

    ledger_path = os.path.join(project, ".asimov", "ledger.yml")
    events = asimov_utils.get_approved_event_data(
        project=project,
        ledger_path=ledger_path,
        subjects=subjects,
        analyses_criteria=analyses_criteria,
    )

    click.echo(f"Found {len(events)} events with approved analyses.")
    if not events:
        click.echo("No events found matching the criteria.")
        return

    f = SmallMultiples.plot(events, properties=properties, priors=priors)
    f.savefig(output, dpi=bp.get("dpi", 300))
    click.echo(f"Saved to {output}")


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
    macros, table = latex.properties_latex_macros(cat_data['events'])

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
