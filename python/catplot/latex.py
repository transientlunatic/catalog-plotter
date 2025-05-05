"""
Utilities to produce latex-format representations of the catalogue.
"""

from .analysis import calculate_property

def sanitise_macro_name(name):
    """
    Convert a string to something which can be used for latex macro names.
    """
    name = name.replace('_','').replace('1', 'one').replace('2', 'two')

def event_value_pair(event, value):

    if precision.get(property) == ".1f":
        macro = "{{{0}}}{{{1:.1f}}}"
    elif precision.get(property) == ".2f":
        macro = "{{{0}}}{{{1:.2f}}}"
    elif precision.get(property) == "g":
        macro = "{{{0}}}{{{1:g}}}"
    elif precision.get(property) == "round3":
        value = int(value.round(-1))
        macro = "{{{0}}}{{{1}}}"
    else:
        macro = "{{{0}}}{{{1:g}}}"
        
    return macro.format(event, value)


def properties_latex_macros(event_data, analyses):

    properties = ['total_mass_source',
                  'chirp_mass_source',
                  'mass_1_source',
                  'mass_2_source',
                  r'chi_eff',
                  'luminosity_distance',
                  'redshift',
                  #'final_mass_source_non_evolved', 'final_spin_non_evolved'
                  ]
    
    precision = {"total_mass_source": ".1f",
                 "chirp_mass_source": ".1f",
                 "mass_1_source": ".1f",
                 "mass_2_source": ".1f",
                 "chi_eff": ".2f",
                 "luminosity_distance": "round3",
                 "redshift": ".1f",
                 "final_mass_source_non_evolved": ".1f",
                 "final_spin_non_evolved": ".2f",
                 }

    macros = ""
    table = ""

    lowers = {}
    uppers = {}
    medians = {}

    for property in properties:
        # This should probably get moved to a function specific to sanitising latex outputs
        p = sanitise_macro_name(property)
        for row, event in enumerate(event_data):
            print(event)
            print(event_data[event][analyses]["metafile"])

            analyses = event_data[event][analyses]["name"]
            datafile = event_data[event][analyses]["metafile"]
            
            with h5.File(datafile) as metafile:
                if analyses is not None and len(analyses) == 1:
                    posterior = metafile[analyses[0]]['posterior_samples']
                elif analyses is None:
                    analysis = list(metafile.keys())
                    analysis.remove("history")
                    analysis.remove("version")
                    analysis = analysis[0]
                    posterior = metafile[analysis]['posterior_samples']

                lower, median, upper = calculate_property(posterior[property])

                lowers[event] = lower
                uppers[event] = upper
                medians[event] = median

                    
        for event in event_data:
            macro = f"\\DeclareRobustCommand{{\{p}lower}}[1]{{\IfEqCase{{#1}}{{"
            for event in event_data:
                macro += event_value_pair(event, lowers[event])
            macro += "}}\n"

            macro += f"\\DeclareRobustCommand{{\{p}upper}}[1]{{\IfEqCase{{#1}}{{"
            for event in event_data:
                macro += event_value_pair(event, uppers[event])
            macro += "}}\n"

            macro += f"\\DeclareRobustCommand{{\{p}median}}[1]{{\IfEqCase{{#1}}{{"
            for event in event_data:
                macro += event_value_pair(event, medians[event])
            macro += "}}\n"

            name = f"\{p}"
            macro += f"\\DeclareRobustCommand{{{name}}}[1]{{\IfEqCase{{#1}}{{"
            for event in event_data:
                macro +=f"{{{event}}}{{${name}median{{{event}}}^{{+{name}upper{{{event}}}}}_{{{name}lower{{{event}}}}}$}}"

            macro += "}}\n"

        macros += f"{macro}\n"

    table = ""
    for event in event_data:
        row = ""
        row += f"\eventname{{{event.split('_')[0][2:]}}}{{{event.split('_')[1]}}} & "
        for property in properties:
            row += f"\{property.replace('_','').replace('1', 'one').replace('2', 'two')}{{{event}}} & "
        row = row[:-1]
        table += row
        table += "\\\\ \n"
    return macros, table
