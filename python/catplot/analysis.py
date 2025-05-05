import numpy as np
import matplotlib
from matplotlib.mlab import GaussianKDE
import matplotlib.patheffects as pe
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import scipy.spatial.distance
import h5py as h5
from copy import copy
from pesummary.utils.credible_interval import hpd_two_sided_credible_interval
import os

def calculate_divergences(event_data, analyses):
    divergences = {}


    properties = ['chirp_mass_source', 'mass_ratio', r'chi_eff', 'chi_p', 'luminosity_distance']
    # ias_property = {"chirp_mass_source": "mchirp", "mass_ratio": "lnq", "chi_eff": None, "chi_p": None, "luminosity_distance": "d_luminosity"}

    ranges = [[0,1000], [0, 1], [-1, 1], [0, 1], [0, 25000]]

    Nevents = len(event_data)
    row = 0

    for row, event in enumerate(event_data):
        divergences[event] = {}
        print(event)
        for col, property in enumerate(properties):
            axis = np.linspace(ranges[col][0], ranges[col][1], 400)
            posteriors = []
            for analysis in analyses:
                datafile = event[analyses]["metafile"]#f"datafiles/{event}_{analysis}.h5"

                if os.path.exists(datafile):
                    try:
                        with h5.File(datafile) as metafile:
                            posterior = metafile['posterior_samples']
                            if (analysis == "4-OGC") and (property == "mass_ratio"):
                                posterior = GaussianKDE(1./posterior[property])(axis)
                            else:
                                posterior = GaussianKDE(posterior[property])(axis)

                    except Exception as e:
                        try:
                            with h5.File(datafile) as metafile:
                                dataset = event_data[event][analysis]['name']
                                posterior = metafile[dataset]['posterior_samples']
                                posterior = GaussianKDE(posterior[property])(axis)
                        except Exception as e:
                            print(e)
                    posteriors += [posterior]
                else:
                    pass
                    #print(datafile)

            try:
                divergences[event][property] = scipy.spatial.distance.jensenshannon(*posteriors)
            except Exception as e:
                #print(property, e)
                divergences[event][property] = np.ma.masked        

    return divergences


def calculate_property(posterior):

    interval = hpd_two_sided_credible_interval(posterior, 90)[0]
    median = np.median(posterior)

    
    return interval[0]-median, median, interval[1]-median
    



def calculate_properties(event_data, analyses):
    divergences = {}


    properties = ['total_mass_source', 'chirp_mass_source', 'mass_1_source', 'mass_2_source',r'chi_eff', 'luminosity_distance',
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
        p = property.replace('_','').replace('1', 'one').replace('2', 'two')
        for row, event in enumerate(event_data):
            print(event)
            print(event_data[event][analyses]["metafile"])
            try:
                with h5.File(event_data[event][analyses]["metafile"]) as metafile:

                    posterior = metafile['posterior_samples']

                    lower, median, upper = calculate_property(posterior[property])

                    lowers[event] = lower
                    uppers[event] = upper
                    medians[event] = median
            except Exception as e:
                with h5.File(event_data[event][analyses]["metafile"]) as metafile:

                    posterior = metafile[event_data[event]['Beyond']['name']]['posterior_samples']

                    lower, median, upper = calculate_property(posterior[property])

                    lowers[event] = lower
                    uppers[event] = upper
                    medians[event] = median

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

    

