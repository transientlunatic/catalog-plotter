from itertools import combinations

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

def calculate_divergences(datafiles, analyses=None):
    """
    Calculate the divergences between all of the analyses in a given set of datafiles and analyses.

    Note
    ----

    Adapting this from some old paper-specific code is very much a work in progress!
    """
    
    divergences = {}
    posteriors = {}

    for datafile in datafiles:
    
        with h5.File(datafile) as metafile:
            if analyses is not None and len(analyses) == 1:
                posterior = metafile[analyses[0]]['posterior_samples']
            elif analyses is None:
                analysis = list(metafile.keys())
                analysis.remove("history")
                analysis.remove("version")
                analysis = analysis[0]
                posterior = metafile[analysis]['posterior_samples']
            posteriors[analysis] = posterior

    ranges = [[0,1000], [0, 1], [-1, 1], [0, 1], [0, 25000]]

    divergences = {}
    for posterior_a, posterior_b in combinations(posteriors):
        properties = set(posterior_a.keys()) & set(posterior_b.keys())
        for property in properties:
            posterior_a_kde = GaussianKDE(posterior_a[property])(axis)
            posterior_b_kde = GaussianKDE(posterior_b[property])(axis)

            try:
                divergences[event][property] = scipy.spatial.distance.jensenshannon(posterior_a, posterior_b)
            except Exception as e:
                #print(property, e)
                divergences[event][property] = np.ma.masked        

    return divergences


def calculate_property(posterior, percentile=90):
    """
    Calculate the median and credible intervals for a given posterior.
    """
    
    interval = hpd_two_sided_credible_interval(posterior, percentile)[0]
    median = np.median(posterior)

    return interval[0]-median, median, interval[1]-median
    

\

    

