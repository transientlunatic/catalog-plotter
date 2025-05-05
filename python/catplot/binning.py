"""
Code for binning and histogramming samples
"""

import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

def combine_analyses(files, parameters, bins):
    """
    this is badly named.
    """
    hist = np.zeros((len(bins[0])-1,len(bins[1])-1))
    for datafile in files:
        with h5.File(datafile) as metafile:
            analysis = list(metafile.keys())
            analysis.remove("history")
            analysis.remove("version")
            analysis = analysis[0]
            posterior = metafile[analysis]['posterior_samples']
            hist += np.histogram2d(posterior[parameters[0]], posterior[parameters[1]], bins=bins)[0]
            
    hist /= np.sum(hist)

    return hist


def contours(datafile, label, parameters, percentile, bins, colour='k'):
    """
    Return percentile contours.
    """
    with h5.File(datafile) as metafile:
        analysis = list(metafile.keys())
        analysis.remove("history")
        analysis.remove("version")
        analysis = analysis[0]
        posterior = metafile[analysis]['posterior_samples']
        hist = gaussian_kde(np.vstack([posterior[parameters[0]], posterior[parameters[1]]]))
        X, Y = np.mgrid[bins[0][0]:bins[0][-1]:len(bins[0])*1j, bins[1][0]:bins[1][-1]:len(bins[1])*1j]
        positions = np.vstack([X.ravel(), Y.ravel()])
        
        c = plt.contour(bins[0], bins[1], np.reshape(hist(positions), (len(bins[0]), len(bins[1]))).T, levels=[percentile], colors=[colour])
        ann = plt.annotate(label, xy=(np.median(posterior[parameters[0]]), np.median(posterior[parameters[1]])),
                           xytext=(np.median(posterior[parameters[0]])+np.random.randn(1), np.median(posterior[parameters[1]])-.5),
                           fontsize=6,
                           arrowprops=dict(facecolor=colour, shrink=0.05),
        )
        
        return c, ann
