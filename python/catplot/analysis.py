
import numpy as np
from pesummary.utils.credible_interval import hpd_two_sided_credible_interval

from ligo.skymap.bayestar import rasterize
from ligo.skymap.io import fits as skymap_fits

import h5py

#from ligo.skymap.postprocess import sky_area
#from ligo.skymap.distance import parameters_to_samples

def calculate_property(posterior, percentile=90):
    """
    Calculate the median and credible intervals for a given posterior.
    """
    
    try:
        interval = hpd_two_sided_credible_interval(posterior, percentile)[0]
    except ValueError:
        interval = [np.nan, np.nan]
    median = np.median(posterior)

    return interval[0]-median, median, interval[1]-median

def skymap_statistics(
        posterior_samples: h5py.Group):
    """
    Return the posterior distributions over sky area and the volume.

    Parameters
    ----------
    posterior_samples : h5py.Group
        The HDF5 group containing posterior samples with 'ra', 'dec', and 'luminosity_distance' datasets.
    
    Returns
    -------
    sky_prob : np.ndarray
        The sky localization probability map.
    vol_samples : np.ndarray
        Samples of the volume distribution.
    """
    
    ra = np.array(posterior_samples["ra"])
    dec = np.array(posterior_samples["dec"])
    dist = np.array(posterior_samples["luminosity_distance"])

    sky_prob, distmu, distsigma, distnorm = rasterize(
    ra, dec, dist,
    nside=512,  # adjust resolution
    prior=None  # or pass per-sample weights if you have them
    )

    # area50 = sky_area(sky_prob, 0.5, nest=True)
    # area90 = sky_area(sky_prob, 0.9, nest=True)

    # levels = np.linspace(0.1, 0.99, 30)
    # areas = [sky_area(sky_prob, c, nest=True) for c in levels]

    ra_s, dec_s, dist_s = parameters_to_samples(
    sky_prob, distmu, distsigma, distnorm,
    n_samples=2_000,
    nest=True
    )

    vol_samples = (4/3) * np.pi * (dist_s**3)  # Mpc^3

    return sky_prob, vol_samples
