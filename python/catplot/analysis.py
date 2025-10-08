
import numpy as np
from pesummary.utils.credible_interval import hpd_two_sided_credible_interval


def calculate_property(posterior, percentile=90):
    """
    Calculate the median and credible intervals for a given posterior.
    """
    
    interval = hpd_two_sided_credible_interval(posterior, percentile)[0]
    median = np.median(posterior)

    return interval[0]-median, median, interval[1]-median