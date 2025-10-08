"""
Data types.
"""

class FloatWithUncertainty:
    """
    Represent a floating point number with independent upper and lower bounds.
    """

    def __init__(self, median, lower=None, upper=None):
        self.median = median
        self.lower = lower
        self.upper = upper
