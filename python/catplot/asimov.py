"""
Code to interact with asimov projects.
"""

import os

def get_results_file(project, event, analysis):
    """
    Find the PESummary metafile for a given event and analysis.
    """

    return os.path.join(project, "results", event, analysis, "posterior_samples.h5")
    
