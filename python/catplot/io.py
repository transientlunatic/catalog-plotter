"""
Utilities to help with reading and writing files.
"""

from contextlib import contextmanager

import yaml
import numpy as np
import h5py as h5


def read_catfile(filename):

    with open(filename) as datafile:

        data = yaml.safe_load(datafile)

    return data


def write_json(catalogue, filename):
    """
    Write a catalogue to a json representation
    """
    raise NotImplementedError

@contextmanager
def read_metafile(filename, analysis=None) -> h5.Group:
    """
    Read an HDF5-format PESummary Metafile, and return the analysis data from it.

    Parameters
    ----------
    filename : str
        Path to the HDF5 metafile.
    analysis : str, optional
        Name of the analysis to read. If None, the first analysis found will be used.

    Yields
    ------
    analysis_data : h5.Group
        The HDF5 group containing the analysis data.

    Raises
    ------
    ValueError
        If no analysis data is found in the metafile.

    Examples
    --------
    >>> with read_metafile("path/to/metafile.h5", analysis="my_analysis") as data:
    ...     print(data.keys())
    ...     print(data["chirp_mass"])
    """
    with h5.File(filename, "r") as metafile_f:
        for key in metafile_f:
            if key in ("history", "version"):
                continue
            node = metafile_f[key]
            yield node
            return
        else:
            raise ValueError("No analysis data found in the metafile.")