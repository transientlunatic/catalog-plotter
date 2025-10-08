"""
Utilities to help with reading and writing files.
"""

import yaml

def read_catfile(filename):

    with open(filename, "r") as datafile:

        data = yaml.safe_load(datafile)

    return data


def write_json(catalogue, filename):
    """
    Write a catalogue to a json representation
    """
    raise NotImplementedError
