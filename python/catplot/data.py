import json

import h5py as h5
import scipy

from .analysis import calculate_property


class CatalogSummary:
    """
    Store summary statistics for a catalogue
    """

    def __init__(self, catdata):
        """
        """

        self.catdata = catdata
        self.data = []
        self._summary_from_samples()

        
    def _summary_from_samples(self):
        for event in self.catdata:
            data = {"name": event['name']}
            data['properties'] = {}
            with h5.File(event['metafile'], 'r') as datafile:
                post = datafile[event['analysis']]['posterior_samples']
                for prop in post.dtype.fields.keys():
                    try:
                        lower, median, upper = calculate_property(post[prop])
                        data['properties'][prop] = {"median": median,
                                                    "lower": lower,
                                                    "upper": upper}
                    except scipy.linalg.LinAlgError:
                        print(event['name'], prop)
            
            self.data.append(data)

    def save_json(self, output):
        with open(output, "w") as output_file:
            json.dump(self.data, output_file)
