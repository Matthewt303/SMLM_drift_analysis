"""
Created on Monday the 25th of March

@author: Matthew Tang
"""
## Import modules

import numpy as np
import matplotlib as plt
from sklearn.cluster import DBSCAN
import warnings

## Classes

class LocalisationData:

    def __init__(self, name):

        self.name = name

        self.localisation_data = []

        self.xydata = []

        self.dbscan_data = []

        self.all_beads_drift_xaxis = []

        self.all_beads_drift_yaxis = []

        self.mean_std_xaxis = float()

        self.mean_std_yaxis = float()

    def load_localisations(self, input_path):

        if isinstance(input_path, str):

            pass

        else:

            raise TypeError('Input path must be a string')

        if input_path[-4:] != '.csv':

            raise NameError('Input file must be a .csv file')

        loc_data = np.genfromtxt(input_path,
                                          dtype = float,
                                          skip_header=1,
                                          delimiter = ',')

        if loc_data.shape[1] == 9 or loc_data.shape[1] == 10:

            pass

        else:

            raise IndexError('The localisation data has too many/few columns')
            print('Please make sure you are using the localisation data from ThunderSTORM')

        self.localisation_data = loc_data

        self.xydata = loc_data[:, 2:4]

        return self.xydata

    def sort_beads(self, epsilon, cluster_size):

        if not isinstance(epsilon, (float, int)):

            raise TypeError('Epsilon must be a float or integer')

        if epsilon <= 0:

            raise ValueError('Epsilon cannot be less than or equal to zero')

        if epsilon < 15:

            warnings.warn('Epsilon is quite small. It is best to ensure it is'
                          'higher than the localisation precision')

        if not isinstance(cluster_size, int):

            raise TypeError('Cluster size must be an integer value')

        if cluster_size < 2:

            raise ValueError('Cluster size must be at least 2')

        xy_localisations = self.xydata.copy()

        if not isinstance(xy_localisations, np.ndarray):

            raise TypeError('Make sure your input data for DBSCAN'
                            'is in an array')

        if xy_localisations.shape[1] != 2:

            raise IndexError('There are too many columns in your input data.'
                             'Please make sure there are only two columns.')

        dbscan = DBSCAN(eps=epsilon, min_samples=cluster_size).fit(
            xy_localisations
        )

        labels = dbscan.labels

        # Set points labelled as noise to -1 instead of 0

        labels[(labels == 0)] = -1

        labelled_data = np.concatenate((xy_localisations, labels), axis=1)

        self.dbscan_data = labelled_data

        return self.dbscan_data

    def measure_drift_xaxis(self, sorted_bead_data):
        pass
