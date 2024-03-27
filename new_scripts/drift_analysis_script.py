"""
Created on Monday the 25th of March

@author: Matthew Tang
"""
## Import modules

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
import matplotlib as mpl
from sklearn.cluster import DBSCAN
import warnings

## Classes

class LocalisationData:

    def __init__(self, name):

        self.name = name

        self.localisation_data = []

        self.xydata = []

        self.dbscan_data = []

        self.labels = []

        self.all_beads_drift_xaxis = []

        self.all_beads_drift_yaxis = []

        self.all_beads_xdrift_std = []

        self.all_beads_ydrift_std = []

        self.mean_std_xaxis = float()

        self.mean_std_yaxis = float()

    def load_localisations(self, input_path):

        # Check if input is string and if file is csv

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

        # Check shape of localisation table, should have 9 or 10 columns if ThunderSTORM was used

        if loc_data.shape[1] == 9 or loc_data.shape[1] == 10:

            pass

        else:

            raise IndexError('The localisation data has too many/few columns')
            print('Please make sure you are using the localisation data from ThunderSTORM')

        self.localisation_data = loc_data

        self.xydata = loc_data[:, 2:4]

        return self.xydata

    def dbscan_beads(self, epsilon, cluster_size):

        # Check if inputs are valid. Epsilon must be float > 0. Cluster size must be int > 0

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

        # Check if localisations from previous step are the correct type.

        xy_localisations = self.xydata.copy()

        if not isinstance(xy_localisations, np.ndarray):

            raise TypeError('Make sure your input data for DBSCAN'
                            'is in an array')

        if xy_localisations.shape[1] != 2:

            raise IndexError('There are too many columns in your input data.'
                             'Please make sure there are only two columns.')

        # Carry out DBSCAN

        dbscan = DBSCAN(eps=epsilon, min_samples=cluster_size).fit(
            xy_localisations
        )

        labels = dbscan.labels_

        self.labels = labels

        # Set points labelled as noise to -1 instead of 0

        labels[(labels == 0)] = -1

        labelled_data = np.concatenate((xy_localisations, labels[:, np.newaxis]), axis=1)

        self.dbscan_data = labelled_data

        return self.dbscan_data

    def measure_drift(self, clustered_bead_data):

        # Check input data shape
        
        if clustered_bead_data.shape[1] != 3:
            
            raise IndexError('Clustered bead data must have three columns.')
        
        # Separate beads

        bead_labels = np.unique(clustered_bead_data[:, 2])

        # Calculate standard deviation for each bead

        for label in bead_labels:

            bead_coords = clustered_bead_data[(clustered_bead_data[:, 2] == label)]

            # Calculate x drift

            x_drift = bead_coords[:, 0] - bead_coords[0, 0]

            x_std = np.std(x_drift[1:])

            # Calculate y drift

            y_drift = bead_coords[:, 1] - bead_coords[0, 1]

            y_std = np.std(y_drift[1:])

            # Append to attributes

            self.all_beads_drift_xaxis.append(x_drift)

            self.all_beads_drift_yaxis.append(y_drift)

            self.all_beads_xdrift_std.append(x_std)

            self.all_beads_ydrift_std.append(y_std)

        # Check number of std devs matches number of labels FIX ERROR HERE
        
        if self.all_beads_ydrift_std.shape[0] != bead_labels.shape[0]:
            
            raise IndexError('Number of standard deviations should match number of labels.'
                             'Check your y-values again.')
        
        if self.all_beads_xdrift_std.shape[0] != bead_labels.shape[0]:
            
            raise IndexError('Number of standard deviations should match number of labels.'
                             'Check your x-values again.')
        
        self.mean_std_xaxis = np.mean(np.vstack(self.all_beads_xdrift_std))
        
        self.mean_std_yaxis = np.mean(np.vstack(self.all_beads_ydrift_std))

    def save_data(self, outpath):

        if not isinstance(outpath, str):

            raise TypeError('Output path is not a string.')

        np.savetxt(outpath + '/sorted_beads.txt', self.dbscan_data, fmt='%.5e',
                   header='Beads localisations sorted by label \n'
                          'x[nm] y[nm] label \n'
                          'Mean drift in x = ' + str(self.mean_std_xaxis) + 'nm \n'
                          'Mean drift in y = ' + str(self.mean_std_yaxis) + 'nm \n')

    def plot_bead_trajectory(self, outpath):

        # Extract x and y coordinates

        x = self.all_beads_drift_xaxis

        y = self.all_beads_drift_yaxis

        # Frame id is a unique value assigned to each localisation for a particular bead

        frame_ids = []

        # Extract frame ids for each bead

        for bead in x:

            ids = np.arange(1, len(bead) + 1)

            frame_ids.append(ids)

        # Compile all x, y, and frame id values

        x_vals, y_vals, frames = np.vstack(x), np.vstack(y), np.vstack(frame_ids) / 2

        # Check all arrays are the same size

        if x_vals.shape[0] == y_vals.shape[0] == frames.shape[0]:

            pass

        else:

            raise ValueError('Sizes of x-values/y-values/frames are not equal')

        # Plot scatterplot

        mpl.rcParams['font.family'] = 'sans-serif'
        mpl.rcParams['font.size'] = 13

        plt.figure(figsize=(10, 10), dpi=500)
        plt.scatter(x_vals, y_vals, c=frames, cmap=plt.cm.RdYlBu, marker='o',s=30, alpha=0.8)
        plt.box(True)
        
        # Add color bar
        
        colorbar = plt.colorbar()
        colorbar.set_ticks(np.arange(np.min(frames), np.max(frames), 20))
        colorbar.update_ticks()

        colorbar.ax.tick_params(width=1, length=3, labelsize=14)
        colorbar.set_label('Time (mins)', fontsize = 16)
        
        # Set plot parameters
        
        ax = plt.gca()
        ax.set_aspect('equal', adjustable='box')
        ax.tick_params(axis='y', which='major', length=6)
        ax.tick_params(axis='x', which='major', length=6)
        ax.xaxis.set_major_locator(MultipleLocator(10))
        ax.yaxis.set_major_locator(MultipleLocator(10))
        ax.xaxis.label.set_color('black')
        ax.yaxis.label.set_color('black')
        
        ax.spines['bottom'].set_color('black')
        ax.spines['top'].set_color('black')
        ax.spines['right'].set_color('black')
        ax.spines['left'].set_color('black')
        ax.spines['bottom'].set_linewidth(2.0)
        ax.spines['top'].set_linewidth(2.0)
        ax.spines['right'].set_linewidth(2.0)
        ax.spines['left'].set_linewidth(2.0)
        
        ax.set_xlim([-100, 100])
        ax.set_ylim([-100, 100])
        ax.set_xlabel('x (nm)', labelpad=12, fontsize=24)
        ax.set_ylabel('y (nm)', labelpad=12, fontsize=24)

        plt.savefig(outpath + '/drift_trajectory.tif')
        plt.savefig(outpath + '/drift_trajectory.png')


## Initialize and import data

print('Loading stuff...')

bead1 = LocalisationData('run1')

data = 'C:/Users/mxq76232/PycharmProjects/SMLM_drift_analysis/.venv/SMLM_drift_analysis/example_data/bead_locs_3.csv'

out = 'C:/Users/mxq76232/Desktop'

eps = 100

clusters = 200

print('Loaded!')

## Load localisations

print('Loading data, please wait')

bead1.load_localisations(data)

print('Data loaded!')

## DBSCAN

print('Carrying out DBSCAN, please wait.')

bead1.dbscan_beads(eps, clusters)

print('Done')

## Measure drift and save data

print('Performing drift analysis')

bead1.measure_drift(bead1.dbscan_data)

bead1.save_data(out)

print('Data saved!')

## Plot

bead1.plot_bead_trajectory
