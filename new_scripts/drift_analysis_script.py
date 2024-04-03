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

class BeadData:

    """
    A class to represent the localisation data arising from an SMLM image acquisition.

    Attributes
    ----------
    name: str
        a user-specified name for the dataset.
        Recommended to be the same name as TifStack
    localisation data: numpy array
        The localisation table arising from sub-pixel localisation
        of an image acquisition series.
    xydata: numpy array
        The xy localisations from the localisation table.
    dbscan_data: numpy array
        The xy localisations and the label assigned to each localisation from
        DBSCAN.
    all_beads_drift_xaxis: list
        List of arrays where each array corresponds to the x-drift of a bead
    all_beads_drift_yaxis: list
        Same as above but for y-drift
    all_beads_drift_xdrift_std: list
        List of standard deviations of the x-drift for each bead
    all_beads_drift_ydrift_std: list
        Same as above but for standard deviations of the y-drift
    mean_std_xaxis: float
        Mean of the x-drift standard deviations of all beads
    mean_std_yaxis: float
        Same as above but for y-axis

    Methods
    ------
    load_localisations(input_path):
        Loads localisation table given the input path

    dbscan_beads(epsilon, cluster_size):
        Hierarchical clustering of beads using DBSCAN, requires
        a minimum cluster radius (epsilon) and minimum cluster size (cluster_size)

    measure_drift(clustered_bead_data):
        From the clustered bead data, all localisations for a particular label
        are extracted. The coordinates of the first localisation are then subtracted
        and the resulting array is stored in a list.

        The standard deviations along x and y are also calculated from the resulting array
        and stored in a list. This process is repeated for all labels. Once it is complete,
        the mean of standard deviations along x and y are computed and stored.

    plot_bead_trajectories(outpath):
        Loops through each array (bead localisation data) stored in self.all_beads_drift_xaxis
        and generates a unique value for each localisation. Repeats for each bead. The unique
        values are converted to the experimental timecourse with units of minutes.

        A scatterplot is then drawn for all bead localisations. The color scheme is such that
        red represents the beginning of the experiment while blue represents the end.
    """

    def __init__(self, name):

        self.name = name

        self.localisation_data = []

        self.xydata = []

        self.dbscan_data = []

        self.all_beads_drift_xaxis = []

        self.all_beads_drift_yaxis = []

        self.all_beads_xdrift_std = []

        self.all_beads_ydrift_std = []

        self.mean_std_xaxis = float()

        self.mean_std_yaxis = float()

    def load_localisations(self, input_path):

        """
        This method loads localisation data given the input path.

        :param input_path: file path to localisation data
        :type input_path: string

        :return self.xydata: xy localisations of beads
        """

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

        """
        Carries out DBSCAN (doi: 10.5555/3001460.3001507) of bead localisations
        given epsilon and a cluster size.

        :param epsilon: minimum cluster radius in nanometers
        :type epsilon: float or integer
        :param cluster_size: minimum number of points that constitute a cluster
        :type cluster_size: integer

        :return self.dbscan_data: xy localisations with labels from DBSCAN
        """

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

        """
        This method measures the drift along x and y for each bead, calculates
        the standard deviation for each bead and then calculates the mean
        for all standard deviations.

        :param clustered_bead_data: localisation data with a label assigned to each localisation
        :type clustered_bead_data: numpy array

        """

        # Check input data shape
        
        if clustered_bead_data.shape[1] != 3:
            
            raise IndexError('Clustered bead data must have three columns.')
        
        # Separate beads

        bead_labels = np.unique(clustered_bead_data[:, 2])

        # Calculate standard deviation for each bead

        for label in bead_labels:

            # Ignore beads labelled as noise

            if int(label) == -1:

                pass

            else:

                bead_coords = clustered_bead_data[(clustered_bead_data[:, 2] == label)]

                # Calculate x drift

                x_drift = (bead_coords[:, 0] - bead_coords[0, 0]).reshape(
                    bead_coords[:, 0].shape[0], 1)

                x_std = np.std(x_drift[1:])

                # Calculate y drift

                y_drift = (bead_coords[:, 1] - bead_coords[0, 1]).reshape(
                    bead_coords[:, 1].shape[0], 1)

                y_std = np.std(y_drift[1:])

                # Append to attributes

                self.all_beads_drift_xaxis.append(x_drift)

                self.all_beads_drift_yaxis.append(y_drift)

                self.all_beads_xdrift_std.append(x_std)

                self.all_beads_ydrift_std.append(y_std)

        # Check number of std devs matches number of labels
        
        if len(self.all_beads_ydrift_std) != bead_labels.shape[0] - 1:
            
            raise IndexError('Number of standard deviations should match number of labels.'
                             'Check your y-values again.')
        
        if len(self.all_beads_xdrift_std) != bead_labels.shape[0] - 1:
            
            raise IndexError('Number of standard deviations should match number of labels.'
                             'Check your x-values again.')
        
        self.mean_std_xaxis = np.mean(np.vstack(self.all_beads_xdrift_std))
        
        self.mean_std_yaxis = np.mean(np.vstack(self.all_beads_ydrift_std))

    def save_data(self, outpath):

        """
        This method saves the output from DBSCAN and the mean drift along
        the x- and y-axes as a .txt file.

        :param outpath: output path where data will be saved
        :type outpath: string
        :return:
        """

        if not isinstance(outpath, str):

            raise TypeError('Output path is not a string.')

        np.savetxt(outpath + '/sorted_beads.txt', self.dbscan_data, fmt='%.5e',
                   header='Beads localisations sorted by label \n'
                          'Mean drift in x = ' + str(self.mean_std_xaxis) + 'nm \n'
                          'Mean drift in y = ' + str(self.mean_std_yaxis) + 'nm \n'
                          'x[nm] y[nm] label \n'
                   )

    def plot_bead_trajectory(self, outpath):

        """
        This method plots a scatterplot and saves it to local storage, as specified
        by the outpath as a .png and a .tif. The scatterplot consists of the drift
        values along x and y for all beads, as calculated by measure_drift. The colors
        represent the experimental timecourse, with red indicating the start and blue
        indicating the end.


        :param outpath: output path where data will be saved
        :type outpath: string
        """

        # Extract x and y coordinates

        x = self.all_beads_drift_xaxis

        y = self.all_beads_drift_yaxis

        # Frame id is a unique value assigned to each localisation for a particular bead

        frame_ids = []

        # Extract frame ids for each bead

        for bead in x:

            ids = np.arange(1, len(bead) + 1).reshape(len(bead), 1)

            frame_ids.append(ids)

        # Compile all x, y, and frame id values

        x_vals, y_vals, frames = np.vstack(x), np.vstack(y), (np.vstack(frame_ids) - 1) / 2

        # Check all arrays are the same size

        if x_vals.shape[0] == y_vals.shape[0] == frames.shape[0]:

            pass

        else:

            raise ValueError('Sizes of x-values/y-values/frames are not equal')

        # Plot scatterplot

        mpl.rcParams['font.family'] = 'sans-serif'
        mpl.rcParams['font.size'] = 13

        plt.figure(figsize=(10, 10), dpi=500)
        plt.scatter(x_vals, y_vals, c=frames, cmap=plt.cm.RdYlBu, marker='o',s=5, alpha=0.8)
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
        
        ax.set_xlim([np.min(x_vals) -5, np.max(x_vals) + 5])
        ax.set_ylim([np.min(y_vals) - 5, np.max(y_vals) + 5])
        ax.set_xlabel('x (nm)', labelpad=12, fontsize=24)
        ax.set_ylabel('y (nm)', labelpad=12, fontsize=24)

        plt.savefig(outpath + '/drift_trajectory.tif')
        plt.savefig(outpath + '/drift_trajectory.png')

class SingleMoleculeData:

    def __init__(self, name):

        self.name = name

        self.xydata = []

        self.sm_xdrift = []

        self.sm_ydrift = []

        self.xdrift_std = float()

        self.ydrift_std = float()

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

        if loc_data.shape[1] != 2:

            raise IndexError('Input file must have two columns.')

        self.xydata = loc_data

    def calculate_drift(self):

        x_coords = self.xydata[:, 0].copy()

        y_coords = self.xydata[:, 1].copy()

        x_drift = (x_coords - x_coords[0]).reshape(x_coords.shape[0], 1)

        self.sm_xdrift = x_drift

        y_drift = (y_coords - y_coords[0]).reshape(y_coords.shape[0], 1)

        self.sm_ydrift = y_drift

        self.xdrift_std = np.std(x_drift[1:])

        self.ydrift_std = np.std(y_drift[1:])

    def plot_sm_drift(self, outpath):

        x = self.sm_xdrift.copy()

        y = self.sm_ydrift.copy()

        frames = np.arange(1, x.shape[0] + 1).reshape(x.shape[0], 1)

        mpl.rcParams['font.family'] = 'sans-serif'
        mpl.rcParams['font.size'] = 9

        plt.figure(figsize=(12, 12), dpi=500)
        plt.scatter(x, y, c=frames, cmap=plt.cm.RdYlBu, marker='o', s=30, alpha=0.8)
        plt.box(True)

        # Add color bar

        colorbar = plt.colorbar()
        colorbar.set_ticks(np.arange(np.min(frames), np.max(frames), 20))
        colorbar.update_ticks()

        colorbar.ax.tick_params(width=1, length=3, labelsize=14)
        colorbar.set_label('Time (frame)', fontsize=16)

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

        ax.set_xlim([np.min(x) - 10, np.max(x) + 10])
        ax.set_ylim([np.min(y) - 10, np.max(y) + 10])
        ax.set_xlabel('x (nm)', labelpad=12, fontsize=24)
        ax.set_ylabel('y (nm)', labelpad=12, fontsize=24)

        plt.savefig(outpath + '/singlemolecule_drift_trajectory.tif')
        plt.savefig(outpath + '/singlemolecule_drift_trajectory.png')
