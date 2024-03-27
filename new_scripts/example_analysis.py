## Import

import drift_analysis_script

## Initialize and import data

print('Loading stuff...')

bead1 = drift_analysis_script.LocalisationData('run1')

data = input('Enter path to data: ')

out = input('Enter where you want data to be saved: ')

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
