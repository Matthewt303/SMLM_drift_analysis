## Import

import drift_analysis_script

## Initialize and import data

print('Loading data...')

project_name = input('Please enter project name: ')

bead1 = drift_analysis_script.LocalisationData(project_name)

data = input('Please enter path to localisation data: ')

out = input('Please enter output path: ')

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

bead1.plot_bead_trajectory(out)

print('Script complete!')
