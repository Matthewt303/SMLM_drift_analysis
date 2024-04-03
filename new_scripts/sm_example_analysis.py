## Import

""" Module providing drift analysis classes
and methods
"""

import drift_analysis_script

## Initialize and import data

print('Loading data...')

project_name = input('Please enter project name: ')

# Initialize sm data
molecule1 = drift_analysis_script.SingleMoleculeData(project_name)

data = input('Please enter path to localisation data: ')

out = input('Please enter output path: ')

print('Loaded!')

## Load localisations

print('Loading data')

molecule1.load_localisations(data)

print('Loaded!')

## Measure drift

print('Measuring drift')

molecule1.calculate_drift()

print('Done')

## Plot

print('Plotting and saving data')

molecule1.plot_sm_drift(out)

print('Done')
