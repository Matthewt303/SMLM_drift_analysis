## Import

import drift_analysis_script

## Initialize and import data

print('Loading data...')

project_name = input('Please enter project name: ')

# Initialize sm data
molecule1 = drift_analysis_script.SingleMoleculeData(project_name)

data = input('Please enter path to localisation data: ')

out = input('Please enter output path: ')

print('Loaded!')
