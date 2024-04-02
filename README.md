# SMLM_drift_analysis
A series of python scripts used to analyse the sample drift of fluorescent beads and single molecules. These scripts were used in a project aiming to quantify and characterise sample drift from a novel reinforced optical cage system microscope. 

# About 
Drift can be a major hinderance in image acquisition during single-molecule localization molecule experiments. To this end, this project aimed to develop a novel optical setup to reduce sample drift, thereby improving the fidelity and reliability of super-resolution images from single-molecule localisation microscopy. In order to quantify drift, we collected widefield image time series' of fluorescent beads, determined the beads' sub-pixel localisations using ThunderSTORM, before applying these scripts to the localisation tables obtained from ThunderSTORM. We also quantified drift from stochastic optical reconstruction microscopy (STORM) by using the localisations of ten separate single-molecules.

While these scripts were used specifically for quantifying stability in the reinforced optical cage system, we anticipate that they can be broadly used to quantify lateral drift in any microscope. However, this project does not contain any methods to analyse axial drift but may be updated in the future to include such features.

# Drift quantification
We quantified drift by first implementing DBSCAN which assigns a label to each bead in a time series. DBSCAN was implemented to account for the varying number of bead localisations that arose from less-than-optimal signal-to-noise ratio. More specifically, each bead ought to have 240 localisations; however, variations in SNR resulted in beads having less than that, e.g. 238 localisations or 239 localisations. Therefore, DBSCAN was used to unambiguously assign each bead across the time series.

After the beads had been identified, the x- and y-localisations were extracted separately, and the initial x- and y-locations (location from frame 1) were subtracted from all x- and y-locations across time, thus obtaining the displacements of each bead along the x- and y-axis. For each bead's displacements, the standard deviations along x and y were calculated. Finally, the standard deviations across all beads was averaged to obtain the drift along x and y. For visualisation, the x and y displacements of all beads over time were plotted as a scattergraph. This procedure was repeated for 10 separate datasets to obtain the overall drift quantification of the microscope. 

# Repository structure
This repository contains several folders:
- 'example_data' contains two example datasets. These are localisation tables of bead image time series' exported from ThunderSTORM. Any input data ought to have the some format and structure as these two files.
- 'original_scripts' contains three python scripts. These were used during the project to quickly assess the results from experiments. Their use is currently deprecated and is not recommended to use them.
- 'new_scripts' contains updated python scripts that carry out the same functions as those in original_scripts. However, the updated scipts are more organised, and so it is strongly recommended to use them.
  - 'drift_analysis_script.py' contains the classes and methods for analysis, with one class for bead localisation data and another for single-molecule localisation data.
  - 'example_analysis.py' contains an example of how to use drift_analysis_script.py for drift quantification.

# Installation
Note: work in progress, please check the development branch for an updated version of the code.
