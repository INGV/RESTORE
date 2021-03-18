
![alt text](https://github.com/angystallone/Seismology_Stuff/blob/main/figures/RESTORE_logo.png?raw=true)


Written by Angela Stallone with help from Giuseppe Falcone

For any comment, question or suggestion write to:
<angela.stallone@ingv.it>

This project has been founded by the Seismic Hazard Center
(Centro di Pericolosità Sismica, CPS, at the Istituto Nazionale di Geofisica e Vulcanologia, INGV)

**Latest revision: March, 2021** 

|---------------------------------------------------------------------------------------------------------|

Zenodo repository:

<https://doi.org/10.5281/zenodo.3952182>

<h2>About</h2>

RESTORE is a powerful tool to takle the short-term aftershock incompleteness issue (STAI).
It is based on a stochastic gap-filling procedure, which reconstructs the missing events in the space-time-magnitude domain based on empirical earthquake properties. The subsets of the catalog that are affected by the STAI issue are automatically detected.


**RESTORE files**:

*Run_RESTORE.py* loads the seismic catalog and the input parameters, then runs the script *RESTORE.py*

*RESTORE.py* is the main module


**Synthetic_Test**:

*Run_Synthetic_Test.py* runs the synthetic test (uses *ETAS_incomplete.txt* as input dataset)

*ETAS_complete.txt* is the synthetic dataset (before STAI modeling)

*ETAS_incomplete.txt* is the synthetic dataset (after STAI modeling)

Compare the replenished catalog with *ETAS_complete.txt*, to check how the missing events are reconstructed by RESTORE

<h2>Requirements</h2>

Python 3

statsmodels >= 0.12.1

Basemap

mc_lilliefors (download it here: <http://doi.org/10.5281/zenodo.4162496>)






