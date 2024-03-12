
![alt text](https://github.com/angystallone/Seismology_Stuff/blob/main/figures/RESTORE_logo.png?raw=true)


Written by Angela Stallone with help from Giuseppe Falcone

For any comment, question or suggestion write to:
<angela.stallone@ingv.it>

This project has been founded by the Seismic Hazard Center
(Centro di Pericolosità Sismica, CPS, at the Istituto Nazionale di Geofisica e Vulcanologia, INGV)

To cite: Stallone, A., & Falcone, G. (2021). Missing Earthquake Data Reconstruction in the Space‐Time‐Magnitude Domain. Earth and Space Science, 8(8), e2020EA001481. [https://doi.org/10.1029/2020EA001481]

**Latest revision: March, 2021** 

|---------------------------------------------------------------------------------------------------------|

Zenodo repository:

<https://doi.org/10.5281/zenodo.3952182>

<h2>About</h2>

RESTORE is a Python tool tackling the short-term aftershock incompleteness issue (STAI).
It is based on a stochastic gap-filling procedure which reconstructs the missing events in the space-time-magnitude domain based on empirical earthquake properties. 
The subsets of the catalog affected by the STAI issue are automatically detected.


**RESTORE files**:

*Run_RESTORE.py* loads the seismic catalog and the input parameters provided by the user, then runs the script *RESTORE.py*

*RESTORE.py* is the main module

*input_file.txt* file containing the input parameters


**Synthetic_Test**: [v 2.0.0]

*Run_Synthetic_Test.py* runs the synthetic test (uses *ETAS_incomplete.txt* as input dataset)

*ETAS_complete.txt* is the synthetic dataset (before STAI modeling)

*ETAS_incomplete.txt* is the synthetic dataset (after STAI modeling)

Compare the replenished catalog with *ETAS_complete.txt*, to check how the missing events are reconstructed by RESTORE

<h2>Requirements</h2>

mc_lilliefors (download it <a href="https://gitlab.com/marcus.herrmann/mc-lilliefors">here</a>)


