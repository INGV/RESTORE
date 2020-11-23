<h1>RESTORE: REal catalogs STOchastic REplenishment</h1>

Written by Angela Stallone with help from Giuseppe Falcone

For any comment, question or suggestion write to:
<angela.stallone@ingv.it>

This project has been founded by the Seismic Hazard Center
(Centro di Pericolosità Sismica, CPS, at the Istituto Nazionale di Geosica e Vulcanologia, INGV)

!!! Latest revision: November, 2020 !!!

To cite:
Stallone A., Falcone G. 2020. *Missing earthquake data reconstruction in the space-time-magnitude domain.*
Preprint on <https://essoar.org> (2020) DOI: 10.1002/essoar.10504916.1

Zenodo repository:
<https://doi.org/10.5281/zenodo.3952182>

<h2>About</h2>

RESTORE is a powerful tool to takle the short-term aftershock incompleteness issue (STAI).
It is based on a stochastic gap-filling procedure, which reconstructs the missing events in the space-time-magnitude domain based on empirical earthquake properties. The subsets of the catalog that are affected by the STAI issue are automatically detected.


**RESTORE files**:

*Run_RESTORE.py* loads the seismic catalog, sets the input parameters and runs *RESTORE.py*

*RESTORE.py* is the main module

*Synthetic_example_ETAS.py* runs the synthetic test (uses *ETAS_incomplete.txt* as input dataset)


**Dataset included**:

*ETAS_complete.txt* is the synthetic dataset (before STAI modeling)

*ETAS_incomplete.txt* is the synthetic dataset (after STAI modeling)

Compare the replenished catalog with *ETAS_complete.txt*, to check how the missing events are reconstructed by RESTORE

<h2>Requirements</h2>

statsmodels >= 0.12.1

matplotlib < 3.0.0



