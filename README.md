**RESTORE**: **RE**al catalogs **STO**chastic **RE**plenishment

RESTORE is a powerful tool to takle the short-term aftershock incompleteness issue (STAI).
It is based on a stochastic gap-filling procedure, which reconstructs the missing events in the space-time-magnitude domain based on empirical earthquake properties.
The subsets of the catalog that are affected by the STAI issue are automatically detected.


**RESTORE files**:

*Run_RESTORE.py* loads the seismic catalog, sets the input parameters and runs *RESTORE.py*

*RESTORE.py* is the main module


**Dataset included**:

*ETAS_complete.txt* is the synthetic dataset (before STAI modeling)

*ETAS_incomplete.txt* is the synthetic dataset (after STAI modeling)

Compare the replenished catalog with *ETAS_complete.txt*, to check how the missing events are reconstructed by RESTORE


|----------------------------------------|

How to cite this program in publications:

<a href="https://zenodo.org/badge/latestdoi/281066620"><img src="https://zenodo.org/badge/281066620.svg" alt="DOI"></a>

