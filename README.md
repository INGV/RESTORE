
![alt text](https://github.com/angystallone/Seismology_Stuff/blob/main/figures/RESTORE_logo.png?raw=true)


<h2>About</h2>

`RESTORE` is a Python tool tackling the short-term aftershock incompleteness issue (STAI).
It is based on a stochastic gap-filling procedure which reconstructs the missing events in the space-time-magnitude domain based on empirical earthquake properties. 
The subsets of the catalog affected by the STAI issue are automatically detected.

<h2>To run</h2>

Set the input parameters in the `input_file.txt` file:

- size: moving-window size (in number of events per window) --> 1000 by default (Mignan and Woessner, 2012).
- step: moving-window step (in number of events per step) --> 250 by default (Mignan and Woessner, 2012).
- st_dev_multiplier: multiplies the Mc standard deviation (Mc Sigma), controls the confidence level for STAI gaps identification; STAI gaps are windows where mc >= mc_ok + n*sigma, where n = st_dev_multiplier; increasing the value of st_dev_multiplier results in a more conservative approach in detecting temporary deviations of Mc.
- Sigma: smoothing distance of the Gaussian kernel, controls the spread of the smoothing; smaller values will result in sharper and more localized smoothing.
- sbin: bin in the latitude and longitude direction (degrees), controls the grid resolution
- fault_length_multiplier: multiplies the fault length rupture, controls the areal extent of the subcatalog where inferences about Mc trend with time are made; smaller values will result in higher resolution of STAI gaps detection, as local seismicity is less diluted.
- t_end_quiet: ending time of the seismically quiescent period.
- b: b-value of the Gutenberg-Richter law (alternatively, it can be estimated with the function provided in RESTORE).
- alpha: significance level for the Lilliefors test.
- mc: reference value for the magnitude of completeness, if set by the user (by default, it is estimated with the function provided in RESTORE).
- depth_distribution: [`scipy.optimize.curve_fit`] distribution to fit to hypocenter depths, available options are: normal, poisson, lognormal, beta, bimodal.
- p0: [`scipy.optimize.curve_fit`] initial guess for the parameters of the hypocenter depth distribution: 'mu', 'sigma' (normal), 'mu' (poisson), 'mu', 'sigma' (lognormal), 'a', 'b' (beta), 'mu1', 'sigma1', 'A1', 'mu2', 'sigma2', 'A2' (bimodal).


Run `RESTORE` using the following command:

```bash
python Run_RESTORE.py
```

<h2>Synthetic_Test [v 2.0.0 only]</h2>

`Run_Synthetic_Test.py` --> runs the synthetic test (uses `ETAS_incomplete.txt` as input dataset)

`ETAS_complete.txt` --> synthetic dataset (before STAI modeling)

`ETAS_incomplete.txt` --> the synthetic dataset (after STAI modeling)

Compare the replenished catalog with `ETAS_complete.txt`, to check how the missing events are reconstructed by RESTORE

<h2>Required external modules</h2>

`mc_lilliefors` (download it <a href="https://gitlab.com/marcus.herrmann/mc-lilliefors">here</a>)

<h2>How to cite</h2>

If you use `RESTORE` in your research, please cite using the following citation:

```bash
@software{Stallone_RESTORE,
author = {Stallone, Angela and Falcone, Giuseppe},
title = {{RESTORE}},
url = {https://github.com/INGV/RESTORE}
}
```

For any comment, question or suggestion write to:
<angela.stallone@ingv.it>


<h2>Acknowledgements</h2>

This project has been founded by the Seismic Hazard Center
(Centro di Pericolosità Sismica, CPS, at the Istituto Nazionale di Geofisica e Vulcanologia, INGV)




