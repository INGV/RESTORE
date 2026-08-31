# RESTORE

`RESTORE` is a Python tool tackling the short-term aftershock incompleteness issue (STAI). It is based on a stochastic gap-filling procedure which reconstructs the missing events in the space-time-magnitude domain based on empirical earthquake properties. The subsets of the catalog affected by the STAI issue are automatically detected.

![alt text](https://github.com/angystallone/Seismology_Stuff/blob/main/figures/RESTORE_logo.png?raw=true)

## Installation

Clone the repository:

```bash
git clone https://github.com/INGV/RESTORE.git
cd RESTORE
```

Create and activate a virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate    # Windows: .venv\Scripts\activate
```

Install RESTORE:

```bash
pip install .
```

Run RESTORE example:

```bash
restore --config example/amatrice/input_file.txt
```

Print the installed RESTORE version
```bash
restore --version
```

To generate multiple RESTORE realizations, specify the desired values of `sizes`, `steps`, and `seeds` in the input file. 
RESTORE will run all combinations of these values.

## Input Catalog Format [ZMAP] and magnitude binning

Input catalogs must follow the ZMAP format, with columns following this order:

1. longitude
2. latitude
3. year
4. month
5. day
6. magnitude
7. hypocentral depth in km
8. hour
9. minute
10. second

Magnitudes must already be binned at `ΔM = 0.1`.

## Outputs

Different seeds return different replenished catalogs. 
For a given window size, window step and seed combination, the output files are written in `outputs/<example>/realizations/size_<size>_step_<step>/seed_<seed>/`.

`Replenished_catalog.txt` contains the ten input columns plus:

- `flag = 0`: observed event
- `flag = 1`: reconstructed event

The runner writes:

- `realization_metrics.csv`: contains one row for each `(size, step, seed)` realization and summarizes the representative metrics across all realizations
- `representative_realization.txt`: the realization that is closest to the central behavior of the tested ensemble (see the note below)
- for each realization: `settings.txt` reports the input parameters for that run
- for each realization: `metadata.json` with code/environment/file provenance

**NOTE**: Each successful `(size, step, seed)` combination represents one RESTORE realization. 
RESTORE defines the central ensemble behavior as the median value of each representative metric across all successful realizations. 
It then calculates the normalized distance of each realization from the center and selects the realization with the smallest distance.
Distances are normalized by the interquartile range of each representative metric across the ensemble so that metrics with different numerical scales can be compared. 
The final distance is a weighted mean: number of reconstructed events, magnitude distribution, spatial distribution, depth distribution, and STAI characteristics are given equal total weight.

The representative metrics are:

- number of reconstructed events
- median reconstructed magnitude
- median reconstructed latitude and longitude
- median reconstructed hypocentral depth
- total STAI duration
- number of STAI gaps

## How RESTORE works

For each realization, RESTORE performs the following steps:

1. **Read and validate the input catalog.**
   Note: the catalog must be in ZMAP format and magnitudes must already be discretized at `ΔM = 0.1`.

2. **Identify the largest earthquake.**
   The event with the largest magnitude is identified as responsible for the STAI.

3. **Define the STAI region.**
   RESTORE estimates the rupture length of the largest event and constructs a rectangular region around it. The size of this region is controlled by `fault_length_multiplier` (provided by the user) and represents the area affected by the STAI.

4. **Create the spatially filtered catalog.**
   Only earthquakes located inside the STAI region are used to estimate the temporal variation of magnitude completeness, detect STAI gaps, and reconstruct missing events. The file `STAI-Filtered_Catalog.txt` contains this subset. 
   Note: this implies that the provided input catalog can cover an arbitrary area. RESTORE automatically extracts the meaningful subset based on `fault_length_multiplier`.

5. **Estimate the reference magnitude of completeness.**
   RESTORE uses earthquakes occurring before `t_end_quiet` (i.e., those occurring during the 'seismic quiescent period'), inside the STAI region, to estimate the reference completeness magnitude `Mc*`, unless the user supplies their own value.

6. **Estimate `Mc` through time.**
   The filtered catalog is analyzed with a moving window approach. RESTORE estimates `Mc` for each window using the Lilliefors test.
   Note: starting from v.3.0.0, RESTORE implements the Lilliefors test with truncated-exponential dithering, as proposed by [Stallone and Spassiani (2026)](https://doi.org/10.26443/seismica.v5i1.2257).

7. **Detect STAI gaps.**
   A time interval is classified as affected by STAI when the moving-window `Mc` exceeds the reference value by the selected confidence threshold. The latter is a (user-defined) multiple of the standard deviation estimated via bootstrap: the lower the multiple, the more RESTORE is sensitive to transient variations in `Mc`.

8. **Estimate missing magnitudes.**
   For each STAI step, RESTORE estimates the expected number of earthquakes and related magnitudes using the Gutenberg-Richter (GR) relation. If the Gutenberg-Richter b-value is None, RESTORE uses the default value b = 1.0.

9. **Simulate missing-event properties.**
   Occurrence times, epicentral coordinates, and hypocentral depths are generated stochastically from the observed properties of the filtered catalog.

10. **Build the replenished catalog.**
    The output catalog contains:
    * original events with magnitude greater than or equal to the reference `Mc`;
    * reconstructed events generated inside the STAI region.

    The final column identifies event origin:
    * `flag = 0`: original event;
    * `flag = 1`: reconstructed event.

11. **Generate the requested ensemble.**
    RESTORE runs for all requested combinations of `size`, `step`, and `seed`. Each `(size, step, seed)` combination is treated as one stochastic realization.

12. **Select a representative realization.**
    By inspecting all successful realizations, RESTORE calculates the median value of each representative metric and defines the central ensemble behavior. Each realization is compared with this center using a normalized distance. The realization with the smallest distance is reported as representative.


## Notes for the user

RESTORE produces plausible stochastic realizations of missing seismicity, not reconstructions of the events that actually occurred. Users are therefore encouraged to generate an ensemble of realizations that explores both the stochastic variability and the variability associated with the choice of moving window size and step.
RESTORE uses the distance from the central ensemble behavior as criterion for identifying a representative realization, but users are free to adopt alternative criteria.
