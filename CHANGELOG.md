# Changelog

## 3.0.0

### Added

- Support for multiple RESTORE realizations, with descriptive statistics for individual realizations.
- Automatic selection of a representative realization, defined as the successful run that is closest to the central behavior of the ensemble.
- `settings.txt` file reporting the input parameters for each realization.
- `metadata.json` file reporting the RESTORE version, Python version, input file path and file provenance.
- `restore` command line through `pyproject.toml`.
- Map showing the full input catalog and the spatial region used for STAI analysis.

### Changed

- Input magnitudes are now required to be binned at `ΔM = 0.1`. Any conversion to this magnitude binning is left to the user.
- If no `b-value` is provided, RESTORE now uses the default value `b = 1.0`.
- Magnitude of completeness is now estimated with the Lilliefors test using truncated exponential dithering ([Stallone and Spassiani (2026)](https://doi.org/10.26443/seismica.v5i1.2257)).
- Stochastic reproducibility is now controlled through explicit random seeds and `default_rng` generator.
- Improved messages, documentation, and function docstrings.

### Fixed

- Spatial smoothing now clips numerical values to the valid domain before applying `arccos`.
- Spatial sampling now validates that smoothed rates are finite and non-zero.
- Fitted hypocentral depth density is now validated for finite, positive, and non-zero values before sampling.
- Beta hypocentral depth distribution is now fitted in `[0, 1]` coordinates and then converted back to physical depth units for sampling and plotting.
- Fixed minor bugs.