# RTP–AFNI ROI optimization


This repository contains code supporting analyses of correspondence between
real-time fMRI (RTP) signals and post-processed AFNI signals during real-time
fMRI neurofeedback.

The repository provides documented implementations of the central analysis
workflow, including candidate-ROI generation, static and dynamic ROI
optimization, procedural-null simulations, and participant-level statistical
analyses.

The public scripts focus on the scientific analysis procedures and key
parameters. 


## Repository contents

### `01_generate_candidate_rois.sh`

Generates the candidate spherical ROIs used for static and dynamic optimization.

For each neurofeedback target, 810 candidate 7-mm-radius spherical ROIs are
generated on a 1-mm grid (9 × 10 × 9 locations).

Candidate centers in the AFNI coordinate convention used by the analysis:

**Amygdala**
- X = 17–25
- Y = 1–10
- Z = −20 to −12
- Original AFNI ROI center = [21, 5, −16]

**Parietal (HIPS)**
- X = 38–46
- Y = 44–53
- Z = 44–52
- Original AFNI ROI center = [42, 48, 48]



### `02_roi_optimization.m`

Implements the Original, Static, and Dynamic ROI analyses.

**Original ROI**

The mean AFNI time series is extracted from the prespecified 7-mm spherical ROI
and correlated with the corresponding RTP time series.

**Static optimization**

All 810 candidate ROIs are evaluated, and the candidate producing the highest
RTP–AFNI Pearson correlation is selected.

For the full-run analysis, selection is based on the complete 260-TR run.

For condition-specific analyses, Static optimization is performed separately
using the concatenated Happy and Rest time points.

**Dynamic optimization**

For the full-run analysis:

- window length = 20 TRs;
- windows advance by one TR;
- 241 overlapping windows are evaluated in each 260-TR run;
- the candidate ROI with the highest RTP–AFNI correlation is selected within
  each window;
- overlapping winning ROI time series are averaged to reconstruct a
  260-TR dynamically optimized AFNI signal.

For condition-specific analyses, Dynamic optimization is performed separately
within each actual 20-TR Happy or Rest block so that optimization does not cross
condition boundaries.


### `03_primary_statistical_analysis.R`

Implements the participant-level statistical analyses.

Correlations are processed hierarchically:

1. Pearson correlations are calculated separately for each available run.
2. Run-level correlations are Fisher Z-transformed.
3. Available runs are averaged within each visit.
4. Visit-level values are averaged within each participant.
5. Participant-level Fisher-Z values are used for inferential analyses.
6. Back-transformed Pearson-r values are used for descriptive reporting and
   visualization.

The statistical models evaluate RTP–AFNI correspondence across ROI strategy
(Original, Static, Dynamic), condition (Full Run, Happy, Rest), and
neurofeedback target.


### `04_procedural_null_simulation.m`

Generates the procedural-null distributions used to estimate the magnitude of
correlation gain that can arise from applying the ROI-search procedure to noisy
signals containing shared variance.

Key simulation parameters:

- 260 TRs per simulated run;
- shared latent signal;
- simulated RTP signal = 0.6 × shared signal + independent Gaussian noise;
- AFNI shared-signal scaling = 0.7;
- voxel-wise Poisson spatial weights with lambda = 0.2;
- enhanced shared-signal weighting within the target ROI using Poisson
  boosting with lambda = 2;
- weights normalized to a maximum of 1;
- 4-mm FWHM spatial Gaussian smoothing;
- 2-mm isotropic simulated voxel resolution;
- 810 candidate ROIs per target;
- 20-TR Dynamic optimization windows.

One hundred independent simulations are generated per target for the full-run
analysis.

A separate set of 100 simulations per target is generated for the
condition-specific analysis, with each simulation providing both Happy and Rest
estimates.

The simulation is intended as a procedural null for ROI-search-related
correlation inflation. 


### `05_procedural_null_analysis.R`

Compares empirical optimization gains with the corresponding target- and
condition-specific procedural-null distributions.

Both empirical and simulated correlations are Fisher Z-transformed before
optimization gains are calculated.

Static optimization gain:

`Z(Static) - Z(Original)`

Dynamic optimization gain:

`Z(Dynamic) - Z(Static)`

For each target and condition, the analysis reports the number and proportion
of the 100 simulated gains that are greater than or equal to the observed
empirical mean gain.

These exceedance counts are treated as descriptive procedural benchmarks and
are not interpreted as simulation-based p-values.

Raw Pearson-r gains are retained for descriptive interpretation and
visualization.


## Task structure

Each analyzed run contains 260 TRs after removal of the initial stabilization
volumes.

The task consists of repeating 20-TR blocks of:

1. Rest
2. Happy
3. Count

For condition-specific analyses:

- Happy blocks: 2, 5, 8, and 11 (80 TRs total)
- Rest blocks: 1, 4, 7, 10, and 13 (100 TRs total)

Count blocks remain included in the full-run analysis but are not included in
the Happy- or Rest-specific analyses.


## Software versions

The analyses were conducted using:

- Turbo BrainVoyager 4.4
- AFNI 18.0.22
- MATLAB R2023a
- R 4.5.0

Principal R packages:

- tidyverse 2.0.0
- afex 1.4.1
- emmeans 1.11.0
- R.matlab 3.7.0


## Reproducibility

The scripts provided here are streamlined, documented implementations of the
scientific analysis workflow used in the manuscript.
