# Cicci_2027_CommsMed

Core code of a patient-specific ventricular electrophysiology pipeline: from a biventricular and torso mesh to the calibration of the model parameters on a 12-lead ECG.

This code accompanies the manuscript "Personalising ventricular electrophysiology models from CT and ECG for 3D activation imaging in healthy cases" (Communications Medicine, [CITATION AND DOI TO BE ADDED]).

Version: [TO BE ADDED]. Licence: [TO BE ADDED].

## What the pipeline does

1. Prepares the meshes and the universal ventricular coordinates (UVCs).
2. Samples the model parameters (Latin hypercube sampling) and runs the lead-field and electrophysiology simulations.
3. Computes the 12-lead ECGs from the body surface potentials (BSP) and builds a PCA basis.
4. Runs a local sensitivity analysis to choose the parameters to calibrate.
5. Builds a target, either from a simulation or from a clinical 12-lead ECG.
6. Trains Gaussian process emulators (GPEs) and runs Bayesian history matching, wave after wave.
7. Selects the best-fit parameter set (smallest RMSE between emulated and target ECG).

The single-cell model (ToRORd_dynCl) has its own global (Sobol) sensitivity analysis, described in step 0 of `workflow.md`.

## Folder content

| Path | Content |
| --- | --- |
| `config.py` | The two paths you must set (see below). |
| `workflow.md` | Step-by-step instructions, with the settings of each script. |
| `libEP/` | Mesh, sampling, simulation and ECG scripts (numbered in the order of use) and `libEP/utils`. |
| `libEP/data/` | Parameter files, parameter ranges, label files, CARPentry parameter files and tags. |
| `libGPE/` | `gpe_training.py`, `history_matching.py`, `sensitivity_analysis.py` and `libGPE/utils`. |
| `libGPE/copy_into_gperks/` | Two modified GPErks files (see "GPErks"). |

## What is not included

- CARPentry, meshtool and mguvc (see "Requirements").
- GPErks (see "GPErks").
- Patient data: the meshes, the electrode positions, the local activation times and the clinical ECGs. The expected formats are described in "Input data".
- Job scripts for a computing cluster. `libEP/4bis_run_simulations.py` runs the simulations on one machine with `mpiexec`. For a cluster, write your own job script from the `carp.pt` commands in that file.
- Tests.

## Requirements

- Python 3.12.7.
- CARPentry release 2024 (academic licence), with meshtool and mguvc. They are not distributed here. Obtain them from the CARPentry authors.
- GPErks, at upstream commit `dd8b7af` (see "GPErks").
- Python packages: GPErks, SALib, gpytorch, matplotlib, natsort, numpy, pandas, pytictoc, scipy, seaborn, scikit-learn, torch, torchmetrics.
  Versions used: [TO BE ADDED]. NumPy must be below 2 (the pinned GPErks uses `numpy.infty`). `setuptools` must be below 81 (the pinned GPErks imports `pkg_resources`). pandas 2.1 or later is needed (`DataFrame.map`).

## GPErks

GPErks has no version tags on GitHub, and its `setup.py` says 0.1.2 in many commits. The GPErks used for this work is the upstream commit `dd8b7af` (18 September 2024). Install that commit, not the latest one.

Two files are modified versions of GPErks files (MIT licence, notice in each file and in `libGPE/copy_into_gperks/LICENSE_GPErks.txt`):

- `gsa_mod.py` uses SALib (first-order and total-order indices only) and has its own box plot and heat map functions, which save the figures.
- `history_matching_mod.py` plots the inputs in several figures (up to 19 inputs each), can mark a target point, and has a new `plot_wave_diagonal`.

Copy both files into the folder `GPErks/perks/` of your GPErks installation, next to the original files. They do not replace anything.

## Third-party code

The functions in `libGPE/utils/msfun.py` were extracted from the library of Marina Strocchi: <https://github.com/CEMRG-publications/Strocchi_2023_PLoSCompBio>. Please cite:

Strocchi M, Longobardi S, Augustin CM, Gsell MAF, Petras A, Rinaldi CA, Vigmond EJ, Plank G, Oates CJ, Wilkinson RD, Niederer SA. Cell to whole organ global sensitivity analysis on a four-chamber heart electromechanics model using Gaussian processes emulators. PLoS Comput Biol. 2023;19(6):e1011257. [doi:10.1371/journal.pcbi.1011257](https://doi.org/10.1371/journal.pcbi.1011257)

## Setup

1. Install the requirements above.
2. Open `config.py` and set
   - `CARP_INSTALL_DIR`: the folder where CARPentry is installed (it must contain `bin/carp.pt`, `meshtool/meshtool` and `carputils`);
   - `OUTPUT_DIR`: the folder where all results are written.
   Each script stops with a clear message if a path is not set. Use absolute paths. A relative path is taken relative to the folder of `config.py`, not to the folder you launch the script from.
3. All other settings are variables at the top of `main()` in each script. Edit them there, as described in `workflow.md`.

The scripts write all results under `OUTPUT_DIR`.

## Input data

Patient-derived data are not released. You must supply:

1. **Torso mesh** (in `carp_txt` format). The tag numbers of the regions, and the settings used to segment the myocardium, are set in `libEP_DIR + "/data/tag_HEART.json"` (the bath tags are in `libEP_DIR + "/data/tag_BATH.json"`).
2. **Electrodes** (in `pts` format):
   - First line: the number of points (10).
   - Then one line per electrode (V1, V2, V3, V4, V5, V6, RA, LA, RL, LL) with x, y and z, in the same units as the mesh.
3. **Clinical ECG** (in `dat` format):
   - 13 columns: time, LI, LII, LIII, aVR, aVL, aVF, V1, V2, V3, V4, V5, V6.
   - Time 0 is the onset of the QRS complex.
   - The sampling interval is 2 ms.
   - The number of samples must be the same as in the simulated ECGs, otherwise the projection on the PCA basis fails.
4. **[OPTIONAL]** Local activation times (in `dat` format): one value per node of the biventricular mesh. It is given with the argument `prepacing_lats` of `setup.write_par_tissue` for the pre-pacing of the cells. Without it, no pre-pacing is applied.

The location of each file is given in `workflow.md`.

## How to run

Follow `workflow.md`. In short:

| Step | Scripts |
| --- | --- |
| 0. Single-cell model and its sensitivity analysis | `libEP/0_single_cell.py`, `libGPE/gpe_training.py` (with `single_cell = True`), `libGPE/sensitivity_analysis.py` |
| 1. Meshes and UVCs | `libEP/1_preprocessing.py` |
| 2. Sampling, simulations, PCA basis | `libEP/2_sampling.py`, `3_sample_assembling_pool.py`, `4_simulation_assembling.py`, `4bis_run_simulations.py`, `5_extract_outputs_synthetic.py` |
| 3. Local sensitivity analysis | `libEP/linear_fitting_sampling.py`, `3_sample_assembling_pool.py`, `4_simulation_assembling.py`, `4bis_run_simulations.py`, `5_extract_outputs_synthetic.py`, then `linear_fitting_screening.py` |
| 4.1 Target from a simulation | `libEP/2_sampling_calibration.py`, `3_sample_assembling_pool.py`, `4_simulation_assembling.py`, `4bis_run_simulations.py`, `5_extract_outputs_synthetic.py`, then `libGPE/gpe_training.py` (with `train_gpes = False`) |
| 4.2 Target from a clinical ECG | `libEP/5_extract_outputs_clinical.py`, then `libGPE/gpe_training.py` (with `train_gpes = False`) |
| 5. QRS-complex calibration | `libEP/2_sampling_calibration.py`, `3_sample_assembling_pool.py`, `4_simulation_assembling.py`, `4bis_run_simulations.py`, `5_extract_outputs_synthetic.py`, `libGPE/gpe_training.py`, `libGPE/history_matching.py`, `libEP/2_sampling_simul_bhm.py` (repeated for each wave), then `libEP/6_select_best_fit.py` |
| 6. T-wave calibration | `libEP/6_select_best_fit.py` (on the result of step 5), then the scripts of step 5 |

Notes:

- The default sample folder differs between stages (for example `samples_LHS`, `samples_fitting` and `samples_initial_space`). Set `samplefolder` to the same value in all the scripts of a stage, as in `workflow.md`.
- The input-parameter label files must be in Python `sorted()` order (capital letters before lower case). The provided files already are. The scripts stop with an error if a label file has empty lines, spaces or duplicates.
- CARPentry writes `lead_field_egms.dat` (the BSP at the 10 electrodes) and `ECG_lead_field` (13 columns: time and the 12 leads). These names come from CARPentry and are kept.
- `libGPE/history_matching.py` writes plots with the prefix `subset1_` (for example `subset1_impl_<wave>.png`), because the inputs are split into subsets.

## Contact

Ludovica Cicci, l.cicci@imperial.ac.uk. Cardiac Electromechanics Research Group (CEMRG), Imperial College London.
