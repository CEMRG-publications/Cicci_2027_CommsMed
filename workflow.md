# Workflow - Cicci_2027_CommsMed

Status: all steps are updated and checked against the scripts. Nothing has been run end to end.
Parts marked **[TO BE UPDATED]** wait for the requirements file.

## General rules

- Paths are set once in `config.py`: `CARP_INSTALL_DIR` (CARPentry folder) and `OUTPUT_DIR` (folder for all results).
- All other settings are variables at the top of `main()` in each script.
- Examples below use `basefolder = OUTPUT_DIR + "/<meshname>"`.
- The pre-processing (step 1) creates `<meshname>_subendolayer`, `<meshname>_FEC` and `<meshname>_gregion` from the base mesh `<meshname>`.
- In `3_sample_assembling_pool.py`, `meshname` is the base mesh name `<meshname>` (for example `KCL1`), never a derived mesh. The script adds `_subendolayer` and `_gregion` by itself. In `4_simulation_assembling.py` and `4bis_run_simulations.py`, `meshname` is the full name of the mesh used in the simulations (for example `KCL1_gregion`).

## Requirements

- CARPentry release 2024 (academic licence), including meshtool and mguvc.
- Python 3.12.7 and the packages listed in the README **[TO BE UPDATED: requirements file]**.
- Torso mesh in `carp_txt` format, with the tag numbers set in `libEP/data/tag_HEART.json`.
- Electrode pts-file(s) should be created ad-hoc with the following format:
    - first line: number of points 
        - for the 12-lead ECG: 10
    - second line onwards: one line per electrode with x, y, z (same units as the mesh), in the correct order;
        - for a 12-lead ECG: V1, V2, V3, V4, V5, V6, RA, LA, RL, LL
    - NOTE: second-to-last electrode used as ground.
- Optional: a file with local activation times (one value per node of the biventricular mesh) for the pre-pacing of the cells, given with the argument `prepacing_lats` of `setup.write_par_tissue`. Without it, no pre-pacing is applied.

## Part A - 0D: single-cell model

The single-cell model (e.g., ToRORd_dynCl) and its global sensitivity analysis. This part does not depend on Part B.

### 0. Single-cell model and its global sensitivity analysis

This step is independent of the steps below.

- **a.** Run `0_single_cell.py` (set `num_proc`, the number of samples simulated at the same time). It samples the ionic model parameters (in `<basefolder>/samples`), runs the single-cell simulations and writes the outputs in `<basefolder>/qoi`.
- **b.** Run `libGPE/gpe_training.py` with `single_cell = True`. It trains the GPEs in `<basefolder>/gpe/GPE_Matern_GL`.
- **c.** Run `libGPE/sensitivity_analysis.py`. It reads those GPEs and writes the Sobol indices and the figures in `<basefolder>/gsa`.

---

## Part B - 3D: biventricular and torso model

The meshes, the simulations of the biventricular and torso model, the ECGs and the calibration.

### 1. Pre-process meshes and generate UVCs

- **a.** Set correct mesh tags in `libEP/data/tag_HEART.json` and `libEP/data/tag_BATH.json`
- **b.** Run `1_preprocessing.py` for current mesh `<meshname>`
    - the input mesh (in `carp_txt` format) is read from `OUTPUT_DIR + "/mesh/<meshname>"`
    - all outputs will be saved in `OUTPUT_DIR + "/" + <meshname>`

### 2. Compute PCA basis/bases

- **a.** Run `2_sampling.py` with
    - `paramfilename = "simulation_parameters"`
    - `labelfilename_scalar = "xlabels_scalar"`
    - `labelfilename_gradient = "xlabels_gradient"`
    - `rangefilename_scalar = "ranges_scalar"`
    - `rangefilename_gradient = "ranges_gradient"`
    - `N_samples = 10 * ( len(xlabels_scalar) + 3*len(xlabels_gradient) )`
    - `samplefolder = basefolder + "/samples_LHS"`
- **b.** Optional: a preliminary single-cell check of the samples, to discard samples that are highly likely to fail (APD, VmMax or VmMin outside a viable range). It needs CARPentry (`bench.pt`). Set `check_single_cell = True` (and `num_proc`) in `2_sampling.py` in step **a**. If the sample folder already exists, the script skips the sampling and runs only the check. The simulations already completed are not run again. The accepted samples are listed in `accepted_samples.txt` (the later steps skip the others), and the samples with a missing output in `missing_statistics.txt`. To run the check again, remove `accepted_samples.txt` first.
- **c.** Run `3_sample_assembling_pool.py` with
    - `labelfilename_scalar = "xlabels_scalar"`
    - `samplefolder = basefolder + "/samples_LHS"`
    - `SElayer = True`
    - `heterogeneity = True`
    - `labelfilename_im = "ToRORd_dynCl"`
- **d.** Run `4_simulation_assembling.py` with
    - `ionicmodelname = "ToRORd_dynCl"`
    - `labelfilename_im = "ToRORd_dynCl"`
    - `labelfilename_scalar = "xlabels_scalar"`
    - `samplefolder = basefolder + "/samples_LHS"`
    - `propagation = "RE+"`
    - `clusterfolder = basefolder` for a local run, otherwise the absolute path of the folder on the cluster
- **e.** Move or copy `basics/mesh` and `basics/data` (saved by step **d** in `basefolder + "/basics"`) and the samples folder of the current step (`<samplefolder>`) into `clusterfolder`. Run the simulations (lead field first, then EP) for every sample:
    - locally: run `4bis_run_simulations.py` with 
        - `samplefolder = basefolder + "/samples_LHS"`
        - `simulationfolder = clusterfolder + "/simulations_LHS"`
    - on a cluster: `clusterfolder` (step **d**) is the path on the cluster. Write your own job script with the `carp.pt` commands used in `4bis_run_simulations.py`. The folder with the samples keeps its name on the cluster.
- **NOTE:** simulations whose `lead_field_egms.dat` is missing or empty are treated as "failed" and skipped by step **f**.
- **f.** Run `5_extract_outputs_synthetic.py` with
    - `paramfilename = "simulation_parameters"`
    - `labelfilename_scalar = "xlabels_scalar"`
    - `labelfilename_gradient = "xlabels_gradient"`
    - `compute_PCA = True`
    - for the QRS complex:
        - `ECGsegment = "QRScomplex"` (sets `qoi_list = ["12leads_truncated"]`)
        - `qoifolder = basefolder + "/PCAbasis_<meshname>_150ms_99%"`
        - `tend_truncated = 150` ms
    - for the T-wave:
        - `ECGsegment = "Twave"` (sets `qoi_list = ["12leads_REP"]`)
        - `qoifolder = basefolder + "/PCAbasis_<meshname>_REP_99%"`
        - `tend = 500` ms

  The PCA basis is saved in `qoifolder`, which is named so that the later steps find it (with `compute_PCA = False`).

### 3. Perform Local Sensitivity Analysis (LSA)

- **a.** Run `linear_fitting_sampling.py` with
    - i. `paramfilename = "simulation_parameters_fitting_biv"`
    - ii. `labelfilename = "xlabels_fitting_biv"`
    - iii. `rangefilename = "ranges_fitting_biv"`
    - iv. `samplefolder = basefolder + "/samples_fitting"`
- **b.** Run `3_sample_assembling_pool.py` with
    - i. `labelfilename_scalar = "xlabels_fitting_biv"`
    - ii. `labelfilename_im = "ToRORd_dynCl"`
    - iii. `samplefolder = basefolder + "/samples_fitting"`
    - iv. `heterogeneity = True`
- **c.** Check that ALL samples have been accepted. If not, find out why and modify the parameters of the rejected ones accordingly. At the end, all samples MUST have been accepted.
- **d.** Run `4_simulation_assembling.py` with
    - i. `labelfilename_scalar = "xlabels_fitting_biv"`
    - ii. `samplefolder = basefolder + "/samples_fitting"`
    - iii. `propagation = "RE+"`
    - iv. `lonfolder_ref = ""`
- **e.** Run the simulations and download the results (as in steps 2e and 2f), with
    - i. `samplefolder = basefolder + "/samples_fitting"` in `4bis_run_simulations.py`
    - ii. `simulationfolder = basefolder + "/simulations_fitting"` in `4bis_run_simulations.py`, so that the results of step 2 are not overwritten (simulation numbers start from 0 in every stage)
- **f.** Run `5_extract_outputs_synthetic.py` with
    - i. `paramfilename = "simulation_parameters_fitting_biv"`, `labelfilename_scalar = "xlabels_fitting_biv"`, `labelfilename_gradient = ""`, `samplefolder = basefolder + "/samples_fitting"` and `simulationfolder = basefolder + "/simulations_fitting"`
    - ii. for the QRS complex:
        1. `PCAbasisfolder = basefolder + "/PCAbasis_KCL1_150ms_99%"`
        2. `qoifolder = basefolder + "/qoi_fitting_150ms"`
        3. `ECGsegment = "QRScomplex"`, `tend_truncated = 150` ms
    - iii. for the T-wave:
        1. `PCAbasisfolder = basefolder + "/PCAbasis_KCL1_REP_99%"`
        2. `qoifolder = basefolder + "/qoi_fitting_REP"`
        3. `ECGsegment = "Twave"`, `tend = 500` ms
- **g.** Run `linear_fitting_screening.py` with
    - i. `cutoff = 0.95`
    - ii. `method = "uniform_variance"`
    - iii. for the QRS complex:
        1. `ECGname = "ECG_12leads_truncated"`
        2. `WLname = "WL_12leads_truncated"`
        3. `PCAbasisfolder = basefolder + "/PCAbasis_KCL1_150ms_99%"`
        4. `qoifolder = basefolder + "/qoi_fitting_150ms"`
    - iv. for the T-wave:
        1. `ECGname = "ECG_12leads_REP"`
        2. `WLname = "WL_12leads_REP"`
        3. `PCAbasisfolder = basefolder + "/PCAbasis_KCL1_REP_99%"`
        4. `qoifolder = basefolder + "/qoi_fitting_REP"`

### 4. Generate a target from a simulation

- **a.** Run `linear_fitting_sampling.py` with
    - i. `paramfilename = "simulation_parameters_fitting_biv"`
    - ii. `labelfilename = "xlabels_fitting_biv"`
    - iii. `rangefilename = "ranges_fitting_biv"`
    - iv. `samplefolder = basefolder + "/target/samples"`

  The number of samples is 2 * the number of parameters (the minimum and the maximum of each parameter). It is set in the script.
- **b.** Run `3_sample_assembling_pool.py` and `4_simulation_assembling.py` (`propagation = "RE+"`), with `samplefolder = basefolder + "/target/samples"`.
- **c.** Run the simulations and download the results, with `samplefolder = basefolder + "/target/samples"` and `simulationfolder = basefolder + "/target/simulations_candidates"` in `4bis_run_simulations.py`.
- **d.** Choose the "best" ECG and rename it to be sample/simulation(s) "0" in the target folder, that is
    - i. `target/samples/0`
    - ii. `target/simulations/SIM_LF_0` and `target/simulations/SIM_EP_0` (taken from `simulations_candidates`)
- **e.** Define the label files of the parameters to be calibrated (one label per line, in alphabetical order):
    - i. `xlabels_calibration_QRS.txt`: parameters selected from the screening results in `qoi_fitting_150ms` (see `final_ranking_cumvar` and the cut-off figure)
    - ii. `xlabels_calibration_Twave.txt`: parameters selected in the same way from `qoi_fitting_REP`, without the parameters already in the QRS-complex list

  The files used for the paper are provided in `libEP/data/xlabels`. To calibrate a different selection, replace them with your own files of the same format. Copy both files into `basefolder + "/target/samples"`.
- **f.** Run `5_extract_outputs_synthetic.py` with `paramfilename = "simulation_parameters_fitting_biv"`, `labelfilename_gradient = ""`, `samplefolder = basefolder + "/target/samples"` and `simulationfolder = basefolder + "/target/simulations"`, and
    - i. for the QRS complex:
        1. `labelfilename_scalar = "xlabels_calibration_QRS"`
        2. `PCAbasisfolder = basefolder + "/PCAbasis_KCL1_150ms_99%"`
        3. `qoifolder = basefolder + "/target/qoi_QRScomplex"`
        4. `ECGsegment = "QRScomplex"`, `tend_truncated = 150`
    - ii. for the T-wave:
        1. `labelfilename_scalar = "xlabels_calibration_Twave"`
        2. `PCAbasisfolder = basefolder + "/PCAbasis_KCL1_REP_99%"`
        3. `qoifolder = basefolder + "/target/qoi_Twave"`
        4. `ECGsegment = "Twave"`, `tend = 500`

  The columns of the matrix X written by this script follow the alphabetical order of the labels. Since the provided label files are already in alphabetical order, no reordering of the labels is needed.
- **g.** For each qoi (`ECGsegment = "QRScomplex"` and `"Twave"`), run `libGPE/gpe_training.py` with `train_gpes = False` to construct the data folder and `outputs.csv` directly in the target qoi folder (`OUTPUT_DIR + "/KCL1/target/qoi_<segment>"`, set by `target_folder` in the script). No GPE is trained and no copying is needed.

### 5. Generate a target from a clinical 12-lead ECG

- **a.** Prepare the clinical ECG as the text file `ECG_from_BSP.dat` with 13 columns: time, LI, LII, LIII, aVR, aVL, aVF, V1, V2, V3, V4, V5, V6. Time 0 must be the onset of the QRS complex. The sampling interval must be 2 ms, and the number of samples must equal that of the simulated ECGs (otherwise the projection on the PCA basis fails).
- **b.** Copy `ECG_from_BSP.dat` into `basefolder + "/target"`.
- **c.** Run `5_extract_outputs_clinical.py` with
    - i. for the QRS complex:
        1. `ECGsegment = "QRScomplex"` (`PCAbasisfolder = basefolder + "/PCAbasis_KCL1_150ms_99%"`)
        2. `qoifolder = basefolder + "/target/qoi_QRScomplex"`
    - ii. for the T-wave:
        1. `ECGsegment = "Twave"` (`PCAbasisfolder = basefolder + "/PCAbasis_KCL1_REP_99%"`)
        2. `qoifolder = basefolder + "/target/qoi_Twave"`
- **d.** Use the label files `xlabels_calibration_QRS.txt` and `xlabels_calibration_Twave.txt` (see step 4e) and copy them into `basefolder + "/target"`.
- **e.** For each qoi, run `libGPE/gpe_training.py` with `train_gpes = False` to construct the data folder and `outputs.csv` directly in the target qoi folder (`OUTPUT_DIR + "/KCL1/target/qoi_<segment>"`, set by `target_folder` in the script). No GPE is trained and no copying is needed. A clinical target has no X folder: the script skips the inputs, so no change of the code is needed.
- **f.** The labels must be in alphabetical order, as the columns of X in step 4f. The provided label files already are. (`reorder_labels.py`, used for this in the old workflow, is not part of this release.)
- **g.** Copy the label file `xlabels_calibration_QRS.txt` (or `xlabels_calibration_Twave.txt`) to the data folder of the corresponding target qoi folder, and change the name into `xlabels.txt`. `history_matching.py` and `6_select_best_fit.py` read it from there.

### Pipeline: QRS-complex calibration

Here: `<history_matching>` = `basefolder + "/History_Matching/QRScomplex"`.

In all the steps below, set `clusterfolder = <history_matching>`. The simulation results of each stage go in their own folder, so that stages do not overwrite each other.

1. Run `2_sampling_calibration.py` with
    - `hmfolder = <history_matching>`
    - `paramfilename = "simulation_parameters_fitting_biv"`
    - `labelfilename = "xlabels_calibration_QRS"`
    - `rangefilename = "ranges_fitting_biv"`
    - `N_samples = len(xlabels) * 20` (already set in the script)
    - `samplefolder = <history_matching> + "/samples_initial_space"`

   The script copies the label files, including `xlabels_empty.txt` and `ToRORd_dynCl.txt`, into the sample folder.

   NOTE: `xlabels_calibration_QRS` must be defined based on the sensitivity analysis results.
2. Run `3_sample_assembling_pool.py` with
    - `labelfilename_scalar = "xlabels_calibration_QRS"`
    - `labelfilename_im = "xlabels_empty"`
    - `samplefolder = <history_matching> + "/samples_initial_space"`
    - `SElayer = True`
    - `heterogeneity = True`
3. Run `4_simulation_assembling.py` with
    - `ionicmodelname = "TT2"` (for the first waves; use `ToRORd_dynCl` in the last waves)
    - `labelfilename_im = "xlabels_empty"`
    - `labelfilename_scalar = "xlabels_calibration_QRS"`
    - `samplefolder = <history_matching> + "/samples_initial_space"`
    - `propagation = "RE-"`
    - `hps_ref = False`
    - `cvfolder_ref = ""`
    - `lonfolder_ref = ""`
    - `clusterfolder = <history_matching>`
4. Move or copy `basics/mesh` and `basics/data` (saved by step 3 in `basefolder + "/basics"`) into `<history_matching>`. Then change `<history_matching> + "/data/pfile_BASE.par"` as
    - `dt = 25`
    - `tend = 150`
5. Locally, run `4bis_run_simulations.py` with (on a cluster, see step 2e of Part B)
    - `clusterfolder = <history_matching>`
    - `samplefolder = <history_matching> + "/samples_initial_space"`
    - `simulationfolder = <history_matching> + "/simulations_initial_space"`
6. Run `5_extract_outputs_synthetic.py` with
    - `ECGsegment = "QRScomplex"` (`tend_truncated = 150`)
    - `paramfilename = "simulation_parameters_fitting_biv"`
    - `labelfilename_scalar = "xlabels_calibration_QRS"`
    - `labelfilename_gradient = ""`
    - `compute_PCA = False` (the PCA basis computed in step 2f of Part B is read from `basefolder + "/PCAbasis_<meshname>_150ms_99%"`)
    - `samplefolder = <history_matching> + "/samples_initial_space"`
    - `simulationfolder = <history_matching> + "/simulations_initial_space"`
    - `qoifolder = <history_matching> + "/qoi_initial_space"`
7. Train the GPEs with `libGPE/gpe_training.py`, with
    - `ECGsegment = "QRScomplex"`
    - `waveno = 0` for the first training (initial space only). For the later waves, `waveno` = the last wave included in the training data. The script uses the data of `qoi_initial_space` and `qoi_wave_1` to `qoi_wave_<waveno>`, and writes the GPEs in `<history_matching> + "/gpe/initial_space"` or `<history_matching> + "/gpe/wave_<waveno>"` (in a subfolder `GPE_Matern_GL`).
8. Perform the history matching with `libGPE/history_matching.py` (`ECGsegment = "QRScomplex"`, `waveno = 1` for the first wave). Wave n uses the GPEs of wave n-1 (`initial_space` for wave 1). It writes `X_simul_<n>.txt`, `X_test_<n>.txt` and `wave_<n>.json` in `<history_matching> + "/output"`, and the figures in `<history_matching> + "/figures"`.
9. Run `2_sampling_simul_bhm.py` with
    - `hmfolder = <history_matching>`
    - `waveno = n`
    
    It writes the new samples in `<history_matching> + "/samples_wave_<n>"` and copies the label files into that folder.

10. Repeat steps 2 to 9 for each wave n, with
    - `samplefolder = <history_matching> + "/samples_wave_<n>"` (steps 2, 3, 5 and 6)
    - `simulationfolder = <history_matching> + "/simulations_wave_<n>"`
    - `qoifolder = <history_matching> + "/qoi_wave_<n>"`
    - no new copy of `basics` and no new change of `pfile_BASE.par` (step 4): `mesh` and `data` are already in `<history_matching>`

    and then `gpe_training.py` with `waveno = n` and `history_matching.py` with `waveno = n + 1`.

Notes:

- Remember that `gpe_training.py` adds the new qoi folder when `waveno` is increased.
- For the last wave of the history matching, use `n_simuls_factor = 100` in `history_matching.py` (the number of new samples is then `len(xlabels) * 100`).

When the history matching is finished, select the best-fit point with `6_select_best_fit.py` (see step 1 of the T-wave calibration).

### Pipeline: T-wave calibration

Here: `<history_matching>` = `basefolder + "/History_Matching/Twave"`.

In all the steps below, set `clusterfolder = <history_matching>`. The simulation results of each stage go in their own folder, so that stages do not overwrite each other.

1. Create the folder `<history_matching>`.
2. Define the parameter file `data/simulation_parameters_fitting_biv_best-fit_QRS.json`: a copy of `data/simulation_parameters_fitting_biv.json` in which the values of the parameters listed in `xlabels_calibration_QRS.txt` are replaced by their best-fit values from the QRS-complex calibration (a point chosen among the non-implausible points of the last history matching wave). It is easier to use one of the parameter files from the QRS-complex waves and just change the last few lines.

   Select the best-fit point (smallest RMSE between the emulated and the target ECG) with `6_select_best_fit.py`, with
    - `ECGsegment = "QRScomplex"`
    - `final_waveno` = the last wave of the history matching
    - `K = 100`

   It writes `best_NIMP_wave_<n>.json` in the output folder of the QRS-complex history matching.
3. Run `2_sampling_calibration.py` with
    - `hmfolder = <history_matching>`
    - `paramfilename = "simulation_parameters_fitting_biv_best-fit_QRS"`
    - `labelfilename = "xlabels_calibration_Twave"` (the file is provided in `libEP/data/xlabels`)
    - `rangefilename = "ranges_fitting_biv"`
    - `N_samples = len(xlabels) * 20` (already set in the script)
    - `samplefolder = <history_matching> + "/samples_initial_space"`
    - `check_single_cell = True` (single-cell check, see step 4) and `num_proc` (number of samples simulated at the same time)

   The script copies the label files, including `xlabels_empty.txt` and `ToRORd_dynCl.txt`, into the sample folder.
4. The single-cell check of the samples is run in step 3 with `check_single_cell = True` (T-wave calibration only). Look at `samples_initial_space/missing_statistics.txt` (if any). To resume after a failure, remove `accepted_samples.txt` and run `2_sampling_calibration.py` again with the same settings: the sampling is skipped and the simulations already completed are not run again.

5. **[OPTIONAL]** Define a sample folder for the calibrated QRS-complex parameters, to compute HPS, conduction velocities, fibres and lead field only once, and save it to `<history_matching> + "/ref_QRS"`.
    - Make the directory `<history_matching> + "/ref_QRS/0"`
    - Copy `simulation_parameters_fitting_biv_best-fit_QRS.json` to `<history_matching> + "/ref_QRS/0/param.json"`
    - Copy `xlabels_calibration_QRS.txt` (needed to compute fibres, CVs and HPS), `xlabels_empty.txt` and `ToRORd_dynCl.txt` to `<history_matching> + "/ref_QRS"`
    - Run `3_sample_assembling_pool.py` with
        - `labelfilename_scalar = "xlabels_calibration_QRS"`
        - `labelfilename_im = "xlabels_empty"`
        - `samplefolder = <history_matching> + "/ref_QRS"`
        - `SElayer = True`
        - `heterogeneity = True`
    - Run `4_simulation_assembling.py` with
        - `ionicmodelname = "ToRORd_dynCl"`
        - `labelfilename_im = "xlabels_empty"`
        - `labelfilename_scalar = "xlabels_calibration_QRS"`
        - `samplefolder = <history_matching> + "/ref_QRS"`
        - `propagation = "RE-"`
        - `hps_ref = False`
        - `cvfolder_ref = ""`
        - `lonfolder_ref = ""`
        - `clusterfolder = <history_matching>`
    - Move or copy `basics/mesh` and `basics/data` (saved by the run of `4_simulation_assembling.py` above in `basefolder + "/basics"`) into `<history_matching>`.
    - Compute the lead field of the reference folder: run `4bis_run_simulations.py` locally (on a cluster, see step 2e of Part B) with
        - `clusterfolder = <history_matching>`
        - `samplefolder = <history_matching> + "/ref_QRS"`
        - `simulationfolder = <history_matching> + "/simulations_ref_QRS"`

        The folder has only the sample 0. The EP simulation is also run and can be ignored. 
    - Rename `<history_matching> + "/simulations_ref_QRS/SIM_LF_0"` as `<history_matching> + "/ref_QRS/SIM_LF_ref"`.
6. Run `3_sample_assembling_pool.py` with
    - `labelfilename_scalar = "xlabels_calibration_Twave"`
    - `labelfilename_im = "ToRORd_dynCl"`
    - `samplefolder = <history_matching> + "/samples_initial_space"`
    - `SElayer = True`
    - `heterogeneity = True`

   Compute ONLY the ADJ files (HPS, CVs and fibres can be computed only once and then referenced to `ref_QRS`).

7. Run `4_simulation_assembling.py` with
    - `ionicmodelname = "ToRORd_dynCl"`
    - `labelfilename_im = "ToRORd_dynCl"`
    - `labelfilename_scalar = "xlabels_calibration_Twave"`
    - `samplefolder = <history_matching> + "/samples_initial_space"`
    - `propagation = "RE+"`
    - `hps_ref = True` (if that is the case)
    - `cvfolder_ref = <history_matching> + "/ref_QRS/0/CV"` (if that is the case)
    - `lonfolder_ref = <history_matching> + "/ref_QRS/0"` (if that is the case)
    - `clusterfolder = <history_matching>`

   Move or copy `basics/mesh` and `basics/data` (saved by this step in `basefolder + "/basics"`) into `<history_matching>`, unless this was already done in step 5. Then check that `<history_matching> + "/data/pfile_BASE.par"` has `dt = 5` and `tend = 500`.

8. Locally, run `4bis_run_simulations.py` with (on a cluster, see step 2e of Part B)
    - `clusterfolder = <history_matching>`
    - `samplefolder = <history_matching> + "/samples_initial_space"`
    - `simulationfolder = <history_matching> + "/simulations_initial_space"`
9. Run `5_extract_outputs_synthetic.py` with
    - `ECGsegment = "Twave"` (`tend = 500`)
    - `paramfilename = "simulation_parameters_fitting_biv_best-fit_QRS"`
    - `labelfilename_scalar = "xlabels_calibration_Twave"`
    - `labelfilename_gradient = ""`
    - `compute_PCA = False` (the PCA basis computed in step **2f** of Part B is read from `basefolder + "/PCAbasis_<meshname>_REP_99%"`)
    - `samplefolder = <history_matching> + "/samples_initial_space"`
    - `simulationfolder = <history_matching> + "/simulations_initial_space"`
    - `qoifolder = <history_matching> + "/qoi_initial_space"`
10. Train the GPEs with `libGPE/gpe_training.py`, with
    - `ECGsegment = "Twave"`
    - `waveno = 0` for the first training (initial space only). For the later waves, `waveno` = the last wave included in the training data (see step 7 of the QRS-complex calibration).
11. Perform the history matching with `libGPE/history_matching.py` (`ECGsegment = "Twave"`, `waveno = 1` for the first wave). It writes its files in `<history_matching> + "/output"` and `<history_matching> + "/figures"`, as in step 8 of the QRS-complex calibration.
12. Run `2_sampling_simul_bhm.py` with
    - `hmfolder = <history_matching>`
    - `waveno = n`
    - `paramfilename = "simulation_parameters_fitting_biv_best-fit_QRS"`
    - `labelfilename = "xlabels_calibration_Twave"`

    It writes the new samples in `<history_matching> + "/samples_wave_<n>"` and copies the label files into that folder.

13. Repeat steps 6 to 12 for each wave n, with
    - `samplefolder = <history_matching> + "/samples_wave_<n>"` (steps 6, 7, 8 and 9)
    - `simulationfolder = <history_matching> + "/simulations_wave_<n>"`
    - `qoifolder = <history_matching> + "/qoi_wave_<n>"`
    - no new copy of `basics` (step 7): `mesh` and `data` are already in `<history_matching>`

   and then `gpe_training.py` with `waveno = n` and `history_matching.py` with `waveno = n + 1`.

Notes:

- Remember that `gpe_training.py` adds the new qoi folder when `waveno` is increased.
- For the last wave of the history matching, use `n_simuls_factor = 100` in `history_matching.py`.

When the history matching is finished, select the best-fit point with `6_select_best_fit.py` (`ECGsegment = "Twave"`, `final_waveno` = the last wave).
