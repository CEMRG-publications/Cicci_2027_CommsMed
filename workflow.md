# Workflow - Cicci_2027_CommsMed

Status: all steps are updated and checked against the scripts. Nothing has been run end to end.
Parts marked **[TO BE UPDATED]** wait for a decision (`2bis_check_imp.py` with `run_bench_pool.py`, requirements file).

## General rules

- Paths are set once in `config.py`: `CARP_INSTALL_DIR` (CARPentry folder) and `OUTPUT_DIR` (folder for all results).
- All other settings are variables at the top of `main()` in each script.
- Examples below use `basefolder = OUTPUT_DIR + "/KCL"` and `meshname = "KCL1"`.
- **Paths are being checked script by script.** The paths in this workflow are those of the scripts tested so far: `0_single_cell.py`, `libGPE/gpe_training.py` and `libGPE/sensitivity_analysis.py`. The paths of the other scripts have not been checked yet.
- Each script that reads samples must use the same `samplefolder` as the script that wrote them (for example `basefolder + "/samples_LHS"`).

## Requirements

- CARPentry release 2024 (academic licence), including meshtool and mguvc.
- Python 3.12.7 and the packages listed in the README **[TO BE UPDATED: requirements file]**.
- Torso mesh in `carp_txt` format in `OUTPUT_DIR + "/torsomesh"`, with the tag numbers set in `1_preprocessing.py`.
- Electrode file `OUTPUT_DIR + "/electrodes/electrodes.pts"`:
    - first line: number of points (10); then one line per electrode with x, y, z (same units as the mesh);
    - order of the electrodes: V1, V2, V3, V4, V5, V6, RA, LA, RL, LL.
- Optional: a file with local activation times (one value per node of the biventricular mesh) for the pre-pacing of the cells, given with the argument `prepacing_lats` of `setup.write_par_tissue`. Without it, no pre-pacing is applied.

## 0. Single-cell model (ToRORd_dynCl) and its global sensitivity analysis

This step is independent of the steps below.

- **a.** Run `0_single_cell.py`. It samples the ionic model parameters, runs the single-cell simulations and writes the outputs in `OUTPUT_DIR + "/ToRORd_dynCl/qoi"`.
- **b.** Run `libGPE/gpe_training.py` with `single_cell = True`. It trains the GPEs in `OUTPUT_DIR + "/ToRORd_dynCl/GPE_Matern_GL"`.
- **c.** Run `libGPE/sensitivity_analysis.py`. It reads those GPEs and writes the Sobol indices and the figures in `OUTPUT_DIR + "/ToRORd_dynCl/GSA"`.

## 1. Generate meshes and UVCs

- **a.** Run `1_preprocessing.py`.
- **b.** Convert the mesh to vtk format and check that everything is correct.

## 2. Compute PCA basis/bases

- **a.** Run `2_sampling.py` with
    - i. `paramfilename = "simulation_parameters"`
    - ii. `labelfilename_scalar = "xlabels_scalar"`
    - iii. `labelfilename_gradient = "xlabels_gradient"`
    - iv. `rangefilename_scalar = "ranges_scalar"`
    - v. `rangefilename_gradient = "ranges_gradient"`
    - vi. `N_samples = 10 * ( len(xlabels_scalar) + 3*len(xlabels_gradient) )` (already set in the script)
    - vii. `samplefolder = basefolder + "/samples_LHS"`
- **b.** **[TO BE UPDATED]** Run `2bis_check_imp.py` in combination with `run_bench_pool.py` to check APD, VmMax and VmMin values for single cells.
- **c.** Run `3_sample_assembling_pool.py` (set `samplefolder = basefolder + "/samples_LHS"`).
- **d.** Run `4_simulation_assembling.py` with
    - i. `propagation = "RE+"`
    - ii. `lonfolder_ref = ""`
    - iii. `samplefolder = basefolder + "/samples_LHS"`
    - iv. `clusterfolder = basefolder` for a local run, otherwise the path of the folder on the cluster

  The bath tags are read from `data/tag_BATH.json`. The job scripts are not generated. See step e.
- **e.** Run the simulations (lead field first, then EP) for every sample:
    - locally: run `4bis_run_simulations.py` (set `mainfolder`, `samplefolder = basefolder + "/samples_LHS"`, `N_samples` and `num_proc`). The results are written in `mainfolder + "/simulations"` (variable `simresult`);
    - on a cluster: copy the folders `mesh`, `data` and `samples_LHS` (from `basefolder`) to the cluster, set `clusterfolder` in step d to the path on the cluster, and write your own job script with the `carp.pt` commands used in `4bis_run_simulations.py`. The folder with the samples keeps its name on the cluster.
- **f.** Download the results (the simulations folder), deleting unnecessary files.
- **g.** Failed simulations: step h skips a simulation whose `lead_field_egms.dat` (the BSP file written by CARPentry) is missing or empty.
- **h.** Run `5_extract_outputs_synthetic.py` with `paramfilename = "simulation_parameters"`, `labelfilename_scalar = "xlabels_scalar"` and `labelfilename_gradient = "xlabels_gradient"`, and
    - i. for the QRS complex:
        1. `ECGsegment = "QRScomplex"` (sets `qoi_list = ["12leads_truncated"]`)
        2. `PCAbasisfolder = ""`
        3. `qoifolder = basefolder + "/PCAbasis_KCL1_150ms_99%"`
        4. `tend_truncated = 150` ms
    - ii. for the T-wave:
        1. `ECGsegment = "Twave"` (sets `qoi_list = ["12leads_REP"]`)
        2. `PCAbasisfolder = ""`
        3. `qoifolder = basefolder + "/PCAbasis_KCL1_REP_99%"`
        4. `tend = 500` ms

Removed from the released pipeline: plotting of the ECGs and the screening of unfeasible simulations by total activation time (`sweep_ECG.py`, `accepted_samples`).

## 3. Perform Local Sensitivity Analysis (LSA)

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
    - ii. `simresult = basefolder + "/simulations_fitting"` in `4bis_run_simulations.py`, so that the results of step 2 are not overwritten (simulation numbers start from 0 in every stage)
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

## 4. Generate a target from a simulation

- **a.** Run `linear_fitting_sampling.py` with
    - i. `paramfilename = "simulation_parameters_fitting_biv"`
    - ii. `labelfilename = "xlabels_fitting_biv"`
    - iii. `rangefilename = "ranges_fitting_biv"`
    - iv. `samplefolder = basefolder + "/target/samples"`

  The number of samples is 2 * the number of parameters (the minimum and the maximum of each parameter). It is set in the script.
- **b.** Run `3_sample_assembling_pool.py` and `4_simulation_assembling.py` (`propagation = "RE+"`), with `samplefolder = basefolder + "/target/samples"`.
- **c.** Run the simulations and download the results, with `samplefolder = basefolder + "/target/samples"` and `simresult = basefolder + "/target/simulations_candidates"` in `4bis_run_simulations.py`.
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

## 5. Generate a target from a clinical 12-lead ECG

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

## Pipeline: QRS-complex calibration

Settings checked against the scripts; the pipeline has not been run end to end.

Folders used below: `HM = basefolder + "/History_Matching/QRScomplex"`.
The simulation results of each stage go in their own folder (variable `simresult` in `4bis_run_simulations.py` and variable `simulationfolder` in `5_extract_outputs_synthetic.py`), so that stages do not overwrite each other.

1. Create the folder `HM`.
2. Run `2_sampling_calibration.py` with
    - a. `basefolder = OUTPUT_DIR + "/KCL/History_Matching/QRScomplex"`
    - b. `paramfilename = "simulation_parameters_fitting_biv"`
    - c. `labelfilename = "xlabels_calibration_QRS"` (the file is provided in `libEP/data/xlabels`)
    - d. `rangefilename = "ranges_fitting_biv"`
    - e. `N_samples = len(xlabels) * 20` (already set in the script)
    - f. `samplefolder = HM + "/samples_initial_space"`

   If more samples are needed, run it a second time with `samplefolder = HM + "/samples_initial_space_additional"` and put the two sets together. The script copies the label files, including `xlabels_empty.txt` and `ToRORd_dynCl.txt`, into the sample folder.
3. Run `3_sample_assembling_pool.py` with
    - a. `labelfilename_scalar = "xlabels_calibration_QRS"`
    - b. `labelfilename_im = "xlabels_empty"`
    - c. `samplefolder = HM + "/samples_initial_space"`
    - d. `heterogeneity = True`
4. Run `4_simulation_assembling.py` with
    - a. `ionicmodelname = "TT2"`
    - b. `labelfilename_im = "xlabels_empty"`
    - c. `labelfilename_scalar = "xlabels_calibration_QRS"`
    - d. `samplefolder = HM + "/samples_initial_space"`
    - e. `storagefolder = basefolder`
    - f. `propagation = "RE-"`
    - g. `hps_ref = False`
    - h. `cvfolder_ref = ""`
    - i. `lonfolder_ref = ""`
5. Run the simulations (see step 2e of the first part, with `samplefolder = HM + "/samples_initial_space"` and `simresult = HM + "/simulations_initial_space"`). Before that, change `basefolder + "/data/pfile_BASE.par"` (copied there by `4_simulation_assembling.py`) as
    - a. `dt = 25`
    - b. `tend = 150`
6. Download the results.
7. Run `5_extract_outputs_synthetic.py` with
    - a. `ECGsegment = "QRScomplex"` (`tend_truncated = 150`)
    - b. `paramfilename = "simulation_parameters_fitting_biv"`
    - c. `labelfilename_scalar = "xlabels_calibration_QRS"`
    - d. `labelfilename_gradient = ""`
    - e. `samplefolder = HM + "/samples_initial_space"`
    - f. `simulationfolder = HM + "/simulations_initial_space"`
    - g. `PCAbasisfolder = basefolder + "/PCAbasis_KCL1_150ms_99%"`
    - h. `qoifolder = HM + "/qoi_initial_space"`
8. Train the GPEs with `libGPE/gpe_training.py`, with
    - a. `ECGsegment = "QRScomplex"`
    - b. `waveno = 0` for the first training (initial space only). For the later waves, `waveno` = the last wave included in the training data. The script uses the data of `qoi_initial_space` and `qoi_wave_1` to `qoi_wave_<waveno>`, and writes the GPEs in `OUTPUT_DIR/KCL1/History_Matching/gpe/QRScomplex/initial_space` or `wave_<waveno>` (in a subfolder `GPE_Matern_GL`).
9. Perform the history matching with `libGPE/history_matching.py` (`ECGsegment = "QRScomplex"`, `waveno = 1` for the first wave). Wave n uses the GPEs of wave n-1 (`initial_space` for wave 1). It writes `X_simul_<n>.txt`, `X_test_<n>.txt` and `wave_<n>.json` in `HM + "/output"`, and the figures in `HM + "/figures"`.
10. Run `2_sampling_simul_bhm.py` with `waveno = n`. It writes the new samples in `HM + "/samples_wave_<n>"` and copies the label files into that folder.
11. Repeat steps 3 to 10 (excluding the change of `pfile_BASE.par` in step 5) for each wave n, with
    - a. `samplefolder = HM + "/samples_wave_<n>"` (steps 3, 4, 5 and 7)
    - b. `simresult` and `simulationfolder` = `HM + "/simulations_wave_<n>"`
    - c. `qoifolder = HM + "/qoi_wave_<n>"`

   and then `gpe_training.py` with `waveno = n` and `history_matching.py` with `waveno = n + 1`.

Notes:

- Remember that `gpe_training.py` adds the new qoi folder when `waveno` is increased.
- In the last waves, use `ionicmodelname = "ToRORd_dynCl"` in step 4 **[TO BE CONFIRMED: also `labelfilename_im = "ToRORd_dynCl"` in steps 3 and 4]**.
- For the last wave of the history matching, use `n_simuls_factor = 100` in `history_matching.py` (the number of new samples is then `len(xlabels) * 100`).

When the history matching is finished, select the best-fit point with `6_select_best_fit.py` (see the T-wave calibration, step 2).

## Pipeline: T-wave calibration

Settings checked against the scripts; the pipeline has not been run end to end.

Folders used below: `HM = basefolder + "/History_Matching/Twave"`.

1. Create the folder `HM`.
2. Define the parameter file `data/simulation_parameters_fitting_biv_best-fit_QRS.json`: a copy of `data/simulation_parameters_fitting_biv.json` in which the values of the parameters listed in `xlabels_calibration_QRS.txt` are replaced by their best-fit values from the QRS-complex calibration (a point chosen among the non-implausible points of the last history matching wave).
   It is easier to use one of the parameter files from the HM waves and just change the last few lines.
   Select the best-fit point (smallest RMSE between the emulated and the target ECG) with `6_select_best_fit.py`, with
    - i. `K = 100`
    - ii. `ECGsegment = "QRScomplex"`, `final_waveno` = the last wave of the history matching

   It writes `best_NIMP_wave_<n>.json` in the output folder of the history matching.
3. Run `2_sampling_calibration.py` with
    - a. `basefolder = OUTPUT_DIR + "/KCL/History_Matching/Twave"`
    - b. `paramfilename = "simulation_parameters_fitting_biv_best-fit_QRS"`
    - c. `labelfilename = "xlabels_calibration_Twave"` (the file is provided in `libEP/data/xlabels`)
    - d. `N_samples = len(xlabels) * 20` (already set in the script)
    - e. `samplefolder = HM + "/samples_initial_space"`
4. Define a sample folder for the calibrated QRS-complex parameters to compute fibres and lead field only once. It must be next to the sample folders (`HM + "/ref_QRS"`), because the scripts look for it there.
    - a. Make the directory `HM + "/ref_QRS/0"`
    - b. Copy `simulation_parameters_fitting_biv_best-fit_QRS.json` to `HM + "/ref_QRS/0/param.json"`
    - c. Copy `xlabels_calibration_QRS.txt` (needed to compute fibres, CVs and HPS) and `ToRORd_dynCl.txt` to `HM + "/ref_QRS"`
    - d. Run `3_sample_assembling_pool.py` with
        1. `labelfilename_scalar = "xlabels_calibration_QRS"`
        2. `samplefolder = HM + "/ref_QRS"`
    - e. Run `4_simulation_assembling.py` with
        1. `labelfilename_scalar = "xlabels_calibration_QRS"`
        2. `samplefolder = HM + "/ref_QRS"`
    - f. Compute the lead field for this folder: run `4bis_run_simulations.py` with `N_samples = 1`, `samplefolder = HM + "/ref_QRS"` and `simresult = HM + "/simulations_ref_QRS"` (the EP simulation is also run and can be ignored)
    - g. Rename `HM + "/simulations_ref_QRS/SIM_LF_0"` as `HM + "/ref_QRS/SIM_LF_ref"`
5. **[TO BE UPDATED]** Run `2bis_check_imp.py` (in combination with `run_bench_pool.py`).
6. Run `3_sample_assembling_pool.py` with
    - a. `labelfilename_im = "ToRORd_dynCl"`
    - b. `labelfilename_scalar = "xlabels_calibration_Twave"`
    - c. `samplefolder = HM + "/samples_initial_space"`

   Compute ONLY the ADJ files (HPS, CVs and fibres can be computed only once and then referenced to `ref_QRS`).
7. Run `4_simulation_assembling.py` with
    - a. `labelfilename_im = "ToRORd_dynCl"`
    - b. `labelfilename_scalar = "xlabels_calibration_Twave"`
    - c. `samplefolder = HM + "/samples_initial_space"`
    - d. `propagation = "RE+"`
    - e. `hps_ref = True` (if that is the case)
    - f. `cvfolder_ref = HM + "/ref_QRS/0/CV"` (if that is the case)
    - g. `lonfolder_ref = HM + "/ref_QRS/0"` (if that is the case)
8. Run the simulations (see step 2e of the first part, with `simresult = HM + "/simulations_initial_space"`) and download the results. Check first that `basefolder + "/data/pfile_BASE.par"` has `dt = 5` and `tend = 500` (the QRS-complex calibration changes them).
9. Run `5_extract_outputs_synthetic.py` with
    - a. `ECGsegment = "Twave"` (`tend = 500`)
    - b. `paramfilename = "simulation_parameters_fitting_biv_best-fit_QRS"`
    - c. `labelfilename_scalar = "xlabels_calibration_Twave"`
    - d. `labelfilename_gradient = ""`
    - e. `samplefolder = HM + "/samples_initial_space"`
    - f. `simulationfolder = HM + "/simulations_initial_space"`
    - g. `PCAbasisfolder = basefolder + "/PCAbasis_KCL1_REP_99%"`
    - h. `qoifolder = HM + "/qoi_initial_space"`
10. Train the GPEs with `libGPE/gpe_training.py`, with `ECGsegment = "Twave"` and `waveno` as in step 8 of the QRS-complex calibration.
11. Perform the history matching with `libGPE/history_matching.py` (`ECGsegment = "Twave"`).
12. Run `2_sampling_simul_bhm.py` with
    - a. `basefolder = OUTPUT_DIR + "/KCL/History_Matching/Twave"`
    - b. `paramfilename = "simulation_parameters_fitting_biv_best-fit_QRS"`
    - c. `labelfilename = "xlabels_calibration_Twave"`
13. Repeat from step 6 for each wave (with `samplefolder = HM + "/samples_wave_<n>"`, and the simulation and qoi folders of the wave, as in step 11 of the QRS-complex calibration).
