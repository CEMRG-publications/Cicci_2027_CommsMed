# Workflow - Cicci_2027_CommsMed

Status: Parts A, B and C are updated and checked against the scripts.
Parts marked **[TO BE UPDATED]** wait for the requirements file.

## General rules

- The user sets two paths, once, in `config.py`: `CARP_INSTALL_DIR` (CARPentry folder) and `OUTPUT_DIR` (folder for all results). `libEP_DIR` (the `libEP` folder) is also defined in `config.py`, and is not changed.
- All other settings are variables at the top of `main()` in each script.
- Examples below use `basefolder = OUTPUT_DIR + "/<meshname>"`.
- Set `meshname` (the base mesh name, for example `KCL1`) at the top of each script. The scripts build `basefolder` from it.
- The pre-processing (step **1**) creates `<meshname>_subendolayer`, `<meshname>_FEC` and `<meshname>_gregion` from the base mesh `<meshname>`. The later scripts add the suffix by themselves: `3_sample_assembling_pool.py` uses `_subendolayer` and `_gregion`, and `4_simulation_assembling.py` and `4bis_run_simulations.py` use `_gregion`.

## Requirements

- CARPentry release 2024 (academic licence), including meshtool and mguvc.
- Python 3.12.7 and the packages listed in the README **[TO BE UPDATED: requirements file]**.
- Torso mesh in `carp_txt` format, with the tag numbers set in `libEP_DIR + "/data/tag_HEART.json"`.
- Electrode pts-file(s) should be created ad-hoc and saved as `libEP_DIR + "/data/electrodes/<meshname>_electrodes.pts"`, with the following format:
    - first line: number of points
        - for the 12-lead ECG: 10
    - second line onwards: one line per electrode with x, y, z (same units as the mesh), in the correct order;
        - for a 12-lead ECG: V1, V2, V3, V4, V5, V6, RA, LA, RL, LL

    **NOTE:** second-to-last electrode used as ground.
- **[OPTIONAL]** A file with local activation times (one value per node of the biventricular mesh) for the pre-pacing of the cells, given with the argument `prepacing_lats` of `setup.write_par_tissue`. Without it, no pre-pacing is applied.

## Part A - 0D: single-cell model

The single-cell model (e.g., ToRORd_dynCl) and its global sensitivity analysis. This part does not depend on Part B.

### 0. Single-cell model and its global sensitivity analysis

This step is independent of the steps below.

- **a.** Run `0_single_cell.py` with
    - `num_proc = 8` (number of samples simulated at the same time)

    It samples the ionic model parameters (in `<basefolder> + "/samples"`), runs the single-cell simulations and writes the outputs in `<basefolder> + "/qoi"`.
- **b.** Run `libGPE/gpe_training.py` with
    - `single_cell = True`

    It trains the GPEs in `<basefolder> + "/gpe/GPE_Matern_GL"`.
- **c.** Run `libGPE/sensitivity_analysis.py`. It reads those GPEs and writes the Sobol indices and the figures in `<basefolder> + "/gsa"`.

---

## Part B - 3D: biventricular and torso model

The meshes, the simulations of the biventricular and torso model, the ECGs and the calibration.

### 1. Pre-process meshes and generate UVCs

- **a.** Set correct mesh tags in `libEP_DIR + "/data/tag_HEART.json"` and `libEP_DIR + "/data/tag_BATH.json"`
- **b.** Run `1_preprocessing.py` with
    - `meshname = <meshname>` (current mesh)

    **NOTE:** the input mesh must be in `carp_txt` format.

### 2. Compute PCA basis/bases

- **a.** Run `2_sampling.py` with
    - `paramfilename = "simulation_parameters"`
    - `labelfilename_scalar = "xlabels_scalar"`
    - `labelfilename_gradient = "xlabels_gradient"`
    - `rangefilename_scalar = "ranges_scalar"`
    - `rangefilename_gradient = "ranges_gradient"`
    - `N_samples = 10 * ( len(xlabels_scalar) + 3*len(xlabels_gradient) )`
    - `samplefolder = basefolder + "/samples_LHS"`
- **b.** **[OPTIONAL]** A preliminary single-cell check of the samples, to discard samples that are highly likely to fail (APD, VmMax or VmMin outside a viable range). It needs CARPentry (`bench.pt`). Set `check_single_cell = True` (and `num_proc`) in `2_sampling.py` in step **a**. If the sample folder already exists, the script skips the sampling and runs only the check. The simulations already completed are not run again. The accepted samples are listed in `accepted_samples.txt` (the later steps skip the others), and the samples with a missing output in `missing_statistics.txt`. To run the check again, remove `accepted_samples.txt` first.
- **c.** Run `3_sample_assembling_pool.py` with
    - `labelfilename_scalar = "xlabels_scalar"`
    - `samplefolder = basefolder + "/samples_LHS"`
    - `heterogeneity = True`
    - `labelfilename_im = "ToRORd_dynCl"`
- **d.** Run `4_simulation_assembling.py` with
    - `ionicmodelname = "ToRORd_dynCl"`
    - `labelfilename_im = "ToRORd_dynCl"`
    - `labelfilename_scalar = "xlabels_scalar"`
    - `samplefolder = basefolder + "/samples_LHS"`
    - `propagation = "RE+"`
    - `clusterfolder = basefolder` for a local run, otherwise the absolute path of the folder on the cluster
- **e.** Move or copy `basics/mesh` and `basics/data` (saved by step **2d** in `basefolder + "/basics"`) and the samples folder of the current step (`<samplefolder>`) into `clusterfolder`. Run the simulations (lead field first, then EP) for every sample:
    - locally: run `4bis_run_simulations.py` with
        - `samplefolder = clusterfolder + "/samples_LHS"`
        - `simulationfolder = clusterfolder + "/simulations_LHS"`
    - on a cluster: `clusterfolder` (step **d**) is the path on the cluster. Write your own job script with the `carp.pt` commands used in `4bis_run_simulations.py`. The folder with the samples keeps its name on the cluster.

    **NOTE:** simulations whose `lead_field_egms.dat` is missing or empty are treated as "failed" and skipped by step **f**.

- **f.** Run `5_extract_outputs_synthetic.py` with
    - for the QRS complex:
        - `ECGsegment = "QRScomplex"`
    - for the T-wave:
        - `ECGsegment = "Twave"`
    - `paramfilename = "simulation_parameters"`
    - `labelfilename_scalar = "xlabels_scalar"`
    - `labelfilename_gradient = "xlabels_gradient"`
    - `samplefolder = basefolder + "/samples_LHS"`
    - `simulationfolder = basefolder + "/simulations_LHS"`
    - `compute_PCA = True`

    The PCA basis is saved in `qoifolder` with a default name. Later steps, called with `compute_PCA = False`, will find it using this name.

### 3. Perform Local Sensitivity Analysis (LSA)

- **a.** Run `linear_fitting_sampling.py` with
    - `paramfilename = "simulation_parameters_fitting_biv"`
    - `labelfilename = "xlabels_fitting_biv"`
    - `rangefilename = "ranges_fitting_biv"`
    - `samplefolder = basefolder + "/samples_fitting"`
- **b.** Run `3_sample_assembling_pool.py` with
    - `labelfilename_scalar = "xlabels_fitting_biv"`
    - `samplefolder = basefolder + "/samples_fitting"`
    - `labelfilename_im = "ToRORd_dynCl"`
    - `heterogeneity = True`
- **c.** **Check that ALL samples have been accepted!** If not, find out why and modify the parameters of the rejected ones accordingly. At the end, all samples MUST have been accepted.
- **d.** Run `4_simulation_assembling.py` with
    - `ionicmodelname = "ToRORd_dynCl"`
    - `labelfilename_im = "ToRORd_dynCl"`
    - `labelfilename_scalar = "xlabels_fitting_biv"`
    - `samplefolder = basefolder + "/samples_fitting"`
    - `propagation = "RE+"`
    - `clusterfolder = basefolder` for a local run, otherwise the absolute path of the folder on the cluster
- **e.** Move or copy `basics/mesh` and `basics/data` (saved by step **3d** in `basefolder + "/basics"`) and the samples folder of the current step (`<samplefolder>`) into `clusterfolder`. Run the simulations (lead field first, then EP) for every sample:
    - locally: run `4bis_run_simulations.py` with
        - `samplefolder = clusterfolder + "/samples_fitting"`
        - `simulationfolder = clusterfolder + "/simulations_fitting"`
    - on a cluster: `clusterfolder` (step **d**) is the path on the cluster. Write your own job script with the `carp.pt` commands used in `4bis_run_simulations.py`. The folder with the samples keeps its name on the cluster.

    **NOTE:** **Check that ALL simulations have run successfully!** Simulations whose `lead_field_egms.dat` is missing or empty are treated as "failed" and skipped by step **f**.

- **f.** Run `5_extract_outputs_synthetic.py` with
    - `paramfilename = "simulation_parameters_fitting_biv"`
    - `labelfilename_scalar = "xlabels_fitting_biv"`
    - `labelfilename_gradient = ""`
    - `samplefolder = basefolder + "/samples_fitting"`
    - `simulationfolder = basefolder + "/simulations_fitting"`
    - `compute_PCA = False`
    - for the QRS complex:
        - `ECGsegment = "QRScomplex"`
        - `qoifolder = basefolder + "/qoi_fitting_150ms"`
    - for the T-wave:
        - `ECGsegment = "Twave"`
        - `qoifolder = basefolder + "/qoi_fitting_REP"`
- **g.** Run `linear_fitting_screening.py` with
    - `cutoff = 0.95`
    - `samplefolder = basefolder + "/samples_fitting"`
    - for the QRS complex:
        - `ECGsegment = "QRScomplex"`
        - `qoifolder = basefolder + "/qoi_fitting_150ms"`
    - for the T-wave:
        - `ECGsegment = "Twave"`
        - `qoifolder = basefolder + "/qoi_fitting_REP"`
- **h.** Define the label files of the parameters to be calibrated (one label per line, in alphabetical order):
    - `xlabels_calibration_QRS.txt`: parameters selected from the screening results in `qoi_fitting_150ms` (see `final_ranking_cumvar` and the cut-off figure)
    - `xlabels_calibration_Twave.txt`: parameters selected in the same way from `qoi_fitting_REP`, without the parameters already in the QRS-complex list

    Save both files into a temporary directory `basefolder + "/xlabels_calibration"` (see step **4**).

### 4. Generate a target...

Here: `<targetfolder>` = `basefolder + "/target"`.

#### ... from a simulation

- **4.1a.** Run `2_sampling_calibration.py` with
    - `hmfolder = <targetfolder>`
    - `paramfilename = "simulation_parameters_fitting_biv"`
    - `labelfilename = "xlabels_fitting_biv"`
    - `rangefilename = "ranges_fitting_biv"`
    - `N_samples = 20` (or any other number)
    - `samplefolder = hmfolder + "/samples"`
    - `check_single_cell = True` (single-cell check)
- **4.1b.** Run `3_sample_assembling_pool.py` with
    - `labelfilename_scalar = "xlabels_fitting_biv"`
    - `samplefolder = <targetfolder> + "/samples"`
    - `heterogeneity = True`
    - `labelfilename_im = "ToRORd_dynCl"`
- **4.1c.** Run `4_simulation_assembling.py` with
    - `ionicmodelname = "ToRORd_dynCl"`
    - `labelfilename_im = "ToRORd_dynCl"`
    - `labelfilename_scalar = "xlabels_fitting_biv"`
    - `samplefolder = <targetfolder> + "/samples"`
    - `propagation = "RE+"`
    - `clusterfolder = <targetfolder>` for a local run, otherwise the absolute path of the folder on the cluster
- **4.1d.** Move or copy `basics/mesh` and `basics/data` (saved by step **4.1c** in `basefolder + "/basics"`) and the samples folder of the current step (`<samplefolder>`) into `clusterfolder`. Run the simulations (lead field first, then EP) for every sample:
    - locally: run `4bis_run_simulations.py` with
        - `samplefolder = clusterfolder + "/samples"`
        - `simulationfolder = clusterfolder + "/simulations"`
    - on a cluster: `clusterfolder` (step **4.1c**) is the path on the cluster. Write your own job script with the `carp.pt` commands used in `4bis_run_simulations.py`. The folder with the samples keeps its name on the cluster.

    **NOTE:** simulations whose `lead_field_egms.dat` is missing or empty are treated as "failed" and skipped by step **4.1g**.

- **4.1e.** Compare the ECGs in `SIM_EP_<i>/ECG_lead_field.dat` with your own plotting tool, choose the target and remove the others. Rename to "0", that is:
    - `<targetfolder> + "/samples/<i>"` -> `<targetfolder> + "/samples/0"`
    - `<targetfolder> + "/simulations/SIM_LF_<i>"` -> `<targetfolder> + "/simulations/SIM_LF_0"`
    - `<targetfolder> + "/simulations/SIM_EP_<i>"` -> `<targetfolder> + "/simulations/SIM_EP_0"`
- **4.1f.** Copy the following files from `basefolder + "/xlabels_calibration"` to `<targetfolder> + "/samples"`:
    - `xlabels_calibration_QRS.txt`
    - `xlabels_calibration_Twave.txt`
- **4.1g.** Run `5_extract_outputs_synthetic.py` with
    - for the QRS complex:
        - `ECGsegment = "QRScomplex"`
        - `labelfilename_scalar = "xlabels_calibration_QRS"`
    - for the T-wave:
        - `ECGsegment = "Twave"`
        - `labelfilename_scalar = "xlabels_calibration_Twave"`
    - `paramfilename = "simulation_parameters_fitting_biv"`
    - `labelfilename_gradient = ""`
    - `samplefolder = <targetfolder> + "/samples"`
    - `simulationfolder = <targetfolder> + "/simulations"`
    - `qoifolder = <targetfolder> + "/qoi_" + ECGsegment`
    - `compute_PCA = False`

#### ... from a clinical 12-lead ECG

- **4.2a.** Prepare the clinical ECG as the text file `ECG_from_BSP.dat` with 13 columns: time, LI, LII, LIII, aVR, aVL, aVF, V1, V2, V3, V4, V5, V6. Time 0 must be the onset of the QRS complex. The sampling interval must be 2 ms, and the number of samples must equal that of the simulated ECGs (otherwise the projection on the PCA basis fails).
- **4.2b.** Create target folder `<targetfolder>`, and copy into it the following files:
    - `ECG_from_BSP.dat`
    - `xlabels_calibration_QRS.txt` (from `basefolder + "/xlabels_calibration"`)
    - `xlabels_calibration_Twave.txt` (from `basefolder + "/xlabels_calibration"`)
- **4.2c.** Run `5_extract_outputs_clinical.py` with
    - for the QRS complex:
        - `ECGsegment = "QRScomplex"`
    - for the T-wave:
        - `ECGsegment = "Twave"`

#### (last step for both synthetic and clinical)

- **4.1h. / 4.2d.** Run `libGPE/gpe_training.py` with
    - for the QRS complex:
        - `ECGsegment = "QRScomplex"`
    - for the T-wave:
        - `ECGsegment = "Twave"`
    - `train_gpes = False`
    - `single_cell = False`

    to construct the data folder and `outputs.csv` directly in the target qoi folder.

**NOTE:** each target qoi folder needs the file `data/xlabels.txt` (`<targetfolder> + "/qoi_<segment>/data/xlabels.txt"`). Its labels must match, in the same order, those of the GPEs (`data/xlabels.txt` in the GPE folder, see the calibration pipelines):

- Synthetic target: `gpe_training.py` writes it from the columns of X. Nothing to copy.
- Clinical target: there is no X, so the file is not written. Copy `xlabels_calibration_QRS.txt` to `<targetfolder> + "/qoi_QRScomplex/data/xlabels.txt"` and `xlabels_calibration_Twave.txt` to `<targetfolder> + "/qoi_Twave/data/xlabels.txt"`.

`history_matching.py` and `6_select_best_fit.py` read these files.


---

## Part C - 3D: history matching

### 5. Pipeline: QRS-complex calibration

Here: `<history_matching>` = `basefolder + "/History_Matching/QRScomplex"`.

In all the steps below, set `clusterfolder = <history_matching>`. The simulation results of each stage go in their own folder, so that stages do not overwrite each other.

- **a.** Run `2_sampling_calibration.py` with
    - `hmfolder = <history_matching>`
    - `paramfilename = "simulation_parameters_fitting_biv"`
    - `labelfilename = "xlabels_calibration_QRS"`
    - `rangefilename = "ranges_fitting_biv"`
    - `N_samples = len(xlabels) * 20` (already set in the script)
    - `samplefolder = <history_matching> + "/samples_initial_space"`

    The script copies the label files, including `xlabels_empty.txt` and `ToRORd_dynCl.txt`, into the sample folder.

    **NOTE:** `xlabels_calibration_QRS` must be defined based on the sensitivity analysis results.
- **b.** Run `3_sample_assembling_pool.py` with
    - `labelfilename_scalar = "xlabels_calibration_QRS"`
    - `labelfilename_im = "xlabels_empty"`
    - `samplefolder = <history_matching> + "/samples_initial_space"`
    - `heterogeneity = True`
- **c.** Run `4_simulation_assembling.py` with
    - `ionicmodelname = "TT2"` (for the first waves; use `ToRORd_dynCl` in the last waves)
    - `labelfilename_im = "xlabels_empty"`
    - `labelfilename_scalar = "xlabels_calibration_QRS"`
    - `samplefolder = <history_matching> + "/samples_initial_space"`
    - `propagation = "RE-"`
    - `hps_ref = False`
    - `cvfolder_ref = ""`
    - `lonfolder_ref = ""`
    - `clusterfolder = <history_matching>`
- **d.** Move or copy `basics/mesh` and `basics/data` (saved by step **c** in `basefolder + "/basics"`) into `<history_matching>`. Then change `<history_matching> + "/data/pfile_BASE.par"` as
    - `dt = 25`
    - `tend = 150`
- **e.** Locally, run `4bis_run_simulations.py` with (on a cluster, see step **2e** of Part B)
    - `clusterfolder = <history_matching>`
    - `samplefolder = <history_matching> + "/samples_initial_space"`
    - `simulationfolder = <history_matching> + "/simulations_initial_space"`
- **f.** Run `5_extract_outputs_synthetic.py` with
    - `ECGsegment = "QRScomplex"`
    - `paramfilename = "simulation_parameters_fitting_biv"`
    - `labelfilename_scalar = "xlabels_calibration_QRS"`
    - `labelfilename_gradient = ""`
    - `compute_PCA = False` (the PCA basis computed in step **2f** of Part B is read from `basefolder + "/PCAbasis_<meshname>_150ms_99%"`)
    - `samplefolder = <history_matching> + "/samples_initial_space"`
    - `simulationfolder = <history_matching> + "/simulations_initial_space"`
    - `qoifolder = <history_matching> + "/qoi_initial_space"`
- **g.** Train the GPEs with `libGPE/gpe_training.py`, with
    - `ECGsegment = "QRScomplex"`
    - `waveno = 0` for the first training (initial space only). For the later waves, `waveno` = the last wave included in the training data. The script uses the data of `qoi_initial_space` and `qoi_wave_1` to `qoi_wave_<waveno>`, and writes the GPEs in `<history_matching> + "/gpe/initial_space"` or `<history_matching> + "/gpe/wave_<waveno>"` (in a subfolder `GPE_Matern_GL`).
- **h.** Perform the history matching with `libGPE/history_matching.py` with
    - `ECGsegment = "QRScomplex"`
    - `waveno = 1` for the first wave

    Wave n uses the GPEs of wave n-1 (`initial_space` for wave 1). It writes `X_simul_<n>.txt`, `X_test_<n>.txt` and `wave_<n>.json` in `<history_matching> + "/output"`, and the figures in `<history_matching> + "/figures"`.
- **i.** Run `2_sampling_simul_bhm.py` with
    - `hmfolder = <history_matching>`
    - `waveno = n`

    It writes the new samples in `<history_matching> + "/samples_wave_<n>"` and copies the label files into that folder.

- **j.** Repeat steps **b** to **i** for each wave n, with
    - `samplefolder = <history_matching> + "/samples_wave_<n>"` (steps **b**, **c**, **e** and **f**)
    - `simulationfolder = <history_matching> + "/simulations_wave_<n>"`
    - `qoifolder = <history_matching> + "/qoi_wave_<n>"`
    - no new copy of `basics` and no new change of `pfile_BASE.par` (step **d**): `mesh` and `data` are already in `<history_matching>`

    and then `gpe_training.py` with `waveno = n` and `history_matching.py` with `waveno = n + 1`.

**Notes:**

- Remember that `gpe_training.py` adds the new qoi folder when `waveno` is increased.
- For the last wave of the history matching, use `n_simuls_factor = 100` in `history_matching.py` (the number of new samples is then `len(xlabels) * 100`).

When the history matching is finished, select the best-fit point with `6_select_best_fit.py` (see step **6b** of the T-wave calibration).

### 6. Pipeline: T-wave calibration

Here: `<history_matching>` = `basefolder + "/History_Matching/Twave"`.

In all the steps below, set `clusterfolder = <history_matching>`. The simulation results of each stage go in their own folder, so that stages do not overwrite each other.

- **a.** Create the folder `<history_matching>`.
- **b.** Define the parameter file `data/simulation_parameters_fitting_biv_best-fit_QRS.json`: a copy of `data/simulation_parameters_fitting_biv.json` in which the values of the parameters listed in `xlabels_calibration_QRS.txt` are replaced by their best-fit values from the QRS-complex calibration (a point chosen among the non-implausible points of the last history matching wave). It is easier to use one of the parameter files from the QRS-complex waves and just change the last few lines.

    Select the best-fit point (smallest RMSE between the emulated and the target ECG) with `6_select_best_fit.py`, with
    - `ECGsegment = "QRScomplex"`
    - `final_waveno` = the last wave of the history matching
    - `K = 100`
    - `flag_clinical = True` for a clinical target (the default `False` is for a synthetic target)

    It writes `best_NIMP_wave_<n>.json` in the output folder of the QRS-complex history matching.
- **c.** Run `2_sampling_calibration.py` with
    - `hmfolder = <history_matching>`
    - `paramfilename = "simulation_parameters_fitting_biv_best-fit_QRS"`
    - `labelfilename = "xlabels_calibration_Twave"` (the file is provided in `libEP_DIR + "/data/xlabels"`)
    - `rangefilename = "ranges_fitting_biv"`
    - `N_samples = len(xlabels) * 20` (already set in the script)
    - `samplefolder = <history_matching> + "/samples_initial_space"`
    - `check_single_cell = True` (single-cell check, see step **d**)
    - `num_proc` (number of samples simulated at the same time)

    The script copies the label files, including `xlabels_empty.txt` and `ToRORd_dynCl.txt`, into the sample folder.
- **d.** The single-cell check of the samples is run in step **c** with `check_single_cell = True` (T-wave calibration only). Look at `samples_initial_space/missing_statistics.txt` (if any). To resume after a failure, remove `accepted_samples.txt` and run `2_sampling_calibration.py` again with the same settings: the sampling is skipped and the simulations already completed are not run again.

- **e.** **[OPTIONAL]** Define a sample folder for the calibrated QRS-complex parameters, to compute HPS, conduction velocities, fibres and lead field only once, and save it to `<history_matching> + "/ref_QRS"`.
    - Make the directory `<history_matching> + "/ref_QRS/0"`
    - Copy `simulation_parameters_fitting_biv_best-fit_QRS.json` to `<history_matching> + "/ref_QRS/0/param.json"`
    - Copy `xlabels_calibration_QRS.txt` (needed to compute fibres, CVs and HPS), `xlabels_empty.txt` and `ToRORd_dynCl.txt` to `<history_matching> + "/ref_QRS"`
    - Run `3_sample_assembling_pool.py` with
        - `labelfilename_scalar = "xlabels_calibration_QRS"`
        - `labelfilename_im = "xlabels_empty"`
        - `samplefolder = <history_matching> + "/ref_QRS"`
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
    - Compute the lead field of the reference folder: run `4bis_run_simulations.py` locally (on a cluster, see step **2e** of Part B) with
        - `clusterfolder = <history_matching>`
        - `samplefolder = <history_matching> + "/ref_QRS"`
        - `simulationfolder = <history_matching> + "/simulations_ref_QRS"`

        The folder has only the sample 0. The EP simulation is also run and can be ignored.
    - Rename `<history_matching> + "/simulations_ref_QRS/SIM_LF_0"` as `<history_matching> + "/ref_QRS/SIM_LF_ref"`.
- **f.** Run `3_sample_assembling_pool.py` with
    - `labelfilename_scalar = "xlabels_calibration_Twave"`
    - `labelfilename_im = "ToRORd_dynCl"`
    - `samplefolder = <history_matching> + "/samples_initial_space"`
    - `heterogeneity = True`

    Only the adjustment files (ADJ) are needed here. The HPS, the CVs and the fibres can be computed only once and stored in `ref_QRS`. To save time, comment out the blocks for the EAS, the CV adjustment and the fibres in `3_sample_assembling_pool.py`, including the last call `setup.check_sample_folder(samplefolder)`.

- **g.** Run `4_simulation_assembling.py` with
    - `ionicmodelname = "ToRORd_dynCl"`
    - `labelfilename_im = "ToRORd_dynCl"`
    - `labelfilename_scalar = "xlabels_calibration_Twave"`
    - `samplefolder = <history_matching> + "/samples_initial_space"`
    - `propagation = "RE+"`
    - `hps_ref = True` (if that is the case)
    - `cvfolder_ref = <history_matching> + "/ref_QRS/0/CV"` (if that is the case)
    - `lonfolder_ref = <history_matching> + "/ref_QRS/0"` (if that is the case)
    - `clusterfolder = <history_matching>`

    Move or copy `basics/mesh` and `basics/data` (saved by this step in `basefolder + "/basics"`) into `<history_matching>`, unless this was already done in step **e**. Then check that `<history_matching> + "/data/pfile_BASE.par"` has `dt = 5` and `tend = 500`.

- **h.** Locally, run `4bis_run_simulations.py` with (on a cluster, see step **2e** of Part B)
    - `clusterfolder = <history_matching>`
    - `samplefolder = <history_matching> + "/samples_initial_space"`
    - `simulationfolder = <history_matching> + "/simulations_initial_space"`
- **i.** Run `5_extract_outputs_synthetic.py` with
    - `ECGsegment = "Twave"`
    - `paramfilename = "simulation_parameters_fitting_biv_best-fit_QRS"`
    - `labelfilename_scalar = "xlabels_calibration_Twave"`
    - `labelfilename_gradient = ""`
    - `compute_PCA = False` (the PCA basis computed in step **2f** of Part B is read from `basefolder + "/PCAbasis_<meshname>_REP_99%"`)
    - `samplefolder = <history_matching> + "/samples_initial_space"`
    - `simulationfolder = <history_matching> + "/simulations_initial_space"`
    - `qoifolder = <history_matching> + "/qoi_initial_space"`
- **j.** Train the GPEs with `libGPE/gpe_training.py`, with
    - `ECGsegment = "Twave"`
    - `waveno = 0` for the first training (initial space only). For the later waves, `waveno` = the last wave included in the training data (see step **5g** of the QRS-complex calibration).
- **k.** Perform the history matching with `libGPE/history_matching.py` with
    - `ECGsegment = "Twave"`
    - `waveno = 1` for the first wave

    It writes its files in `<history_matching> + "/output"` and `<history_matching> + "/figures"`, as in step **5h** of the QRS-complex calibration.
- **l.** Run `2_sampling_simul_bhm.py` with
    - `hmfolder = <history_matching>`
    - `waveno = n`
    - `paramfilename = "simulation_parameters_fitting_biv_best-fit_QRS"`
    - `labelfilename = "xlabels_calibration_Twave"`

    It writes the new samples in `<history_matching> + "/samples_wave_<n>"` and copies the label files into that folder.

- **m.** Repeat steps **f** to **l** for each wave n, with
    - `samplefolder = <history_matching> + "/samples_wave_<n>"` (steps **f**, **g**, **h** and **i**)
    - `simulationfolder = <history_matching> + "/simulations_wave_<n>"`
    - `qoifolder = <history_matching> + "/qoi_wave_<n>"`
    - no new copy of `basics` (step **g**): `mesh` and `data` are already in `<history_matching>`

    and then `gpe_training.py` with `waveno = n` and `history_matching.py` with `waveno = n + 1`.

**Notes:**

- Remember that `gpe_training.py` adds the new qoi folder when `waveno` is increased.
- For the last wave of the history matching, use `n_simuls_factor = 100` in `history_matching.py`.

When the history matching is finished, select the best-fit point with `6_select_best_fit.py` (`ECGsegment = "Twave"`, `final_waveno` = the last wave).
