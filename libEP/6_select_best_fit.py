#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file is part of 'Cicci_2027_CommsMed'.

__author__ = Ludovica Cicci
__email__ = l.cicci@imperial.ac.uk
"""

# extra imports
import os
import sys
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt

# make config.py (in the parent folder) importable
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# dependencies
from config import OUTPUT_DIR, check_paths
from utils import IO
from utils import tools
from GPErks.gp.experiment import load_experiment_from_config_file
from GPErks.train.emulator import GPEmulator
from GPErks.serialization.path import posix_path


def main():
    check_paths("OUTPUT_DIR")

    meshname = "KCL1"  # NOTE: mesh basename
    basefolder = OUTPUT_DIR + "/" + meshname

    # current ECG segment: QRScomplex | Twave
    ECGsegment = "QRScomplex"

    # path to history matching basefolder
    hmfolder = basefolder + "/History_Matching/" + ECGsegment

    # final wave number
    final_waveno = 1

    # number of best-fitting samples to save
    K = 100

    # flag for clinical data
    flag_clinical = False

    # name of datafile with clinical ECG data
    if flag_clinical:
        bspfilename = "ECG_from_BSP"

    # path to folder with target simulation
    targetfolder = basefolder + "/target"

    # additional settings
    if "QRScomplex" in ECGsegment:
        # PCA basis folder
        PCAbasisfolder = basefolder + "/PCAbasis_" + meshname + "_150ms_99%"
        PCAbasisfilename = "WL_12leads_truncated"
        # simulation end time point
        tend = 150
    elif "Twave" in ECGsegment:
        # PCA basis folder
        PCAbasisfolder = basefolder + "/PCAbasis_" + meshname + "_REP_99%"
        PCAbasisfilename = "WL_12leads_REP"
        # simulation end time point
        tend = 500

    # path to HM output folder
    outputfolder = hmfolder + "/output"

    # path to folder with HM plots
    plotfolder = hmfolder + "/figures"

    # path to GPE folder
    if final_waveno == 1:
        gpefolder = hmfolder + "/gpe/initial_space"
    else:
        gpefolder = hmfolder + "/gpe/wave_" + str(final_waveno - 1)

    # ---> START
    # IO
    os.makedirs(plotfolder, exist_ok=True)
    valfolder = outputfolder + "/validation_wave_" + str(final_waveno) + "_K" + str(K)
    if not os.path.exists(valfolder):
        os.mkdir(valfolder)
    else:
        raise Exception(f"{valfolder} already exists!")

    # get target ECG
    if flag_clinical:
        ECG = np.genfromtxt(targetfolder + "/" + bspfilename + ".dat", delimiter=" ")
    else:
        ECG = IO.read_data(targetfolder + "/simulations/SIM_EP_0/ECG_lead_field")
    ECG_t = []
    for j in range(1, np.size(ECG, 1)):
        ECG_t.extend(ECG[: int(tend / 2) + 1, j])

    # load matrix of singular vectors WL of size p x L (L = #principal components, p = #ECG-inputs, i.e., timesteps)
    WL = np.genfromtxt(PCAbasisfolder + "/" + PCAbasisfilename + ".txt", dtype=float).T
    L = np.size(WL, axis=1)

    # get target inputs and corresponding labels
    target_inputs_labels = IO.read_text(
        targetfolder + "/qoi_" + ECGsegment + "/data/xlabels"
    )
    # get target PCA modes
    Y_target_df = pd.read_csv(targetfolder + "/qoi_" + ECGsegment + "/outputs.csv")

    # get output names
    qoi_list = list(Y_target_df.columns.values)
    if not L == len(qoi_list):
        raise Exception(f"Error! Inconsistent number of outputs. Check {PCAbasisfolder}/{PCAbasisfilename} and {targetfolder}/qoi_" + ECGsegment + "/outputs.csv")

    # get NIMP input parameters for current wave
    NIMP = IO.read_json(outputfolder + "/wave_" + str(final_waveno))["NIMP"]

    # K cannot be larger than the number of non-implausible points
    if K > len(NIMP):
        print(f"WARNING: K = {K} is larger than the number of NIMP points ({len(NIMP)}). Using K = {len(NIMP)}.")
        K = len(NIMP)

    # load last wave's NIMP points and predict PCA coefficients for each sample
    Y_pred = np.zeros((len(qoi_list), len(NIMP)))
    for q, qoi in enumerate(qoi_list):

        # find current output
        ylabels = IO.read_text(gpefolder + "/GPE_Matern_GL/data/ylabels")
        idx = np.where(np.array(ylabels) == qoi)[0][0]

        # current output GPE folder
        idxfolder = gpefolder + "/GPE_Matern_GL/output/" + str(idx)

        print(f"\nLoading GPE for output {qoi}...\n")

        # load dataset
        file = open(idxfolder + "/dataset.pickle", "rb")
        dataset = pickle.load(file)
        file.close()

        # load trained emulator
        experiment = load_experiment_from_config_file(
            idxfolder + "/experiment.ini", dataset
        )

        # device = "cuda" if torch.cuda.is_available() else "cpu"
        device = "cpu"
        emulator = GPEmulator(experiment, device)

        best_model_file = posix_path(idxfolder, "best_model.pth")
        emulator.load_state(best_model_file)
        emulator.hyperparameters()
        #  NOTICE: GPEs must have been trained using GPErks library (https://github.com/stelong/GPErks)

        # compute GPE-prediction of PCA coefficients for each sample
        Y_pred[q, :], _ = emulator.predict(np.array(NIMP))

    # recontruct GPE-based ECG
    ECG_gpe = np.matmul(Y_pred.T, WL.T)

    # compute root-mean-square error
    RMSE = tools.rmse_vectorised(ECG_t, ECG_gpe)

    # find index with smallest RMSE and save corresponding input parameters
    best_idx_RMSE = np.argmin(RMSE)
    best_NIMP_RMSE = NIMP[best_idx_RMSE]
    best_ECG_gpe_RMSE = ECG_gpe[best_idx_RMSE]

    # create dictionary of best inputs
    best_match_dct_RMSE = {}
    for j in range(len(best_NIMP_RMSE)):
        best_match_dct_RMSE[target_inputs_labels[j]] = best_NIMP_RMSE[j]
    IO.write_json(
        outputfolder + "/best_NIMP_wave_" + str(final_waveno),
        best_match_dct_RMSE,
    )

    # sort best signals according to RMSE
    sort_best_idx_RMSE = np.argsort(RMSE)

    # plot
    fig = plt.figure(figsize=(20, 5))
    plt.plot(
        best_ECG_gpe_RMSE.T,
        color="teal",
        label=f"GPE [{RMSE[best_idx_RMSE]:.2f}]",
    )

    # add target signal
    if "Twave" in ECGsegment:
        ECG_t_mod = []
        for it in range(12):
            ECG_t_mod = ECG_t_mod + [0] * 76 + ECG_t[76 + it * 251 : (it + 1) * 251 + 1]
        plt.plot(ECG_t_mod, "r--", label="Target")
    else:
        plt.plot(ECG_t, "r--", label="Target")
    plt.legend()
    plt.title(f"{meshname} - {ECGsegment} : RMSE best NIMP wave_{str(final_waveno)}")
    plt.savefig(plotfolder + "/" + meshname + "_best_NIMP_wave_" + str(final_waveno))

    # create dictionary of K best inputs for validation
    fig = plt.figure(figsize=(20, 5))
    for k in range(K):
        validation_dct_RMSE = {}
        for j in range(len(best_NIMP_RMSE)):
            validation_dct_RMSE[target_inputs_labels[j]] = NIMP[
                sort_best_idx_RMSE[k]
            ][j]
        IO.write_json(
            valfolder
            + "/"
            + str(k)
            + "_NIMP_wave_"
            + str(final_waveno),
            validation_dct_RMSE,
        )
        plt.plot(ECG_gpe[sort_best_idx_RMSE[k]].T, color="teal")

    # add target
    if tend == 500:
        ECG_t_mod = []
        for it in range(12):
            ECG_t_mod = (
                ECG_t_mod + [0] * 76 + ECG_t[76 + it * 251 : (it + 1) * 251 + 1]
            )
        plt.plot(ECG_t_mod, "r--", label="Target")
    else:
        plt.plot(ECG_t, "r--", label="Target")
    plt.legend()
    plt.title(f"{meshname} - {ECGsegment} : RMSE best {str(K)} NIMPs wave_{str(final_waveno)}: RMSE in [{RMSE[sort_best_idx_RMSE[0]]:.2f},{RMSE[sort_best_idx_RMSE[K-1]]:.2f}]")
    plt.savefig(f"{valfolder}/{meshname}_{ECGsegment}_NIMP_wave_{str(final_waveno)}")

    # print information
    print(f"Best NIMP according to (GPE-approx) ECG signal RMSE: {best_idx_RMSE}\n --> inputs = {best_NIMP_RMSE}")
        
    # <--- END


if __name__ == "__main__":
    main()
