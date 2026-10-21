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
import multiprocessing
import numpy as np

# make config.py (in the parent folder) importable
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# dependencies
from config import OUTPUT_DIR, check_paths
from utils import data_processing
from utils import kfcv_training
from GPErks.serialization.labels import read_labels_from_file


def main():
    check_paths("OUTPUT_DIR")

    # If True, train the GPEs. Otherwise, build data folder for sample "0" and stop
    train_gpes = True
    
    # If True, train GPEs on single-cell outputs. Otherwise, use ECG outputs
    single_cell = False

    if single_cell:
        gpe_folder_basename = OUTPUT_DIR + "/ToRORd_dynCl"
    else:
        gpe_folder_basename = OUTPUT_DIR + "/KCL1/History_Matching"

        # current ECG segment
        ECGsegment = "QRScomplex"

        if not train_gpes:
            # target folder for current case
            target_folder = OUTPUT_DIR + "/KCL1/target"
        else:
            # last wave included in the training data (0: initial space only)
            waveno = 0

    # choose kernel type
    kernel_type = "Matern"

    # choose early stopping criteria. Default: GL
    esc_type = "GL"

    if single_cell:
        qoi_folder_list = [OUTPUT_DIR + "/ToRORd_dynCl/qoi"]
        gpe_folder = gpe_folder_basename + "/GPE_" + kernel_type + "_" + esc_type
    elif train_gpes:
        # path(s) to folder with QoI(s): initial space and waves 1 to waveno
        qoi_basefolder = gpe_folder_basename + "/" + ECGsegment
        qoi_folder_list = [qoi_basefolder + "/qoi_initial_space"] + [
            qoi_basefolder + "/qoi_wave_" + str(k) for k in range(1, waveno + 1)
        ]

        # path to folder where GPE(s) will be saved. NOTE: can't be on an external disk
        wave_name = "initial_space" if waveno == 0 else "wave_" + str(waveno)
        gpe_folder = (
            gpe_folder_basename
            + "/gpe/"
            + ECGsegment
            + "/"
            + wave_name
            + "/GPE_"
            + kernel_type
            + "_"
            + esc_type
        )
    else:
        qoi_folder_list = [target_folder + "/qoi_" + ECGsegment]
        gpe_folder = qoi_folder_list[0]

    # ================================================================
    # Pre-process data
    # ================================================================
    os.makedirs(gpe_folder, exist_ok=True)
    if not os.path.exists(gpe_folder + "/outputs.csv"):
        # merge different datasets (if any) in a single .csv file
        if len(qoi_folder_list) > 1:
            for i in range(len(qoi_folder_list)):
                data_processing.write_csv(
                    qoi_folder_list[i], gpe_folder, ID="_" + str(i)
                )

            print("\nCombining data in a single csv file...")
            data_processing.combine_csv_files(
                gpe_folder, "outputs", len(qoi_folder_list)
            )
            data_processing.combine_csv_files(
                gpe_folder, "parameters", len(qoi_folder_list)
            )
        else:
            data_processing.write_csv(qoi_folder_list[0], gpe_folder)

    data_dir = gpe_folder + "/data"
    if not os.path.exists(data_dir):
        os.mkdir(data_dir)

        # write data in .txt file
        data_processing.write_txt(gpe_folder, data_dir)

    if not train_gpes:
        print("Data folder built in " + gpe_folder + ". No GPE trained.")
        return

    # ================================================================
    # IO folder(s)
    # ================================================================
    output_dir = gpe_folder + "/output"
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # ================================================================
    # GPE(s) training
    # ================================================================

    # extract feature list from .txt file
    idx_feature = np.loadtxt(data_dir + "/features_idx_list.txt", dtype=int)
    if np.size(idx_feature) == 1:
        idx_feature = [int(idx_feature)]
    else:
        idx_feature = list(idx_feature)

    # loop on each (selected) output of interest
    for idx in idx_feature:
        # current QoI folder
        idxfolder = output_dir + "/" + str(idx)
        if not os.path.exists(idxfolder):
            os.mkdir(idxfolder)

        # load dataset
        x_ = np.loadtxt(data_dir + "/X.txt", dtype=float)
        y_ = np.loadtxt(data_dir + "/Y.txt", dtype=float)
        if np.shape(y_[-1]):
            y_ = y_[:, idx]

        # load labels
        xlabels = read_labels_from_file(data_dir + "/xlabels.txt")
        ylabels = read_labels_from_file(data_dir + "/ylabels.txt")[idx]

        print("\n\nCurrent feature: {}".format(ylabels))

        # training (K-fold cross-validation. See kfcv_training.py to specify FOLD number)
        if not os.path.exists(idxfolder + "/best_model.pth"):
            kfcv_training.gpe_training(
                idxfolder, x_, y_, xlabels, ylabels, kernel_type, esc_type
            )
        else:
            print("GPE snapshot/best_model.pth already found. Skipping training.")

    # save summary training for all outputs in one file
    summaryfile = gpe_folder + "/final_summary_training"
    ylabels = read_labels_from_file(data_dir + "/ylabels.txt")
    if not os.path.exists(summaryfile + ".csv"):
        data_processing.write_summary(summaryfile, ylabels, output_dir)
    else:
        data_processing.write_summary(summaryfile + "_updated", ylabels, output_dir)


if __name__ == "__main__":
    # start the cross-validation workers as new processes (with "fork" they hang after the first output)
    multiprocessing.set_start_method("spawn", force=True)
    main()
