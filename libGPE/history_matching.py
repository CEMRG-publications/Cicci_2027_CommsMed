#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file is part of 'Cicci_2027_CommsMed'.

__author__ = Ludovica Cicci
__email__ = l.cicci@imperial.ac.uk
"""

# other imports
import os
import sys

# make config.py (in the parent folder) importable
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# dependencies
from config import OUTPUT_DIR, check_paths
from GPErks.serialization.labels import read_labels_from_file
from utils import bayesian_history_matching


if __name__ == "__main__":
    check_paths("OUTPUT_DIR")

    # current ECG segment
    ECGsegment = "QRScomplex"

    # current wave ID
    waveno = 3

    basefolder = OUTPUT_DIR + "/KCL"

    # path to main folder with GPE(s) for current wave (trained on the data up to the previous wave)
    gpe_wave_name = "initial_space" if waveno == 1 else "wave_" + str(waveno - 1)
    gpe_folder = OUTPUT_DIR + "/gpe/KCL/" + ECGsegment + "/" + gpe_wave_name + "/GPE_Matern_GL"

    # path to Bayesian History Matching (BHM) folder
    bhm_folder = basefolder + "/History_Matching/" + ECGsegment

    # path to folder with target data (built with gpe_training.py, see workflow.txt)
    target_folder = basefolder + "/target/qoi_" + ECGsegment

    # threshold value for implausibility criterion
    cutoff = 3.0

    # number of test points per wave
    n_tests = 100000

    # number of points to simulate per wave: n_simuls_factor * number of parameters (use 100 in the last wave)
    n_simuls_factor = 10

    # ================================================================
    # IO folder(s)
    # ================================================================
    if not os.path.exists(bhm_folder):
        os.makedirs(bhm_folder)

    output_dir = bhm_folder + "/output"
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    plot_dir = bhm_folder + "/figures"
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)

    # ================================================================
    # Bayesian History Matching
    # ================================================================
    # path to txt-file with experimental mean and with experimental std, if any. Otherwise, use empty string "" to set the experimental std equal to 10% of experimental mean
    exp_mean_path = target_folder + "/data/Y"
    exp_std_path = ""

    # load labels for input parameters and output features
    xlabels = read_labels_from_file(gpe_folder + "/data/xlabels.txt")
    xlabels_target = read_labels_from_file(target_folder + "/data/xlabels.txt")
    ylabels = read_labels_from_file(gpe_folder + "/data/ylabels.txt")
    ylabels_target = read_labels_from_file(target_folder + "/data/ylabels.txt")

    # TODO: automatise the process
    if xlabels != xlabels_target:
        raise Exception("ATTENTION! There is no correspondence between:\n -->{}\n and\n -->{}.\nPlease, change {} and corresponding values in X.txt".format(gpe_folder + "/data/xlabels.txt",target_folder + "/data/xlabels.txt",gpe_folder + "/data/xlabels.txt"))

    bayesian_history_matching.history_matching(
        gpe_folder,
        output_dir,
        xlabels,
        ylabels,
        ylabels_target,
        exp_mean_path,
        exp_std_path=exp_std_path,
        waveno=waveno,
        cutoff=cutoff,
        maxno=1,
        n_tests=n_tests,
        n_simuls=len(xlabels) * n_simuls_factor,
        plot_dir=plot_dir,
        targetfile="", #target_folder + "/data/X",
    )
