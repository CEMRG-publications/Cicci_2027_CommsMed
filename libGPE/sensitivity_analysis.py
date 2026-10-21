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

# make config.py (in the parent folder) importable
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# dependencies
from config import OUTPUT_DIR, check_paths
from utils import global_sensitivity_analysis
from utils import msfun
from GPErks.serialization.labels import read_labels_from_file


def main():
    check_paths("OUTPUT_DIR")

    # path to folder where the GPEs of the single-cell model are saved (see gpe_training.py)
    gpe_folder = OUTPUT_DIR + "/ToRORd_dynCl/GPE_Matern_GL"

    # path to folder where GSA results will be saved
    gsa_folder = OUTPUT_DIR + "/gsa"
    os.makedirs(gsa_folder, exist_ok=True)

    # choose cut-off threshold for GSA
    cutoff = 0.75 #0.9

    # mode to rank the parameters according to STi, for i in idx_features 
    mode = "max" # max | sum | mean
    gsa_mode = "STi_median" # "STi_median" | "Si_median"

    # ================================================================
    # IO folder(s)
    # ================================================================
    output_dir = gsa_folder + "/output"
    os.makedirs(output_dir, exist_ok=True)

    plot_dir = gsa_folder + "/figures"
    os.makedirs(plot_dir, exist_ok=True)

    data_dir = gpe_folder + "/data"

    # ================================================================
    # GSA Sobol' indices on each output
    # ================================================================
    # extract feature list from .txt file
    idx_feature = np.loadtxt(data_dir + "/features_idx_list.txt", dtype=int)
    if np.size(idx_feature) == 1:
        idx_feature = [int(idx_feature)]
    else:
        idx_feature = list(idx_feature)

    for idx in idx_feature:
        # current QoI folder
        idxgpefolder = gpe_folder + "/output/" + str(idx)
        idxfolder = output_dir + "/" + str(idx)

        # load labels
        xlabels = read_labels_from_file(data_dir + "/xlabels.txt")
        ylabels = read_labels_from_file(data_dir + "/ylabels.txt")[idx]

        print("\n\nCurrent feature: {}".format(ylabels))

        if not os.path.exists(idxfolder + "/sobol_indices_summary.txt"):
            global_sensitivity_analysis.compute_Sobol_indices(
                idxgpefolder, idxfolder, xlabels, ylabels, n_samples=1024, n_draws=1000
            )
        else:
            print(
                "GSA sobol_indices_summary.txt already found. Skipping computing Sobol indices."
            )

        # plot Sobol indices
        global_sensitivity_analysis.plot_Sobol_indices(idxfolder, xlabels, ylabels)

    # ================================================================
    # Param ranking & plot
    # ================================================================
    print("\nPlot GSA heatmap for each feature.\n")

    # load labels
    xlabels = read_labels_from_file(data_dir + "/xlabels.txt")
    ylabels = read_labels_from_file(data_dir + "/ylabels.txt")
    ylabels = [ylabels[i] for i in idx_feature]

    # load Sobol indices
    ST_all = np.zeros((len(xlabels), len(ylabels)), dtype=float)
    S1_all = np.zeros((len(xlabels), len(ylabels)), dtype=float)
    for i, idx in enumerate(idx_feature):
        S1 = np.loadtxt(output_dir + "/" + str(idx) + "/Si_median.txt")
        S1_all[:, i] = S1
        ST = np.loadtxt(output_dir + "/" + str(idx) + "/STi_median.txt")
        ST_all[:, i] = ST

    '''
    # plot heat map
    msfun.gsa_heat(
        ST_all,
        S1_all,
        xlabels,
        ylabels,
        plot_dir + "/",
        horizontal=True,
        normalise=False,
        xlabels_latex=xlabels,
        ylabels_latex=ylabels,
        fontsize=6,
        height=14,
        width=max(10, 10 * len(ylabels) / 20),
    )
    '''

    # plot input rank, with specified cut-off
    print(
        "\nPlot GSA parameters ranking for {} with {} cut-off.\n".format(
            gsa_mode, cutoff
        )
    )

    # plot parameter ranking
    msfun.gsa_parameters_ranking_S(
        data_dir + "/",
        output_dir + "/",
        gsa_mode=gsa_mode,
        mode=mode,
        uncertainty=False,
    )

    # plot parameter ranking: comulative explained variance above a given cut-off
    msfun.plot_rank_GSA(
        data_dir + "/",
        output_dir + "/",
        criterion=gsa_mode,
        mode=mode,
        figname=plot_dir + "/Rank_" + mode + "_" + gsa_mode + "_" + str(cutoff) + ".png",
        th=0.0,
        normalise=False,
        separate_colors=True,
        cutoff_th=cutoff,
    )

    '''
    # plot parameter ranking: Sobol indices >= cut-off * Sobol index of 1-st input
    msfun.plot_rank_GSA_mod(
        data_dir + "/",
        output_dir + "/",
        criterion="STi_median",
        mode=mode,
        figname=plot_dir
        + "/Rank_" + mode + "_STi_wrt_ST1x"
        + str(round(1 - cutoff, 2))
        + ".png",
        th=0.0,
        normalise=False,
        separate_colors=True,
        cutoff_th=1 - cutoff,
    )
    '''


if __name__ == "__main__":
    main()
