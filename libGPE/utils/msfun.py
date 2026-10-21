#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file is part of 'Cicci_2027_CommsMed'.

__author__ = Ludovica Cicci
__email__ = l.cicci@imperial.ac.uk

The functions in this file were extracted from the library of Marina Strocchi:
https://github.com/CEMRG-publications/Strocchi_2023_PLoSCompBio

Please cite:
Strocchi M, Longobardi S, Augustin CM, Gsell MAF, Petras A, Rinaldi CA, Vigmond EJ,
Plank G, Oates CJ, Wilkinson RD, Niederer SA. Cell to whole organ global sensitivity
analysis on a four-chamber heart electromechanics model using Gaussian processes
emulators. PLoS Comput Biol. 2023;19(6):e1011257. doi:10.1371/journal.pcbi.1011257
"""

# dependencies
from GPErks.serialization.labels import read_labels_from_file

# extra imports
import os
import numpy as np
import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def gsa_parameters_ranking_S(
    loadpath,
    loadpath_sobol,
    gsa_mode="STi_median",
    mode="max",
    threshold_cutoff=0.9,
    output_file=None,
    features_file=None,
    important_params_idx_file=None,
    uncertainty=True,
):
    """
    Ranks parameters using sensitivity indies.

    Args:
        - loadpath: datafolder containing data used for GPE training
        - loadpath_sobol: folder containing GPEs and GSA results
        - gsa_mode: output GSA file to read in [STi, S1]
        - mode: how to rank the parameters [max,mean,sum]
        - threshold_cutoff: to decide which parameters are important
        - output_file: output file to save ranking
        - features_file: features file containing which features to consider
        - important_params_idx_file: where to save the indices of the important parameters
        - uncertainty: if True, the GSA had to be run accounting for GPE uncertainty,
                       e.g. by sampling the posterior distribution of the GPEs N=1000 times
    """

    print("Checking folder structure...")
    to_check = [loadpath + "xlabels.txt", loadpath + "ylabels.txt", loadpath_sobol]
    for f in to_check:
        if not os.path.exists(f):
            raise Exception("Cannot find file " + f)
        else:
            print(f + " found.")

    xlabels = read_labels_from_file(loadpath + "xlabels.txt")
    ylabels = read_labels_from_file(loadpath + "ylabels.txt")

    if features_file is None:
        features_file = loadpath + "features_idx_list.txt"

    if not os.path.exists(features_file):
        raise Exception("Cannot find " + features_file)

    features = np.loadtxt(features_file, dtype=int)
    if len(features.shape) == 0:
        features = [features]
    else:
        features = list(features)

    print("Ranking parameters using features:")
    for idx in features:
        print(ylabels[idx])

    S = np.zeros((len(xlabels), len(features)), dtype=float)
    for i, idx in enumerate(features):
        if not os.path.exists(loadpath_sobol + str(idx) + "/" + gsa_mode + ".txt"):
            raise Exception(
                "Cannot find GSA file "
                + loadpath_sobol
                + str(idx)
                + "/"
                + gsa_mode
                + ".txt"
            )
        S_idx = np.loadtxt(loadpath_sobol + str(idx) + "/" + gsa_mode + ".txt")
        if uncertainty:
            S[:, i] = np.mean(S_idx, axis=0)
        else:
            S[:, i] = S_idx

    S_total = np.zeros((len(xlabels), 1), dtype=float)
    print("Ranking parameters according to their " + mode + " effect...")
    if mode == "mean":
        S_total = np.mean(S, axis=1)
    elif mode == "sum":
        S_total = np.sum(S, axis=1)
    elif mode == "max":
        S_total = np.max(S, axis=1)
    else:
        print("mode not recognised: please choose between mean, max and sum")

    ranked = np.argsort(S_total)
    ranked = ranked[::-1]
    ranked_S = S_total[ranked]

    if output_file is None:
        output_file = loadpath_sobol + "Rank_" + gsa_mode + "_" + mode + ".txt"

    f = open(output_file, "w")
    for i in range(len(xlabels)):
        f.write(xlabels[ranked[i]] + "\t" + str(ranked_S[i]) + "\n")
    f.close()

    print("Normalising ranked sensitivity to compute explained variance...")
    ranked_S_norm = list(np.array(ranked_S) / sum(ranked_S))

    ranked_S_norm_cumulative = []
    for i in range(len(xlabels)):
        ranked_S_norm_cumulative.append(sum(ranked_S_norm[0 : i + 1]))

    if output_file is None:
        output_file = (
            loadpath_sobol + "Rank_" + gsa_mode + "_" + mode + "_ExpVariance.txt"
        )
    else:
        output_file = output_file[:-4] + "_ExpVariance.txt"

    f = open(output_file, "w")
    for i in range(len(xlabels)):
        f.write(
            xlabels[ranked[i]]
            + "\t"
            + str(ranked_S_norm[i])
            + "\t"
            + str(ranked_S_norm_cumulative[i])
            + "\n"
        )
    f.close()

    if important_params_idx_file is not None:
        idx_cutoff = np.where(np.array(ranked_S_norm_cumulative) > threshold_cutoff)[0][
            0
        ]
        idx_param = ranked[range(idx_cutoff + 1)]
        np.savetxt(important_params_idx_file + "_cumulative", idx_param, fmt="%g")


def gsa_heat(
    ST,
    S1,
    xlabels,
    ylabels,
    savepath,
    height=7,
    width=5,
    horizontal=False,
    xlabels_latex=None,
    ylabels_latex=None,
    fontsize=6,
    normalise=False,
):
    """
    Plots a heatmap for a global sensitivity analysis.

    Args:
        - ST: total order effects matrix
        - S1: first order effects matrix
        - xlabels: parameter names
        - ylabels: output names
        - savepath: where to save the figure
        - height: figure height
        - width: figure width
        - correction: if you want to get rid of very small effects (below 0.01)
        - horizontal: if you want an horizontal heatmap or not
        - xlabels_latex: latex parameter names - for paper plotting
        - ylabels_latex: latex output names - for paper plotting
        - fontsize: size of figure font

    """
    if normalise:
        ST = ST / np.sum(ST, axis=0)

    """
    if (xlabels_latex is not None) or (ylabels_latex is not None):
        plt.rc("text", usetex=True)
    else:
        plt.rc("text", usetex=False)
    """ 
    plt.rc("text", usetex=False)
        
    if not horizontal:
        fig, axes = plt.subplots(1, 2, figsize=(width, height))
        rot_angle_x = 45
        rot_angle_y = 45
    else:
        fig, axes = plt.subplots(2, 1, figsize=(height, width))
        rot_angle_x = 90
        rot_angle_y = 0
    cbar_size = 0.25
    cbar_orientation = "vertical"
    cbar_args = {
        "shrink": cbar_size,
        "label": "Sensitivity",
        "orientation": cbar_orientation,
    }

    # First-order effects
    if not horizontal:
        df = pd.DataFrame(data=S1, index=xlabels, columns=ylabels)
    else:
        df = pd.DataFrame(data=np.transpose(S1), index=ylabels, columns=xlabels)

    h1 = sns.heatmap(
        df,
        cmap="rocket_r",
        vmin=0.0,
        vmax=1.0,
        square=True,
        linewidth=0.5,
        linecolor="black",
        cbar=True,
        cbar_kws=cbar_args,
        ax=axes[0],
    )

    # h1.set_yticks(np.arange(df.shape[0])+0.5)
    axes[0].set_title("First-order effects", fontsize=9, fontweight="bold")
    axes[0].tick_params(left=False, bottom=False)
    if xlabels_latex is None:
        h1.set_xticklabels(
            h1.get_xticklabels(), rotation=rot_angle_x, va="top", fontsize=fontsize
        )
    else:
        h1.set_xticklabels(
            xlabels_latex, rotation=rot_angle_x, va="top", fontsize=fontsize
        )

    if ylabels_latex is None:
        h1.set_yticklabels(
            h1.get_yticklabels(), rotation=rot_angle_y, ha="right", fontsize=fontsize
        )
    else:
        h1.set_yticklabels(
            ylabels_latex, rotation=rot_angle_y, va="top", fontsize=fontsize
        )

    # Total effects
    if not horizontal:
        df = pd.DataFrame(data=ST, index=xlabels, columns=ylabels)
        img_name = "heatmap.png"
    else:
        df = pd.DataFrame(data=np.transpose(ST), index=ylabels, columns=xlabels)
        img_name = "heatmap_horizontal.png"

    ht = sns.heatmap(
        df,
        cmap="rocket_r",
        vmin=0.0,
        vmax=1.0,
        square=True,
        linewidth=0.5,
        linecolor="black",
        cbar=True,
        cbar_kws=cbar_args,
        ax=axes[1],
    )

    # ht.set_yticks(np.arange(df.shape[0])+0.5)
    axes[1].set_title("Total effects", fontsize=9, fontweight="bold")
    axes[1].tick_params(left=False, bottom=False)
    if xlabels_latex is None:
        ht.set_xticklabels(
            ht.get_xticklabels(), rotation=rot_angle_x, va="top", fontsize=fontsize
        )
    else:
        ht.set_xticklabels(
            xlabels_latex, rotation=rot_angle_x, va="top", fontsize=fontsize
        )

    if ylabels_latex is None:
        ht.set_yticklabels(
            ht.get_yticklabels(), rotation=rot_angle_y, ha="right", fontsize=fontsize
        )
    else:
        ht.set_yticklabels(
            ylabels_latex, rotation=rot_angle_y, va="top", fontsize=fontsize
        )

    """
    if (xlabels_latex is not None) or (ylabels_latex is not None):
        plt.rc("text", usetex=True)
    else:
        plt.rc("text", usetex=False)
    """
    plt.rc("text", usetex=False)

    plt.tight_layout()
    plt.savefig(savepath + img_name, bbox_inches="tight", dpi=300)


def plot_rank_GSA(
    datapath,
    loadpath,
    rank_file=None,
    criterion="STi_median",
    mode="max",
    figname="",
    normalise=False,
    th=0.0,
    annotate=False,
    figsize=(18, 3),
    fontsize=6,
    xlabels_latex=None,
    separate_colors=False,
    color_important=None,
    cutoff_th=0.9,
):
    """
    Plots the parameter ranking for a .

    Args:
        - datapath: folder with data (xlabels.txt, etc...)
        - loadpath: path where you saved your parameter ranking
        - rank_file: if you want to provide a different parameter ranking file
                     that is not in the loadpath folder
        - criterion: STi or Si e.g. total or first order effects to use for ranking
        - mode: max or mean to rank the parameters
        - figname: name of output figure
        - normalise: if you want to normalise so that the values all sum up to 1
        - th: threshold to determine which parameters are important and which ones aren't
        - annotate: write numbers on top of each bar
        - figsize: size of output figure
        - fontsize: size of figure font
        - xlabels_latex: parameter names in latex for paper plots
        - separate_colors: if you want a different colour for important and unimportant parameters
        - color_important: what colour you want the important parameter bars to be
        - cutoff_th: fraction of output variance explained by the most important parameters
    """

    color = ["#ff8000"]

    # assumes that xlabels and ylabels are the same for all tests
    index_i = read_labels_from_file(datapath + "/xlabels.txt")

    x = np.arange(len(index_i))
    barWidth = 0.25

    fig, ax = plt.subplots(1, 1, figsize=figsize, constrained_layout=True)

    if rank_file is None:
        rank_file = loadpath + "/Rank_" + criterion + "_" + mode + ".txt"

    f = open(rank_file, "r")
    lines = f.readlines()

    r_dct = {}
    for line in lines:
        line_split = re.split(r"\t+", line)
        r_dct[line_split[0]] = float(line_split[1])

    bars = []
    for l in index_i:
        bars.append(r_dct[l])
    idx_sorted = np.argsort(np.array(bars))

    r = [xx + barWidth for xx in x]

    bars_sorted = [bars[idx] for idx in idx_sorted]
    bars_sorted = bars_sorted[::-1]
    if xlabels_latex is not None:
        index_i_sorted = [xlabels_latex[idx] for idx in idx_sorted]
    else:
        index_i_sorted = [index_i[idx] for idx in idx_sorted]
    index_i_sorted = index_i_sorted[::-1]

    bars_sorted_norm = list(np.array(bars_sorted) / sum(bars_sorted))

    bars_sorted_sum = []
    for i in range(len(bars_sorted)):
        bars_sorted_sum.append(sum(bars_sorted_norm[0 : i + 1]))

    if normalise:
        barplot = bars_sorted_norm
    else:
        barplot = bars_sorted

    plt.xticks(x + barWidth, index_i_sorted, rotation=90, fontsize=fontsize)
    ax.tick_params(axis="both", labelsize=fontsize)

    cutoff_param = np.where(np.array(bars_sorted_sum) > cutoff_th)[0][0]
    if color_important is None:
        color_important = "darkorange"
    color_unimportant = "lightgray"
    colors = [
        color_important,
    ] * (cutoff_param + 1) + [
        color_unimportant,
    ] * (len(bars_sorted_sum) - cutoff_param - 1)

    if separate_colors:
        bars = ax.bar(r, barplot, width=barWidth, edgecolor="white", color=colors)
    else:
        bars = ax.bar(r, barplot, color=color, width=barWidth, edgecolor="white")

    if "STi" in criterion:
        ax.set_ylabel("$ST$", fontsize=fontsize)
    else:
        ax.set_ylabel("$S_1$", fontsize=fontsize)

    if th > 0:
        ax.plot(
            [-2 * barWidth, len(index_i) + 2 * barWidth],
            [th, th],
            color="black",
            linestyle="--",
        )

    plt.legend()

    y_max = np.ceil(np.max(np.array(barplot)) * 10) / 10
    cutoff = r[cutoff_param] + 2 * barWidth
    ax.plot([cutoff, cutoff], [0, y_max], color="black", linestyle="--")

    if annotate:
        for i, bar in enumerate(bars):
            yval = bar.get_height()
            ax.text(
                bar.get_x(),
                yval + 0.005,
                str(round(bars_sorted_sum[i] * 100)) + "%",
                fontsize=fontsize,
            )

    if figname == "":
        plt.show()
    else:
        plt.savefig(figname, bbox_inches="tight", dpi=300)
