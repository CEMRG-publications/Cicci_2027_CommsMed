#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file is part of 'Cicci_2027_CommsMed'.

__author__ = Ludovica Cicci
__email__ = l.cicci@imperial.ac.uk
"""

# dependencies

# extra imports
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA as PCA




def plot_ecg_linear_fitting(xlabels, ECG):
    """Plot and compare ECG for each input.

    Parameters
    ----------
    xlabels: list
        Input labels.
    ECG: array
        ECG signal. Must be defined accordingly to "linear_fitting_sampling.py" (see also "linear_fitting_screening.py").
    """

    Nparam = len(xlabels)
    for k in range(Nparam):
        plt.figure(figsize=(15, 9))
        plt.plot(ECG[2 * k, :], label="1")
        plt.plot(ECG[2 * k + 1, :], label="-1")
        plt.title(xlabels[k])
        plt.legend()

    return




def plot_cumvar(xlabels, var_array, plotfolder, L_max=None):
    """Plot and compare cumulative variances.

    Parameters
    ----------
    xlabels: list
        Input labels.
    var_array: numpy array
        Variance(s) of a polynomial of a random variable, for each input.
    L_max: int
        Maximum number of outputs to account for.
    """

    # account for all outputs
    if not L_max:
        L_max = np.size(var_array, axis=1)

    # rank according to absolute cumulative variance
    order_list = np.argsort(np.sum(var_array[:, :L_max], axis=1))
    var_array = var_array[order_list, :]
    xlabels = np.array(xlabels)[order_list].tolist()

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 25))
    inputs = [int(i) for i in range(len(xlabels))]
    start = np.zeros(np.shape(var_array[:, 0]))
    for l in range(L_max):
        if L_max < 20:
            col = plt.cm.tab20c(range(L_max))
        else:
            col = plt.cm.GnBu(np.linspace(0, 1, L_max))
        ax.barh(
            inputs,
            var_array[:, l],
            tick_label=xlabels,
            left=start,
            color=col[L_max - l - 1],
            alpha=0.8,
        )
        start += var_array[:, l]
        ax.set_title(
            "Cumulative sensitivity indices for {} principal components".format(L_max),
            fontsize=15,
        )
    plt.savefig(plotfolder + "/cumulative_variance.png", bbox_inches="tight")

    return


def plot_ranked_cutoff(xlabels, var, cutoff=0.9, figname=""):
    """Plot ranked normalised variances.

    Parameters
    ----------
    xlabels: list
        Input labels.
    var: array
        Variable variances.
    cutoff: float
        Cut-off value for most important parameters (i.e., those accounting for cutoff*100 % of output variance). Must be in [0,1]. Default: 0.9.
    figname: str
        Full path to figure to be saved. If "", do not save. Default: "".
    """

    if cutoff < 0 or cutoff > 1:
        raise Exception("Cut-off value must be between 0 and 1.")

    # rank inputs in descending according to the cumulative variance
    order_list = np.argsort(-var)
    xlabels_ranked = np.array(xlabels)[order_list].tolist()
    var_ranked = var[order_list]

    # normalise variance
    var_ranked_norm = var_ranked / np.sum(var_ranked)

    # keep only most important parameters
    index = np.where(np.cumsum(var_ranked_norm) >= cutoff)[0][0]

    # plot
    colors = ["lightgray"] * len(xlabels_ranked)
    for ind in range(index + 1):
        colors[ind] = "teal"  # "darkorange"
    fig, ax = plt.subplots(figsize=(len(xlabels_ranked) / 3, len(xlabels_ranked) / 9))
    plt.bar(xlabels_ranked, var_ranked_norm, color=colors)
    plt.xticks(rotation=90)

    if figname:
        figtitle = os.path.basename(figname)
        plt.title(
            "{}: {} selected parameters ({}% variance)".format(
                figtitle, int(index + 1), int(cutoff * 100)
            ),
            fontsize=16,
        )
        plt.savefig(figname + ".png", bbox_inches="tight")
    else:
        plt.title(
            "{}% cut-off on normalised cumulative variance".format(int(cutoff * 100)),
            fontsize=16,
        )

    return




def plot_final_ranking(
    final_criterium,
    xlabels,
    criteria_dct,
    cutoff=0.90,
    barcol="tab:blue",
    figname="",
):
    """Plot parameters and corresponding normalised variance (according to final criterium, e.g., cumvar or area-diff) and screen out non-influential parameters.

    Parameters
    ----------
    final_criterium: numpy array
        Parameter-variances based on final criterium (default order).
    xlabels: list
        List of labels (default order).
    criteria_dct: dict
        Dictionary with each ranking criterium and corresponding (default order) parameter-variances.
    cutoff: float.
        Cut-off value to decide which parameters are important. Must be in (0,1). Default: 0.90.
    barcol: str
        Bars colour in plot. Default: "tab:blue".
    figname: str
        Full path to figure to be saved. If "", do not save. Default: "".
    """

    # compute ranking based on final criterium
    final_order = np.argsort(-final_criterium)
    final_criterium_ranked = final_criterium[final_order]
    final_xlabels_ranked = np.array(xlabels)[final_order].tolist()

    # cut-off least important parameters based on normalised final criterium
    index = np.where(
        np.cumsum(final_criterium_ranked / np.sum(final_criterium)) > cutoff
    )[0][0]

    # plot
    facecol = [barcol[0]] * len(final_xlabels_ranked)
    for ind in range(index, len(final_xlabels_ranked)):
        facecol[ind] = "white"

    # order all criteria based on final ranking
    ordered_criteria_dct = {}
    for i, key in enumerate(criteria_dct):
        c = criteria_dct[key]
        ordered_criteria_dct[key] = c[final_order]

    # plot with cutoff
    df = pd.DataFrame(ordered_criteria_dct, index=final_xlabels_ranked)
    ax = df.plot.bar(
        width=0.75,
        color=barcol,
        rot=90,
        figsize=(len(final_xlabels_ranked) / 3, len(final_xlabels_ranked) / 9),
    )
    plt.bar(final_xlabels_ranked, final_criterium_ranked, color=facecol, alpha=0.3)
    ax.set_title(
        "{}% cut-off: non-influential parameters from {} onwards".format(
            int(cutoff * 100), final_xlabels_ranked[index]
        )
    )

    if figname:
        plt.savefig(figname + ".png", bbox_inches="tight")

    return


def plot_PCA_modes(X, paramname, paramvar, figname=""):
    """
    Parameters
    ----------
    X : numpy array
        Database of ECGs on which to perform PCA.
    paramname : str
        Name of current input parameter. If empty string is used, plots individual modes with 95% confidence interval.
    paramvar : float or list of floats
        Variance of current output(s) wrt current parameter.
    figname: str
        Full path to figure to be saved. If "", do not save. Default: "".
    """

    # perform a principal component analysis
    pca = PCA(n_components=0.99, svd_solver="full")
    obj = pca.fit(X)

    # check paramvar
    if type(paramvar) == float or type(paramvar) == int:
        paramvar = obj.n_components_ * [paramvar]
    if len(paramvar) > obj.n_components_:
        raise Exception("Wrong number of principal components. Check PCA basis.")

    # extract PCA quantities
    svec_l = obj.components_
    mean = obj.mean_
    sval_l = obj.singular_values_
    N_pca = len(paramvar)

    # threshold value for plots
    if paramname:
        plt.figure(figsize=(16, 5))
        t = np.arange(np.size(X, 1))
        colors = plt.cm.Spectral(np.linspace(0, 1, N_pca))
        for i in range(N_pca):
            principal_axis = svec_l[i]
            if paramvar[i] > 0.1:
                plt.fill_between(
                    t,
                    mean,
                    mean - np.sqrt(paramvar[i]) * principal_axis,
                    color=colors[i],
                    alpha=0.5,
                    label="{}".format(i + 1),
                )
            else:
                plt.fill_between(
                    t,
                    mean,
                    mean - np.sqrt(paramvar[i]) * principal_axis,
                    color=colors[i],
                    alpha=0.5,
                )
            plt.fill_between(
                t,
                mean,
                mean + np.sqrt(paramvar[i]) * principal_axis,
                color=colors[i],
                alpha=0.5,
            )
            plt.xlabel("Time [ms]")
            plt.title("Modes of variation wrt {}".format(paramname))
            plt.ylim([-2, 2])
            plt.legend(mode="expand", ncol=40, fontsize=8, handletextpad=0.15)
        plt.plot(t, mean, "k--", label="mean")
        plt.tight_layout()

        if figname:
            plt.savefig(figname + ".png", bbox_inches="tight")
    else:
        for i in range(N_pca):
            plt.figure(figsize=(16, 5))
            t = np.arange(np.size(X, 1))
            std = np.sqrt(sval_l[i])
            principal_axis = svec_l[i]
            plt.fill_between(
                t,
                mean,
                mean - 1.96 * std * principal_axis,
                alpha=0.5,
                color="tab:red",
                label="mean - 1.96*std",
            )
            plt.fill_between(
                t,
                mean,
                mean + paramvar[i] * std * principal_axis,
                color="tab:blue",
                alpha=0.5,
                label="mean + 1.96*std",
            )
            plt.xlabel("Time [ms]")
            plt.title("PC mode {}".format(i + 1))
            plt.ylim([-2, 2])
            plt.legend(mode="expand", ncol=2, fontsize=8, handletextpad=0.15)
            plt.plot(t, mean, "k--", label="mean")
            plt.tight_layout()

        if figname:
            plt.savefig(figname + "_mode_" + str(i + 1) + ".png", bbox_inches="tight")

    return








