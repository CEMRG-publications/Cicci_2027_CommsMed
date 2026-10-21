#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file is part of 'Cicci_2027_CommsMed'.

__author__ = Ludovica Cicci
__email__ = l.cicci@imperial.ac.uk
"""

# dependencies
from utils import IO

# extra imports
import os
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA as PCA


def compute_12lead_ecg(filename, data):
    """Compute 12-lead ECG from the electrical signal recorded at 10 electrodes locations (V1,V2,V3,V4,V5,V6,RA,RL,LA,LL).

    Parameters
    ----------
    filename: str
        Name of the ECG file.
    data: dict
        Body surface potentials (BSP) at the 10 electrodes.
    """

    # initialization
    ECG = {}

    # time
    ECG["t"] = range(len(data["V1"]))

    # defines the Wilson Central Terminal
    wct = data["LA"] + data["RA"] + data["LL"]

    # precordeal leads
    ECG["V1"] = data["V1"] - wct / 3
    ECG["V2"] = data["V2"] - wct / 3
    ECG["V3"] = data["V3"] - wct / 3
    ECG["V4"] = data["V4"] - wct / 3
    ECG["V5"] = data["V5"] - wct / 3
    ECG["V6"] = data["V6"] - wct / 3

    # limb Leads
    ECG["LI"] = data["LA"] - data["RA"]
    ECG["LII"] = data["LL"] - data["RA"]
    ECG["LIII"] = data["LL"] - data["LA"]

    # augmented leads
    ECG["aVR"] = data["RA"] - 0.5 * (data["LA"] + data["LL"])
    ECG["aVL"] = data["LA"] - 0.5 * (data["RA"] + data["LL"])
    ECG["aVF"] = data["LL"] - 0.5 * (data["LA"] + data["RA"])

    # save dat-file
    ecgdata = np.matrix(
        [
            ECG["t"],
            ECG["LI"],
            ECG["LII"],
            ECG["LIII"],
            ECG["aVR"],
            ECG["aVL"],
            ECG["aVF"],
            ECG["V1"],
            ECG["V2"],
            ECG["V3"],
            ECG["V4"],
            ECG["V5"],
            ECG["V6"],
        ]
    )
    IO.write_2D_data(filename, ecgdata.T)

    return ECG


def bsp_to_ecg(filename, data):
    """Convert 10-lead body surface potentials (BSP) into 12-lead ECG.

    Parameters
    ----------
    filename: str
        Name of the ECG file.
    data: numpy array
        BSP signals. Each column corresponds to one specific electrode, in the following order: V1,V2,V3,V4,V5,V6,RA,LA,RL,LL.
    """

    # extract the signal from each electrodes (given in the correct order: V1,V2,V3,V4,V5,V6,RA,LA,RL,LL)
    BSP = {
        "V1": data[:, 0],
        "V2": data[:, 1],
        "V3": data[:, 2],
        "V4": data[:, 3],
        "V5": data[:, 4],
        "V6": data[:, 5],
        "RA": data[:, 6],
        "LA": data[:, 7],
        "RL": data[:, 8],
        "LL": data[:, 9],
    }

    ECG = compute_12lead_ecg(filename, BSP)

    return ECG


def compute_PCA(
    X_list,
    xname,
    threshold_PCA=0.90,
    visualise=True,
    plotfolder="./figure",
):
    """Compute principal component analysis on data.

    Parameters
    ----------
    X_list : list
        Input data for each sample.
    xname: str
        Name of input quantity.
    threshold_PCA : float
        Ratio of variance to be explained by N principal components. Must be between in [0,1].
    visualise : bool
        If True, plot bar plot of explained varience against number of principal components. Default: True.
    plotfolder : str
        Path to folder where PCA plot will be saved. Default: "./figure".
    """

    if threshold_PCA < 0 or threshold_PCA > 1:
        raise Exception(
            "Ratio for the PCA must be in [0,1]. Current value {} NOT accepted.".format(
                threshold_PCA
            )
        )

    # IO
    if not os.path.exists(plotfolder):
        os.mkdir(plotfolder)

    # If 0 < n_components < 1 and svd_solver == "full", select the number of components such that the amount of variance that needs to be explained is greater than the percentage specified by n_components.
    if visualise:
        # full matrix of principal component scores
        pca = PCA(n_components=0.99, svd_solver="full")
        T = pca.fit_transform(np.array(X_list))

        var_ratio = 100 * pca.explained_variance_ratio_
        cum_var_ratio = 100 * np.cumsum(pca.explained_variance_ratio_)
        N_pca = np.where(cum_var_ratio >= 100 * threshold_PCA)[0][0] + 1
        Tl = T[:, :N_pca]

        # extract singular vectors/values
        svec_l = pca.components_[:N_pca, :]
        sval_l = pca.singular_values_[:N_pca]

        print(
            "Number of principal components capturing {}% of variation in {}: {}".format(
                100 * threshold_PCA, xname, N_pca
            )
        )

        # plot
        x = range(1, len(var_ratio) + 1)

        fig, ax1 = plt.subplots()
        ax1.bar(x, var_ratio, color="gray")
        ax1.set_ylabel("Explained variance [%]", color="gray")
        ax1.tick_params(axis="y", labelcolor="gray")

        ax1.set_xlabel("Principal component")
        ax1.set_ylim(0, 100 - 10 * int(int(var_ratio[0]) / 10))

        ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
        plt.plot(x, cum_var_ratio, color="red", marker="X")
        plt.axvline(
            x=N_pca, color="blue", linestyle="--", label="{}%".format(threshold_PCA)
        )
        ax2.set_ylabel("Cumulative Explained Variance [%]", color="red")
        ax2.tick_params(axis="y", labelcolor="red", grid_alpha=0)
        ax2.set_ylim(
            ax1.set_ylim()[0] + 10 * int(int(var_ratio[0]) / 10),
            ax1.set_ylim()[1] + 10 * int(int(var_ratio[0]) / 10),
        )

        plt.title(str(100 * threshold_PCA) + " % explined variance", color="blue")
        plt.savefig(plotfolder + "/PCA_" + xname + "_explained_variance.png")

        # plot singular value decay in a semilogy scale
        plt.figure()
        plt.plot(np.linspace(1, len(sval_l), len(sval_l)), sval_l, "b-o")
        plt.yscale("log")
        plt.xscale("log")
        plt.xlabel("index")
        plt.ylabel("Singular Values")
        plt.title("Decay of PCA singular values")
        plt.savefig(plotfolder + "/PCA_" + xname + "_singval_decay.png")

        # save singular value decay
        IO.write_2D_data(plotfolder + "/PCA_" + xname + "_singval_decay", sval_l)

    else:
        # truncated matrix of principal component scores
        pca = PCA(n_components=threshold_PCA, svd_solver="full")
        Tl = pca.fit_transform(np.array(X_list))
        N_pca = pca.n_components_

        # extract singular vectors/values
        svec_l = pca.components_
        sval_l = pca.singular_values_

        # save singular value decay
        IO.write_2D_data(plotfolder + "/PCA_" + xname + "_singval_decay", sval_l)

        print(
            "Number of principal components capturing {}% of variation in {}: {}".format(
                100 * threshold_PCA, xname, N_pca
            )
        )

    return Tl, svec_l








# vectorised Pearson correlation coefficient computation


# compute weighted Pearson correlation coefficients (based on ChatGPT)
    

# detect peaks in the ECG signal (Q, R, S, T waves)


# compute RMSE given a 1D-array as a target and a 2D-array of samples
def rmse_vectorised(target, samples):
    """
    Compute RMSE between a 1D target array and each row in a 2D samples array.

    Parameters
    ----------
    target: array/list
        1D array of shape (N,)
    samples: array
        2D array of shape (Nsamples, N)
    """
    # Ensure target is a row vector for broadcasting
    target = np.array(target).reshape(1, -1)
  
    # Compute squared differences, mean across features, then square root
    rmse = np.sqrt(np.mean((samples - target) ** 2, axis=1))
    
    return rmse

# compute relative RMSE given a 1D-array as a target and a 2D-array of samples (Bradley et al., 2000)


# compute normalised RMSE given a 1D-array as a target and a 2D-array of samples (Ushenin et al., 2021)