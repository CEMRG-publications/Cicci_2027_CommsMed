#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file is part of 'Cicci_2027_CommsMed'.

__author__ = Ludovica Cicci
__email__ = l.cicci@imperial.ac.uk
"""

# dependencies
from GPErks.utils.random import set_seed
from GPErks.train.emulator import GPEmulator
from GPErks.gp.experiment import load_experiment_from_config_file
from GPErks.serialization.path import posix_path
from GPErks.perks.gsa_mod import SobolGSA
from GPErks.constants import DEFAULT_RANDOM_SEED
from GPErks.constants import HEIGHT, WIDTH

# extra imports
import os, sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pickle


def compute_Sobol_indices(
    loadpath, savepath, xlabels, ylabels, n_samples=1024, n_draws=1000, threshold=0.0001
):
    if not os.path.exists(savepath):
        os.mkdir(savepath)

    # load dataset
    file = open(loadpath + "/dataset.pickle", "rb")
    dataset = pickle.load(file)
    file.close()

    # load trained emulator
    experiment = load_experiment_from_config_file(loadpath + "/experiment.ini", dataset)
    # device = "cuda" if torch.cuda.is_available() else "cpu"
    device = "cpu"
    emulator = GPEmulator(experiment, device)

    best_model_file = posix_path(loadpath, "best_model.pth")
    emulator.load_state(best_model_file)
    emulator.hyperparameters()

    # enforce reproducibility
    seed = DEFAULT_RANDOM_SEED
    set_seed(seed)

    # GPE-based GSA using full posterior distribution samples to derive Sobol'
    # indices distributions
    gsa = SobolGSA(dataset, n=n_samples, seed=seed)
    gsa.estimate_Sobol_indices_with_emulator(emulator, n_draws=n_draws)
    gsa.correct_Sobol_indices(threshold=threshold)

    # plots
    gsa.plot_boxplot(savepath, ylabels)
    gsa.plot_heatmap(savepath, ylabels)

    # save Sobol' indices
    np.savetxt(
        savepath + "/" + "Si_median.txt",
        np.round(np.median(gsa.S1, axis=0), 6).reshape(-1, 1),
        fmt="%.6f",
    )
    np.savetxt(
        savepath + "/" + "STi_median.txt",
        np.round(np.median(gsa.ST, axis=0), 6).reshape(-1, 1),
        fmt="%.6f",
    )
    np.savetxt(savepath + "/" + "index_i.txt", gsa.index_i, fmt="%s")

    with open(savepath + "/" + "sobol_indices_summary.txt", "w") as f:
        sys_out = sys.stdout
        sys.stdout = f
        gsa.summary()
        sys.stdout = sys_out

    return


def plot_Sobol_indices(savepath, xlabels, ylabels):
    S1_mu = np.genfromtxt(savepath + "/" + "Si_median.txt").reshape(-1, 1)
    ST_mu = np.genfromtxt(savepath + "/" + "STi_median.txt").reshape(-1, 1)

    df_S1 = pd.DataFrame(data=S1_mu.T, columns=xlabels)
    df_ST = pd.DataFrame(data=ST_mu.T, columns=xlabels)

    fig, axes = plt.subplots(
        2, 1, figsize=(max(8 * WIDTH, 8 * WIDTH * len(ylabels) / 50), 2 * HEIGHT / 3)
    )

    h1 = sns.heatmap(
        df_S1,
        cmap="rocket_r",
        vmin=0.0,
        vmax=1.0,
        square=True,
        linewidth=0.1,
        cbar_kws={"shrink": 0.8},
        ax=axes[0],
    )
    axes[0].set_title("S1", fontsize=16, fontweight="bold")
    axes[0].tick_params(left=False, bottom=False)
    h1.set_xticklabels(h1.get_xticklabels(), rotation=90, va="top", size=8)
    h1.set_yticklabels(h1.get_yticklabels(), rotation=0, ha="right", size=10)

    ht = sns.heatmap(
        df_ST,
        cmap="rocket_r",
        vmin=0.0,
        vmax=1.0,
        square=True,
        linewidth=0.1,
        cbar_kws={"shrink": 0.8},
        ax=axes[1],
    )
    axes[1].set_title("ST", fontsize=16, fontweight="bold")
    axes[1].tick_params(left=False, bottom=False)
    ht.set_xticklabels(ht.get_xticklabels(), rotation=90, va="top", size=8)
    ht.set_yticklabels(ht.get_yticklabels(), rotation=0, ha="right", size=10)

    plt.suptitle(
        f"{ylabels} global sensitivity to model parameters",
        fontsize=21,
        fontweight="bold",
    )
    fig.tight_layout()
    plt.savefig(savepath + "/" + "plot_gsa_heatmap.png")
    plt.close()

    return
