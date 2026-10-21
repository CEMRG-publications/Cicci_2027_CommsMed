"""
This file is part of 'Cicci_2027_CommsMed'.

Modified version of GPErks/perks/gsa.py (SALib-based version of June 2024).
GPErks: https://github.com/stelong/GPErks
GPErks is released under the MIT License, Copyright (c) 2021 Stefano Longobardi,
Gianvito Taneburgo. The licence text is in LICENSE_GPErks.txt.
Copy this file into GPErks/perks/ (GPErks commit dd8b7af).

Changes with respect to the GPErks file: Sobol indices computed with SALib (first-order
and total-order only), own boxplot and heatmap functions that save the figures.

__author__ = Ludovica Cicci (modifications)
__email__ = l.cicci@imperial.ac.uk
"""

from typing import Callable, Optional

import numpy as np
import pandas as pd
from SALib.analyze import sobol
from SALib.sample import saltelli

from GPErks.constants import HEIGHT, WIDTH
import matplotlib.pyplot as plt
import matplotlib.gridspec as grsp
import seaborn as sns

from GPErks.constants import (
    DEFAULT_GSA_CONF_LEVEL,
    DEFAULT_GSA_N,
    DEFAULT_GSA_N_BOOTSTRAP,
    DEFAULT_GSA_N_DRAWS,
    DEFAULT_GSA_SKIP_VALUES,
    DEFAULT_GSA_THRESHOLD,
    DEFAULT_GSA_Z,
)
from GPErks.gp.data.dataset import Dataset
from GPErks.train.emulator import GPEmulator
from GPErks.utils.array import get_minmax


class SobolGSA:
    """
    Sobol global sensitivity analysis using Saltelli method and
    integrating emulator uncertainty.
    """

    def __init__(
        self,
        dataset: Dataset,
        n: int = DEFAULT_GSA_N,
        seed: Optional[int] = None,
    ):
        self.n = n
        self.seed = seed

        self.d = dataset.input_size
        self.index_i = dataset.x_labels
        self.ylabel = dataset.y_label
        self.minmax = (
            get_minmax(dataset.X_train)
            if (dataset.l_bounds is None and dataset.u_bounds is None)
            else np.hstack(
                (
                    np.array(dataset.l_bounds).reshape(-1, 1),
                    np.array(dataset.u_bounds).reshape(-1, 1),
                )
            )
        )

        self.ST = np.zeros((0, self.d), dtype=float)
        self.S1 = np.zeros((0, self.d), dtype=float)

        self.ST_std = np.zeros((0, self.d), dtype=float)
        self.S1_std = np.zeros((0, self.d), dtype=float)

        self.current_axis: Optional[np.ndarray] = None

    def assemble_Saltelli_space(self):
        problem = {
            "num_vars": self.d,
            "names": self.index_i,
            "bounds": self.minmax,
        }
        X = saltelli.sample(
            problem,
            self.n,
            calc_second_order=False,
            skip_values=DEFAULT_GSA_SKIP_VALUES,
        )
        return problem, X

    def estimate_Sobol_indices_with_emulator(
        self,
        emulator: GPEmulator,
        n_draws: int = DEFAULT_GSA_N_DRAWS,
    ):
        problem, X = self.assemble_Saltelli_space()
        Y = emulator.sample(X, n_draws)
        self._estimate_Sobol_indices_from_evaluations(problem, Y)

    def estimate_Sobol_indices_with_simulator(
        self,
        f: Callable[[np.ndarray], np.ndarray],
    ):
        problem, X = self.assemble_Saltelli_space()
        Y = np.squeeze(f(X)).reshape(1, -1)
        self._estimate_Sobol_indices_from_evaluations(problem, Y)

    def _estimate_Sobol_indices_from_evaluations(self, problem, Y):
        for y in Y:
            S = sobol.analyze(
                problem,
                y,
                calc_second_order=False,
                num_resamples=DEFAULT_GSA_N_BOOTSTRAP,
                conf_level=DEFAULT_GSA_CONF_LEVEL,
                parallel=False,
                n_processors=0,
                seed=self.seed,
            )
            T_Si, first_Si, (_, second_Si) = sobol.Si_to_pandas_dict(S)

            self.ST = np.vstack((self.ST, T_Si["ST"].reshape(1, -1)))
            self.S1 = np.vstack((self.S1, first_Si["S1"].reshape(1, -1)))

            self.ST_std = np.vstack(
                (self.ST_std, T_Si["ST_conf"].reshape(1, -1) / DEFAULT_GSA_Z)
            )
            self.S1_std = np.vstack(
                (
                    self.S1_std,
                    first_Si["S1_conf"].reshape(1, -1) / DEFAULT_GSA_Z,
                )
            )

    def correct_Sobol_indices(self, threshold: float = DEFAULT_GSA_THRESHOLD):
        for S in [self.ST, self.S1]:
            Q1 = np.percentile(S, q=25, axis=0)
            l = np.where(Q1 < threshold)[0]
            S[:, l] = np.zeros((S.shape[0], len(l)), dtype=float)

    def summary(self):
        df_Si = pd.DataFrame(
            data=np.round(np.median(self.S1, axis=0), 6).reshape(-1, 1),
            index=self.index_i,
            columns=["Si"],
        )        
        df_STi = pd.DataFrame(
            data=np.round(np.median(self.ST, axis=0), 6).reshape(-1, 1),
            index=self.index_i,
            columns=["STi"],
        )
        print(df_Si)
        print(df_STi)

    def plot_boxplot(self, savepath, ylabel):
        df_S1 = pd.DataFrame(data=self.S1, columns=self.index_i)
        df_ST = pd.DataFrame(data=self.ST, columns=self.index_i)

        #plt.style.use("seaborn")
        gs = grsp.GridSpec(2, 1)
        fig = plt.figure(figsize=(4 * WIDTH, 2 * HEIGHT))

        ax0 = fig.add_subplot(gs[0, 0])
        sns.boxplot(ax=ax0, data=df_S1)
        ax0.set_ylim(0, 1)
        ax0.set_title("First-order effect", fontweight="bold", fontsize=16, loc="left")
        ax0.set_xticklabels(ax0.get_xticklabels(), rotation=90, size=10)

        ax1 = fig.add_subplot(gs[1, 0])
        sns.boxplot(ax=ax1, data=df_ST)
        ax1.set_ylim(0, 1)
        ax1.set_title("Total effect", fontweight="bold", fontsize=16, loc="left")
        ax1.set_xticklabels(ax1.get_xticklabels(), rotation=90, size=10)

        plt.suptitle(
            f"{ylabel} global sensitivity to model parameters",
            fontsize=21,
            fontweight="bold",
        )
        fig.tight_layout()
        plt.savefig(savepath + "/" + "plot_gsa_boxplot.png")


    def plot_heatmap(self, savepath, ylabel):
        #S1_mean = np.round(np.mean(self.S1, axis=0), 6).reshape(-1, 1)
        #ST_mean = np.round(np.mean(self.ST, axis=0), 6).reshape(-1, 1)
        S1_mean = np.round(np.median(self.S1, axis=0), 6).reshape(-1, 1)
        ST_mean = np.round(np.median(self.ST, axis=0), 6).reshape(-1, 1)

        df_S1 = pd.DataFrame(data=S1_mean.T, columns=self.index_i)
        df_ST = pd.DataFrame(data=ST_mean.T, columns=self.index_i)

        fig, axes = plt.subplots(2, 1, figsize=(8 * WIDTH  * len(self.index_i)/65, 2 * HEIGHT / 3))
        fntsz = 16  * len(self.index_i)/13

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
        axes[0].set_title("S1", fontsize=fntsz, fontweight="bold")
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
        axes[1].set_title("ST", fontsize=fntsz, fontweight="bold")
        axes[1].tick_params(left=False, bottom=False)
        ht.set_xticklabels(ht.get_xticklabels(), rotation=90, va="top", size=8)
        ht.set_yticklabels(ht.get_yticklabels(), rotation=0, ha="right", size=10)

        plt.suptitle(
            f"{ylabel} sensitivity",
            fontsize=5 + fntsz,
            fontweight="bold",
        )
        fig.tight_layout()
        plt.savefig(savepath + "/" + "plot_gsa_heatmap.png")
        