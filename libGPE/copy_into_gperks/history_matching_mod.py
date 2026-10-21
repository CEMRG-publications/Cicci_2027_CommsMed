"""
This file is part of 'Cicci_2027_CommsMed'.

Modified version of GPErks/perks/history_matching.py (GPErks commit dd8b7af).
GPErks: https://github.com/stelong/GPErks
GPErks is released under the MIT License, Copyright (c) 2021 Stefano Longobardi,
Gianvito Taneburgo. The licence text is in LICENSE_GPErks.txt.
Copy this file into GPErks/perks/ (GPErks commit dd8b7af).

Changes with respect to the GPErks file: plot_wave plots the inputs in several figures
(up to 19 inputs each) and can mark a target point; new method plot_wave_diagonal.

__author__ = Ludovica Cicci (modifications)
__email__ = l.cicci@imperial.ac.uk
"""

import matplotlib.gridspec as grsp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import iqr
import os

from GPErks.log.logger import get_logger
from GPErks.utils.array import get_minmax
from GPErks.utils.indices import diff, part_and_select, whereq_whernot
from GPErks.utils.jsonfiles import load_json, save_json

log = get_logger()


class Wave:
    """
    A module to perform Bayesian history matching using trained univariate emulators to
    match target distributions' mean and SD values.
    """

    def __init__(
        self,
        emulator=None,
        Itrain=None,
        cutoff=None,
        maxno=None,
        mean=None,
        var=None,
    ):
        self.emulator = emulator
        self.Itrain = Itrain
        self.cutoff = cutoff
        self.maxno = maxno
        self.mean = mean
        self.var = var
        self.I = None
        self.PV = None
        self.NIMP = None
        self.nimp_idx = None
        self.IMP = None
        self.imp_idx = None

    def compute_impl(self, X):
        n_samples = X.shape[0]
        output_dim = len(self.emulator)

        M = np.zeros((n_samples, output_dim), dtype=float)
        V = np.zeros((n_samples, output_dim), dtype=float)
        for j, emul in enumerate(self.emulator):
            mean, std = emul.predict(X)
            var = np.power(std, 2)
            M[:, j] = mean
            V[:, j] = var

        I = np.zeros((n_samples,), dtype=float)
        PV = np.zeros((n_samples,), dtype=float)
        for i in range(n_samples):
            In = np.sqrt((np.power(M[i, :] - self.mean, 2)) / (V[i, :] + self.var))
            PVn = V[i, :] / self.var

            I[i] = np.sort(In)[-self.maxno]
            PV[i] = np.sort(PVn)[-self.maxno]

        return I, PV

    def find_regions(self, X):
        n_samples = X.shape[0]
        I, PV = self.compute_impl(X)
        l = list(np.where(I < self.cutoff)[0])
        nl = diff(range(n_samples), l)

        self.I = I
        self.PV = PV
        self.nimp_idx = l
        self.NIMP = X[l]
        self.imp_idx = nl
        self.IMP = X[nl]

    def print_stats(self):
        nimp = len(self.nimp_idx)
        imp = len(self.imp_idx)
        tests = nimp + imp
        perc = 100 * nimp / tests

        stats = pd.DataFrame(
            index=["TESTS", "IMP", "NIMP", "PERC"],
            columns=["#POINTS"],
            data=[tests, imp, nimp, f"{perc:.4f} %"],
        )
        print(stats)

    def reconstruct_tests(self):
        n_samples = self.NIMP.shape[0] + self.IMP.shape[0]
        input_dim = self.NIMP.shape[1]
        X = np.zeros((n_samples, input_dim), dtype=float)
        X[self.nimp_idx] = self.NIMP
        X[self.imp_idx] = self.IMP
        return X

    def save(self, filename):
        dct = vars(self)
        excluded_keys = ["emulator"]
        obj_dct = {}
        obj_dct.update({k: dct[k] for k in set(list(dct.keys())) - set(excluded_keys)})
        save_json(obj_dct, filename)

    def load(self, filename):
        obj_dict = load_json(filename)
        for k, v in obj_dict.items():
            setattr(self, k, v)

    def get_nimps(self, n_points):
        nimp = len(self.nimp_idx)
        if n_points >= nimp - 1:
            raise ValueError(
                "Not enough NIMP points to choose from! n_points must be strictly "
                + "less than W.NIMP.shape[0] - 1."
            )
        else:
            _, X = part_and_select(self.NIMP, n_points)
            _, nl = whereq_whernot(self.NIMP, X)
        return X, self.NIMP[nl]

    # Note: the Wave object instance internal structure will be compromised after
    # calling this method: we recommend calling self.copy() and/or self.save()
    # beforehand!
    def augment_nimp(self, n_total_points, scaling=0.1, n_max=2000):
        X = np.copy(self.NIMP)
        lbounds = self.Itrain[:, 0]
        ubounds = self.Itrain[:, 1]

        log.info(
            f"\nRequested points: {n_total_points}\nAvailable points: "
            + f"{X.shape[0]}\nStart searching..."
        )

        count = 0
        a, b = (
            X.shape[0] if X.shape[0] < n_total_points else n_total_points,
            n_total_points - X.shape[0] if n_total_points - X.shape[0] > 0 else 0,
        )
        log.info(
            f"[Iteration: {count:<2}] Found: {a:<{len(str(n_total_points))}} "
            + f"({'{:.2f}'.format(100*a/n_total_points):>6}%) | "
            + f"Missing: {b:<{len(str(n_total_points))}}"
        )

        while X.shape[0] < n_total_points:
            count += 1

            bounds = get_minmax(X)
            scale = scaling * np.array(
                [bounds[i, 1] - bounds[i, 0] for i in range(X.shape[1])]
            )

            temp = np.random.normal(loc=X, scale=scale)
            count2 = 0
            while True:
                count2 += 1
                d1 = temp - lbounds.reshape((1, -1))
                d2 = ubounds.reshape((1, -1)) - temp
                flag = np.logical_or(
                    np.sum(np.sign(d1), axis=1) != temp.shape[1],
                    np.sum(np.sign(d2), axis=1) != temp.shape[1],
                )
                if count2 > n_max:
                    temp = temp[~flag, :]
                    break
                if np.sum(flag) > 0:
                    temp[flag, :] = np.random.normal(loc=X[flag, :], scale=scale)
                    continue
                else:
                    break

            I, _ = self.compute_impl(temp)
            nimp_idx = np.where(I < self.cutoff)[0]
            X = np.vstack((X, temp[nimp_idx]))

            a, b = (
                X.shape[0] if X.shape[0] < n_total_points else n_total_points,
                n_total_points - X.shape[0] if n_total_points - X.shape[0] > 0 else 0,
            )
            log.info(
                f"[Iteration: {count:<2}] Found: {a:<{len(str(n_total_points))}} "
                + f"({'{:.2f}'.format(100*a/n_total_points):>6}%) | "
                + f"Missing: {b:<{len(str(n_total_points))}}"
            )

        log.info("\nDone.")

        nimp = len(self.nimp_idx)
        _, NIMP_aug = part_and_select(X[nimp:], n_total_points - nimp)
        I, PV = self.compute_impl(NIMP_aug)
        self.NIMP = np.vstack((X[:nimp], NIMP_aug))
        self.I = np.concatenate((self.I, I))
        self.PV = np.concatenate((self.PV, PV))
        imp = len(self.imp_idx)
        self.nimp_idx += list(range(nimp + imp, imp + n_total_points))

    def get_trains(self, X, n_points):
        n_samples = X.shape[0]
        if n_points > n_samples:
            raise ValueError(
                "Cannot return more points than totally available points! "
                + "Choose n_points <= X_train.shape[0]."
            )
        elif n_points == n_samples:
            return X
        else:
            I, PV = self.compute_impl(X)
            l = np.argsort(I)[:n_points]
            return X[l]

    def copy(self):
        W = Wave(
            emulator=self.emulator,
            Itrain=self.Itrain,
            cutoff=self.cutoff,
            maxno=self.maxno,
            mean=self.mean,
            var=self.var,
        )
        W.I = np.copy(self.I)
        W.PV = np.copy(self.PV)
        W.NIMP = np.copy(self.NIMP)
        W.nimp_idx = self.nimp_idx.copy()
        W.IMP = np.copy(self.IMP)
        W.imp_idx = self.imp_idx.copy()
        return W

    def plot_wave(self, xlabels=None, display="impl", filepath=None, target=[]):
        X = self.reconstruct_tests()

        input_dim = X.shape[1]

        if xlabels is None:
            xlabels = [f"p{i+1}" for i in range(X.shape[1])]

        if display == "impl":
            C = self.I
            cmap = "jet"
            vmin = 1.0
            vmax = self.cutoff
            cbar_label = "Implausibility measure"

        elif display == "var":
            C = self.PV
            cmap = "bone_r"
            vmin = np.max(
                [
                    np.percentile(self.PV, 25) - 1.5 * iqr(self.PV),
                    self.PV.min(),
                ]
            )
            vmax = np.min(
                [
                    np.percentile(self.PV, 75) + 1.5 * iqr(self.PV),
                    self.PV.max(),
                ]
            )
            cbar_label = "GPE variance / EXP. variance"

        else:
            raise ValueError(
                "Not a valid display option! Can only display implausibilty maps "
                + "('impl') or proportion-of-exp.variance maps ('var')."
            )

        max_input_per_plot = np.min([input_dim,19])
        N_ss = int(np.ceil(input_dim / max_input_per_plot))
        input_dim_real = input_dim
        for subset in range(N_ss):

            if subset == N_ss-1:
                input_dim = input_dim_real - max_input_per_plot*(N_ss-1)
            else:
                input_dim = max_input_per_plot

            height = 9.36111
            width = 5.91667
            fig = plt.figure(figsize=(2 * width, 1.2 * 2 * height / 3))
            gs = grsp.GridSpec(
                input_dim - 1,
                input_dim,
                width_ratios=(input_dim - 1) * [1.0] + [0.1],
            )

            for k in range(input_dim * input_dim):
                i = k % input_dim
                j = k // input_dim

                if i > j:
                    axis = fig.add_subplot(gs[i - 1, j])
                    axis.set_facecolor("xkcd:light grey")

                    im = axis.hexbin(
                        X[:, j + subset*max_input_per_plot],
                        X[:, i + subset*max_input_per_plot],
                        C=C,
                        reduce_C_function=np.min,
                        gridsize=16,
                        cmap=cmap,
                        vmin=vmin,
                        vmax=vmax,
                    )

                    # add target point, if any
                    if target:
                        axis.plot(target[j + subset*max_input_per_plot],target[i + subset*max_input_per_plot],"m*",markersize=max_input_per_plot * 6 / input_dim,label="Target")

                    axis.set_xlim([self.Itrain[j + subset*max_input_per_plot, 0], self.Itrain[j + subset*max_input_per_plot, 1]])
                    axis.set_ylim([self.Itrain[i + subset*max_input_per_plot, 0], self.Itrain[i + subset*max_input_per_plot, 1]])

                    fontsize = int(12 * 5 / input_dim)
                    if i == input_dim - 1:
                        axis.set_xlabel(xlabels[j + subset*max_input_per_plot], fontsize=fontsize)
                    else:
                        axis.set_xticklabels([])
                    if j == 0:
                        axis.set_ylabel(xlabels[i + subset*max_input_per_plot], fontsize=fontsize)
                    else:
                        axis.set_yticklabels([])

                    plt.rc("xtick", labelsize=fontsize)
                    plt.rc("ytick", labelsize=fontsize)

            cbar_axis = fig.add_subplot(gs[:, input_dim - 1])
            cbar = fig.colorbar(im, cax=cbar_axis)
            cbar.set_label(cbar_label, fontsize=12 * 10 / max_input_per_plot)
            fig.tight_layout()
            if filepath is not None:
                filename = "subset" + str(subset+1) + "_" + os.path.basename(filepath)
                filefolder = os.path.dirname(filepath)
                plt.savefig(filefolder + "/" + filename, bbox_inches="tight", dpi=300)
                plt.close()
            else:
                plt.show()

 
    
    def plot_wave_diagonal(self, xlabels=None, filepath=None, target=[]):
        X = self.reconstruct_tests()
    
        input_dim = X.shape[1]
    
        if xlabels is None:
            xlabels = [f"p{i+1}" for i in range(X.shape[1])]
    
        C = self.I
        cmap = "jet"
        vmin = 1.0
        vmax = self.cutoff
        cbar_label = "Implausibility measure"
    
        max_input_per_plot = np.min([input_dim,20])
        N_ss = int(np.ceil(input_dim / max_input_per_plot))
        input_dim_real = input_dim
        for subset in range(N_ss):
    
            if subset == N_ss-1:
                input_dim = input_dim_real - max_input_per_plot*(N_ss-1)
            else:
                input_dim = max_input_per_plot
    
            height = 9.36111
            width = 5.91667
            fig = plt.figure(figsize=(2 * width, 1.2 * 2 * height / 3))
            gs = grsp.GridSpec(
                input_dim,
                input_dim+1,
                width_ratios=(input_dim) * [1.0] + [0.1],
            )
    
            for k in range(input_dim * input_dim):
                i = k % input_dim
                j = k // input_dim
    
                if i >= j: #i > j:
                    #axis = fig.add_subplot(gs[i - 1, j])
                    axis = fig.add_subplot(gs[i, j])
                    axis.set_facecolor("xkcd:light grey")
    
                    im = axis.hexbin(
                        X[:, j + subset*max_input_per_plot],
                        X[:, i + subset*max_input_per_plot],
                        C=C,
                        reduce_C_function=np.min,
                        gridsize=16,
                        cmap=cmap,
                        vmin=vmin,
                        vmax=vmax,
                    )
    
                    # add target point, if any
                    if target:
                        axis.plot(target[j + subset*max_input_per_plot],target[i + subset*max_input_per_plot],"m*",markersize=max_input_per_plot * 6 / input_dim,label="Target")
    
                    axis.set_xlim([self.Itrain[j + subset*max_input_per_plot, 0], self.Itrain[j + subset*max_input_per_plot, 1]])
                    axis.set_ylim([self.Itrain[i + subset*max_input_per_plot, 0], self.Itrain[i + subset*max_input_per_plot, 1]])
    
                    fontsize = int(12 * 5 / input_dim)
                    if i == input_dim - 1:
                        axis.set_xlabel(xlabels[j + subset*max_input_per_plot], fontsize=fontsize)
                    else:
                        axis.set_xticklabels([])
                    if j == 0:
                        axis.set_ylabel(xlabels[i + subset*max_input_per_plot], fontsize=fontsize)
                    else:
                        axis.set_yticklabels([])
    
                    plt.rc("xtick", labelsize=fontsize)
                    plt.rc("ytick", labelsize=fontsize)
    
            cbar_axis = fig.add_subplot(gs[:, input_dim])
            cbar = fig.colorbar(im, cax=cbar_axis)
            cbar.set_label(cbar_label, fontsize=12 * 10 / max_input_per_plot)
            fig.tight_layout()
            if filepath is not None:
                filename = "subset" + str(subset+1) + "_" + os.path.basename(filepath)
                filefolder = os.path.dirname(filepath)
                plt.savefig(filefolder + "/" + filename, bbox_inches="tight", dpi=300)
                plt.close()
            else:
                plt.show()
