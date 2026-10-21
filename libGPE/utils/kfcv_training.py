#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file is part of 'Cicci_2027_CommsMed'.

__author__ = Ludovica Cicci
__email__ = l.cicci@imperial.ac.uk
"""

# dependencies
from GPErks.utils.random import set_seed
from GPErks.constants import HEIGHT, WIDTH
from GPErks.constants import (
    DEFAULT_RANDOM_SEED,
    DEFAULT_TRAIN_SNAPSHOT_DIR,
    DEFAULT_TRAIN_SNAPSHOT_RESTART_TEMPLATE,
    DEFAULT_TRAIN_SNAPSHOT_EPOCH_TEMPLATE,
)
from GPErks.gp.data.dataset import Dataset
from GPErks.gp.experiment import GPExperiment
from GPErks.gp.mean import LinearMean
from GPErks.perks.cross_validation import KFoldCrossValidation
from GPErks.train.early_stop import (
    GLEarlyStoppingCriterion,
    NoEarlyStoppingCriterion,
    PQEarlyStoppingCriterion,
    UPEarlyStoppingCriterion,
)
from GPErks.perks.inference import Inference
from GPErks.serialization.path import posix_path
from GPErks.utils.metrics import IndependentStandardError
from GPErks.train.emulator import GPEmulator
from GPErks.train.snapshot import NeverSaveSnapshottingCriterion

# extra imports
import os, sys
import numpy as np
import torch
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from gpytorch.likelihoods import GaussianLikelihood
from gpytorch.kernels import (
    MaternKernel,
    ScaleKernel,
    RBFKernel,
)
from torchmetrics import R2Score
import pickle

FOLD = 5
RESTART = 5
MAXEPOCHS = 1000
METRICS_LIST = [
    R2Score(),
    IndependentStandardError(),
]  # expand this list with other metrics from torchmetrics, if intrested (use R2Score as last metric)


def plot_inference(inference, savepath):
    fig, axis = plt.subplots(1, 1, figsize=(2 * WIDTH, 2 * HEIGHT / 3))

    idx_sort = np.argsort(
        inference.y_pred_mean
    )  # let's sort predicted values for a better visualisation
    x = np.arange(len(idx_sort))

    ci = 1.96  # 95% confidence interval

    axis.scatter(
        x,
        inference.y_test[idx_sort],
        facecolors="none",
        edgecolors="C0",
        label="observed",
    )
    axis.scatter(
        x,
        inference.y_pred_mean[idx_sort],
        facecolors="C0",
        s=16,
        label="predicted",
    )
    axis.errorbar(
        x,
        inference.y_pred_mean[idx_sort],
        yerr=ci * inference.y_pred_std[idx_sort],
        c="C0",
        ls="none",
        lw=0.5,
        label=f"uncertainty ({ci} STD)",
    )

    axis.set_xticks([])
    axis.set_xticklabels([])
    axis.legend(loc="upper left")

    fig.tight_layout()
    plt.savefig(savepath + "/" + "plot_gpe_inference.png")

    return


def gpe_training(savepath, x_, y_, xlabels, ylabels, kernel_type, esc_type):
    os.chdir(savepath)

    # enforce reproducibility
    seed = DEFAULT_RANDOM_SEED
    set_seed(seed)

    # ================================================================
    # Splitting the dataset
    # ================================================================
    # split dataset in training, validation and test sets
    x_train, x_test, y_train, y_test = train_test_split(
        x_, y_, test_size=0.2, random_state=seed  # used to be test_size=0.1
    )
    x, x_val, y, y_val = train_test_split(
        x_train, y_train, test_size=0.2, random_state=seed
    )

    # ================================================================
    # GPE training with K-fold cross-validation
    # ================================================================
    # build dataset
    dataset = Dataset(
        x,
        y,
        X_val=x_val,
        y_val=y_val,
        x_labels=xlabels,
        y_label=ylabels,
    )

    # define experiment
    likelihood = GaussianLikelihood()
    mean_function = LinearMean(degree=1, input_size=dataset.input_size, bias=True)

    print("Kernel of choice: " + kernel_type)
    if kernel_type == "Matern":
        kernel = ScaleKernel(MaternKernel(ard_num_dims=dataset.input_size))
    elif kernel_type == "RBF":
        kernel = ScaleKernel(RBFKernel(ard_num_dims=dataset.input_size))
    else:
        raise Exception("Kernel {} NOT implemented.")

    experiment = GPExperiment(
        dataset,
        likelihood,
        mean_function,
        kernel,
        n_restarts=RESTART,
        metrics=METRICS_LIST,
        seed=seed,
        learn_noise=True,
    )

    # dump experiment (skeleton) in config file, useful to reporduce the experiment
    config_file = "experiment.ini"
    experiment.save_to_config_file(config_file)

    # k-fold cross-validation training: the validation set is used to stop training at the
    # right moment for EVERY split; the accuracy of the split-specific emulator is obtained
    # by testing against the respective split-specific left-out part of the training dataset
    # device = "cuda" if torch.cuda.is_available() else "cpu"
    device = "cpu"
    devices = [device]
    kfcv = KFoldCrossValidation(experiment, devices, n_splits=FOLD, max_workers=1)

    # optimiser
    optimizer = torch.optim.Adam(experiment.model.parameters(), lr=0.1)

    # early stopping criteria
    print("Early stopping criterium of choice: " + esc_type)
    if esc_type == "GL":
        esc = GLEarlyStoppingCriterion(max_epochs=MAXEPOCHS, alpha=0.1, patience=8)
    elif esc_type == "PQ":
        esc = PQEarlyStoppingCriterion(
            max_epochs=MAXEPOCHS, alpha=1.0, patience=8, strip_length=5
        )
    elif esc_type == "UP":
        esc = UPEarlyStoppingCriterion(
            max_epochs=MAXEPOCHS, strip_length=5, successive_strips=4
        )
    else:
        print("Do not recognise early stopping criteria. Max epoch 500 is used!")
        esc = NoEarlyStoppingCriterion(max_epochs=MAXEPOCHS)

    # output the best model state (a dict) and its training statistics (an object)
    best_model_dct, best_train_stats_dct = kfcv.train(
        optimizer,
        early_stopping_criterion=esc,
    )

    # check training stats at each split
    best_epochs = []
    for _, bts in best_train_stats_dct.items():
        best_epochs.append(bts.best_epoch)

    # best-split emulator (according to highest R2-score) is automatically available within kfcv instance
    inference = Inference(kfcv.emulator)

    # output statistics to text file
    with open("training_summary.txt", "w") as f:
        sys_out = sys.stdout
        sys.stdout = f
        kfcv.summary()
        print("Best epoch {}\n".format(best_epochs))
        inference.summary()
        print("Best-split emulator fitted hyperparameters:")
        kfcv.emulator.hyperparameters()
        sys.stdout = sys_out

    # ================================================================
    # GPE training using the entire dataset
    # ================================================================
    # build a new dataset now including a testing set
    del dataset
    dataset = Dataset(
        x_train,
        y_train,
        X_test=x_test,
        y_test=y_test,
        x_labels=xlabels,
        y_label=ylabels,
    )

    # save dataset in pickle file
    file = open(savepath + "/dataset.pickle", "wb")
    pickle.dump(dataset, file)
    file.close()

    # define experiment
    likelihood = GaussianLikelihood()
    mean_function = LinearMean(degree=1, input_size=dataset.input_size, bias=True)

    if kernel_type == "Matern":
        kernel = ScaleKernel(MaternKernel(ard_num_dims=dataset.input_size))
    elif kernel_type == "RBF":
        kernel = ScaleKernel(RBFKernel(ard_num_dims=dataset.input_size))
    else:
        raise Exception("Kernel {} NOT implemented.")

    experiment = GPExperiment(
        dataset,
        likelihood,
        mean_function,
        kernel,
        n_restarts=RESTART,
        metrics=METRICS_LIST,
    )

    device = "cpu"  # "cuda" if torch.cuda.is_available() else "cpu"
    emulator = GPEmulator(experiment, device)

    optimizer = torch.optim.Adam(experiment.model.parameters(), lr=0.1)

    # making use of knowledge coming from the performed CV,
    # we run the training for an exact number of epochs
    max_epochs = int(np.mean(best_epochs))
    esc = NoEarlyStoppingCriterion(max_epochs)
    snpc = NeverSaveSnapshottingCriterion(
        posix_path(
            savepath,
            DEFAULT_TRAIN_SNAPSHOT_RESTART_TEMPLATE,
        ),
        DEFAULT_TRAIN_SNAPSHOT_EPOCH_TEMPLATE,
    )
    _, best_train_stats = emulator.train(
        optimizer, early_stopping_criterion=esc, snapshotting_criterion=snpc
    )
    # best_train_stats.plot(with_early_stopping_criterion=True)

    emulator.experiment.metrics = METRICS_LIST
    inference = Inference(emulator)

    with open("training_summary.txt", "a") as f:
        sys_out = sys.stdout
        sys.stdout = f
        print("\n*** Final GPE ***")
        inference.summary()
        sys.stdout = sys_out

    plot_inference(inference, savepath)

    # delete kfcv snapshots
    os.system("rm -r " + DEFAULT_TRAIN_SNAPSHOT_DIR)

    return





## TODO
