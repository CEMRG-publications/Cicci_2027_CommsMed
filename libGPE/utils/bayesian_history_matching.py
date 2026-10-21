#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file is part of 'Cicci_2027_CommsMed'.

__author__ = Ludovica Cicci
__email__ = l.cicci@imperial.ac.uk
"""

# dependencies
from GPErks.utils.random import set_seed
from GPErks.constants import DEFAULT_RANDOM_SEED
from GPErks.gp.experiment import load_experiment_from_config_file
from GPErks.train.emulator import GPEmulator
from GPErks.serialization.path import posix_path
from GPErks.utils.array import get_minmax
from GPErks.perks.history_matching_mod import Wave
from GPErks.utils.sampling import Sampler

# extra imports
import pickle
import numpy as np


def history_matching(loadpath, savepath, xlabels, ylabels, bhm_ylabels, exp_mean_path, exp_std_path="", waveno=1, cutoff=3.0, maxno=1, n_tests=1e5, n_simuls=128, plot_dir="./figures", targetfile=""):
    #maxno = 1  # explained below
    ## The univariate GPE of each output feature will give for each point x_i a specific implausibility measure.
    ## With the current implausibility criterion, for each x_i we take the maximum implausibility across all the output
    ## features. With maxno=1, the maximum is calculated across all the output features (i.e., till the last worse
    ## implausibility measure). If maxno=2 --> till the previous-to-last worse implausibility measure and so on.
    ## With this criterion, the worse-performing emulator (the output feature which is the least well captured) will
    ## dominate the entire analysis and thus determine if a point is non-implausible or implausible
    
    # enforce reproducibility
    seed = DEFAULT_RANDOM_SEED
    set_seed(seed)

    print('Running history matching using features ')
    for ylab in bhm_ylabels:
        print(ylab)
    
    # extract output features used for HM from .txt file
    ylabels_dict = {key: idx for idx, key in enumerate(ylabels)}
    idx_feature = [ylabels_dict[key] for _, key in enumerate(bhm_ylabels)]
    if not idx_feature:
        raise Exception("There is no correspondence between the GPE outputs and the selected BHM features!")

    # load target data (mean and standard deviation) to match
    exp_mean = np.loadtxt(exp_mean_path + ".txt", dtype=float)
    if exp_std_path:
        exp_std = np.loadtxt(exp_std_path + ".txt", dtype=float)
    else:
        exp_std = 0.1 * exp_mean
        print("Using experimental std = 10% of experimental mean.")
    exp_var = exp_std**2

    # Load a pre-trained univariate Gaussian process emulator (GPE) for each output feature to match
    emulator_list = []
    for idx in idx_feature:
        print("\nLoading GPE(s)...\n")
        
        # current QoI folder
        idxfolder = loadpath  + "/output/" + str(idx)

        # load dataset
        file = open(idxfolder + "/dataset.pickle", "rb")
        dataset = pickle.load(file)
        file.close()

        # load trained emulator
        experiment = load_experiment_from_config_file(
            idxfolder + "/experiment.ini", dataset
        )
        # device = "cuda" if torch.cuda.is_available() else "cpu"
        device = "cpu"
        emulator = GPEmulator(experiment, device)

        best_model_file = posix_path(idxfolder, "best_model.pth")
        emulator.load_state(best_model_file)
        emulator.hyperparameters()
        #  NOTICE: GPEs must have been trained using GPErks library (https://github.com/stelong/GPErks)

        emulator_list.append(emulator)

    # get spanning range for each of parameters from training dataset
    I = get_minmax(dataset.X_train)

    # instantiate wave object
    W = Wave(
        emulator=emulator_list,
        Itrain=I,
        cutoff=cutoff,
        maxno=maxno,
        mean=exp_mean,
        var=exp_var,
    )
    
    # sample from current not-ruled-out-yet parameter space
    if waveno == 1:
        # initial wave is performed on a big Latin hypercube design using same parameter ranges of the training dataset
        sampler = Sampler(design="lhs", dim=I.shape[0], seed=seed)
        X = sampler.sample(
            int(n_tests),
            l_bounds=list(I[:, 0]),
            u_bounds=list(I[:, 1]),
        )
    else:
        # load TEST dataset from previous wave non-implausible region
        X = np.loadtxt(savepath + "/X_test_"+ str(waveno-1) + ".txt", dtype=float)
        n_tests = X.shape[0]

    # Run one iteration of HM: apply implausibility criterion to detect regions of non-implausible and of implausible points
    W.find_regions(X)
    W.print_stats()  # show statistics about the two obtained spaces
    
    # extract target parameter values from .txt file, if any
    if targetfile:
        target = list(np.genfromtxt(targetfile+".txt"))
    else:
        target = []

    # plots
    W.plot_wave_diagonal(xlabels=xlabels, filepath=plot_dir + "/diag_impl_"+ str(waveno) + ".png", target=target)  # plot current wave (impl. measure plot)
    W.plot_wave(xlabels=xlabels, display="impl", filepath=plot_dir + "/impl_"+ str(waveno) + ".png", target=target)  # plot current wave (impl. measure plot)
    W.plot_wave(xlabels=xlabels, display="var", filepath=plot_dir + "/var_"+ str(waveno) + ".png", target=target)  # check accuracy of GPEs for the current wave

    # To continue on next wave:

    # (0) save a copy of current wave before manipulating its internal structure (as it might be needed later for other purposes)
    W0 = W.copy()
    W0.save(savepath + "/wave_" + str(waveno) + ".json")

    # (1) from current non-implausible region, select points to be simulated and to be used as tests for next wave

    # if current available points are not enough, use 'cloud technique' to generate new NIMP points from existing ones
    n_avail_nimps = len(W0.nimp_idx)
    if n_tests + n_simuls > n_avail_nimps:
        print("\nUsing cloud technique to augment NIMP points.")
        n_total_points = n_tests + n_simuls
        W.augment_nimp(n_total_points)  
    
    # get the requested datasets
    X_simul, X_test = W.get_nimps(n_simuls)

    # save data to files.
    np.savetxt(savepath + "/X_simul_" + str(waveno) + ".txt", X_simul, fmt="%.6e")
    np.savetxt(savepath + "/X_test_" + str(waveno) + ".txt", X_test, fmt="%.6e")
    
    
    # Next steps:
    # (2) simulate selected points (X_simul)
    # (3) add simulated points and respective results to the training dataset used in the previous wave
    # (3) train GPEs on augmented training dataset
    # (4) start a new wave of history matching: the initial parameter space is no more a Latin hypercube design but is the non-implausible region obtained in the previous wave
    
    return