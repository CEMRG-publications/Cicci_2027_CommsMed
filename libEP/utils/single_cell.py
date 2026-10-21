#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file is part of 'Cicci_2027_CommsMed'.

__author__ = Ludovica Cicci
__email__ = l.cicci@imperial.ac.uk
"""


# dependencies
from utils import IO
from utils import carp_func

# extra imports
import os
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA as PCA
from natsort import natsorted


def extract_im_par(samplefolder, labelfolder, imfilename, im="ToRORd_dynCl"):
    """Extract ionic model parameters from each sample.

    Parameters
    ----------
    samplefolder: str
        Path to folder with all sample(s) details.
    labelfolder: str
        Path to the input parameter datafile.
    imfilename: str
        Name of .txt file with labels of all ionic model parameters (without extension).
    im: str
        Name of ionic model. Default: "ToRORd_dynCl".
    """

    # set path to parameters .json file and to input labels .txt file
    imfile = labelfolder + "/" + imfilename

    # check folder structure
    print("\nChecking folder structure...")
    to_check = [imfile + ".txt"]
    for f in to_check:
        if not os.path.exists(f):
            raise Exception("Cannot find {}. Check path to file.".format(f))
        else:
            print(f + " found.")

    # get list of ionic model parameters
    im_par = list(IO.read_text(imfile))

    # find ionic model parameters for each sample
    for i in natsorted(os.listdir(samplefolder)):
        if not ".txt" in i:
            # current simulation paramfile
            paramfile = samplefolder + "/" + str(i) + "/param"
            if not os.path.exists(paramfile + ".json"):
                raise Exception("Cannot find {}.".format(paramfile + ".json"))

            # load dictionary with all current simulation parameters
            par_all_dct = IO.read_json(paramfile)
            par_all_keys = list(par_all_dct.keys())

            # biventricular heterogeneity
            biv_ch = []
            if any("_LV" in key for key in par_all_keys):
                biv_ch.append("_LV")
            if any("_RV" in key for key in par_all_keys):
                biv_ch.append("_RV")
            if not biv_ch:
                biv_ch = [""]

            # find ionic model parameters among scalar and spatially varying parameters (if any)
            im_par_scl = []
            im_par_grad = []
            for key in par_all_keys:
                if key in im_par:
                    im_par_scl.append(key)
                else:
                    ch = biv_ch[0]
                    if key.replace(ch + "_theta", "") in im_par:
                        im_par_grad.append(key.replace(ch + "_theta", ""))

            # convert theta, grad_rho and grad_z into apex/base values at the epi/endocardium
            #   1- apex_epi =: theta
            #   2- apex_endo = theta * grad_rho
            #   3- beta_epi = theta * grad_z
            #   4- beta_endo = theta * grad_rho * grad_z

            # initialisation
            imp_scl_dct = {}
            imp_grad_dct = {}
            for ch in biv_ch:
                imp_grad_dct[ch] = {}
                for p in ["apex_epi", "apex_endo", "base_epi", "base_endo"]:
                    imp_grad_dct[ch][p] = {}

            # save scalar ionic parameter values
            for imp in im_par_scl:
                imp_scl_dct[imp] = par_all_dct[imp]

            # save spatially varying ionic parameter values at apex/base epi/endo
            for imp in list(set(im_par_grad)):
                for ch in biv_ch:
                    imp_grad_dct[ch]["apex_epi"][imp] = par_all_dct[imp + ch + "_theta"]
                    imp_grad_dct[ch]["apex_endo"][imp] = (
                        par_all_dct[imp + ch + "_theta"]
                        * par_all_dct[imp + ch + "_grad_rho"]
                    )
                    imp_grad_dct[ch]["base_epi"][imp] = (
                        par_all_dct[imp + ch + "_theta"]
                        * par_all_dct[imp + ch + "_grad_z"]
                    )
                    imp_grad_dct[ch]["base_endo"][imp] = (
                        par_all_dct[imp + ch + "_theta"]
                        * par_all_dct[imp + ch + "_grad_rho"]
                        * par_all_dct[imp + ch + "_grad_z"]
                    )

            # assemble all parameters in a list
            imp_dct = {}
            if list(set(im_par_grad)):
                for ch in biv_ch:
                    imp_list = []
                    for p in ["apex_epi", "apex_endo", "base_epi", "base_endo"]:
                        # TODO: add "flags=..." only when using ToR-ORd models
                        if im == "ToRORd_dynCl":
                            imp_list.append(
                                "flags=ENDO,"
                                + ",".join(
                                    [
                                        imp + "=" + str(imp_scl_dct[imp])
                                        for imp in list(imp_scl_dct.keys())
                                    ]
                                )
                                + ","
                                + ",".join(
                                    [
                                        imp + "=" + str(imp_grad_dct[ch][p][imp])
                                        for imp in list(imp_grad_dct[ch][p].keys())
                                    ]
                                )
                            )
                        else:
                            imp_list.append(
                                ",".join(
                                    [
                                        imp + "=" + str(imp_scl_dct[imp])
                                        for imp in list(imp_scl_dct.keys())
                                    ]
                                )
                                + ","
                                + ",".join(
                                    [
                                        imp + "=" + str(imp_grad_dct[ch][p][imp])
                                        for imp in list(imp_grad_dct[ch][p].keys())
                                    ]
                                )
                            )
                    imp_dct[ch.replace("_", "")] = imp_list
            else:
                # TODO: add "flags=..." only when using ToR-ORd models
                if im == "ToRORd_dynCl":
                    imp_list = [
                        "flags=ENDO,"
                        + ",".join(
                            [
                                imp + "=" + str(imp_scl_dct[imp])
                                for imp in list(imp_scl_dct.keys())
                            ]
                        )
                    ]
                else:
                    imp_list = [
                        ",".join(
                            [
                                imp + "=" + str(imp_scl_dct[imp])
                                for imp in list(imp_scl_dct.keys())
                            ]
                        )
                    ]
                imp_dct["scl"] = imp_list

            # save output
            benchfolder = samplefolder + "/" + str(i) + "/BENCH"
            if not os.path.exists(benchfolder):
                os.mkdir(benchfolder)

            for _, key in enumerate(imp_dct):
                for j in range(len(imp_dct[key])):
                    if j == 0:
                        IO.write_string_text(
                            benchfolder + "/imp_par_" + key + "_" + str(j),
                            imp_dct[key][j],
                        )
                    elif not imp_dct[key][j] == imp_dct[key][j - 1]:
                        IO.write_string_text(
                            benchfolder + "/imp_par_" + key + "_" + str(j),
                            imp_dct[key][j],
                        )

    return


def check_AP_stats(
    APstatfolder,
    range_dct={
        "APD(n)": [180, 440],
        "VmMax": [10, 55],
        "Tri": [50, 150],
        "VmMin": [-95, -80],
    },  # based on doi.org/10.1016/j.cmpb.2023.107860
):
    """Check action potential statistics.

    Parameters
    ----------
    APstatsolder: str
        Path to folder with action potential statistics datafile(s).
    range_dct: dict
        Viable range for action potential features, as given the the AP_stats datafile. Default: {"APD(n)":[180,440],"VmMax":[10,55],"Tri":[50,150],"VmMin":[-95,-80]}.
    """

    # AP_stats header
    header = [
        "Beat",
        "Pre(P||*)",
        "lc",
        "APD(n)",
        "DI(n)",
        "DI(n-1)",
        "Tri",
        "VmMax",
        "VmMin",
    ]

    reject_dct = {}
    for filename in os.listdir(APstatfolder):
        if ".dat" in filename:
            # file with AP statistics
            APstatfile = APstatfolder + "/" + filename.replace(".dat", "")

            # load list with AP statistics for current sample
            AP_stat_dct = {}
            AP_stat = IO.read_data(APstatfile)

            if len(AP_stat):
                # extract statistics from last beat
                if np.ndim(AP_stat) > 1:
                    AP_stat = AP_stat[-1, :]

                for j, val in enumerate(AP_stat):
                    AP_stat_dct[header[j]] = val

                tmp_dct = {}
                for _, key in enumerate(range_dct):
                    if (
                        AP_stat_dct[key] < range_dct[key][0]
                        or AP_stat_dct[key] > range_dct[key][1]
                    ):
                        tmp_dct.update({key: AP_stat_dct[key]})

                if tmp_dct:
                    reject_dct[filename] = tmp_dct

            else:
                print("{} for sample is empty!".format(filename))
                reject_dct[filename] = None

    return list(reject_dct.keys())


def extract_AP(
    APstatfolder,
    APtracefilename,
    numstim=1,
    tstart=1.0,
    bcl=1000.0,
    cleandir=False,
    visualise=False,
):
    """Extract the action potential trace from last beat.

    Parameters
    ----------
    APstatfolder: str
        Path to folder with action potential statistics.
    APtracefilename: str
        Name of file containing the AP trace.
    numstim: int
        Number of stimuli. Default: 1.
    tstart: double.
        Start stimulation. Default: 1.
    bcl: double
        Basic cycle length. Default: 1000.
    dt_out: double
        Temporal output granularity. Default: 1.
    cleandir: bool
        If True, remove simulation output directory. Default: False.
    visualise: bool
        If True, plot action potential trace.
    """

    # file with AP trace
    APtracefile = APstatfolder + "/" + APtracefilename
    if not os.path.exists(APtracefile + "_AP_last_beat.txt"):
        if not os.path.exists(APtracefile + ".txt"):
            raise Exception(
                "AP trace {} has not been saved.".format(APtracefile + ".txt")
            )

        # read action potential trace from file
        AP = np.genfromtxt(
            APtracefile + ".txt",
            dtype=float,
        )

        # get initial timestep of last nstim beat(s)
        ind = np.where(AP[:, 0] >= int(tstart + (numstim - 1) * bcl))[0]
        lastAP = AP[ind[0] : -1, :]

        # save AP of last beat
        IO.write_array_text(APtracefile + "_AP_last_beat", lastAP)

        if visualise:
            AP = np.genfromtxt(
                APtracefile + "_AP_last_beat.txt",
                dtype=float,
            )
            plt.plot(AP[:, 0] - AP[0, 0], AP[:, 1], "-")
            plt.title("{}".format(APtracefilename))
            plt.show()

        # remove unnecessary files
        if cleandir:
            os.remove(APtracefile + ".txt")

    return


def samples_screening(
    CARP_INSTALL_DIR,
    samplefolder,
    ionicmodel,
    numstim=1,
    tstart=1.0,
    bcl=1000.0,
    dt_out=1.0,
    range_dct={
        "APD(n)": [180, 440],
        "VmMax": [10, 55],
        "Tri": [50, 150],
        "VmMin": [-95, -80],
    },  # based on doi.org/10.1016/j.cmpb.2023.107860
    savetrace=True,
    savestate=False,
    cleandir=False,
    visualise=False,
):
    """Run single-cell action potential simulation and check AP, for each sample.

    Parameters
    ----------
    CARP_INSTALL_DIR: str
        Path to CARPentry installation directory.
    samplefolder: str
        Path to folder with all sample(s) details.
    ionicmodel: str
        Name of the ionic model as given in "bench --list-imp".
    numstim: int
        Number of stimuli. Default: 1.
    tstart: double.
        Start stimulation. Default: 1.
    bcl: double
        Basic cycle length. Default: 1000.
    dt_out: double
        Temporal output granularity. Default: 1.
    range_dct: dict
        Viable range for action potential features, as given the the AP_stats datafile. Default: {"APD(n)":[180,440],"VmMax":[10,55],"Tri":[50,150],"VmMin":[-95,-80]}.
    savetrace: bool.
        If True, save trace of the action potential. Default: True.
    savestate: bool.
        If True, save initial state at the end of the pacing protocol. Default: False.
    cleandir: bool
        If True, remove simulation output directory. Default: False.
    visualise: bool
        If True, plot action potential trace. Defalt: False.
    """

    # find .txt file with accepted parameters. If already exists, rename it.
    readmefile = samplefolder + "/accepted_samples"
    if os.path.exists(readmefile + ".txt"):
        readmefile += "_overwritten"

    accepted_samples = "accepted_samples\n"

    # for each sample, run single cell simulations and decide to accept or reject the sample
    for i in natsorted(os.listdir(samplefolder)):
        if not ".txt" in i:
            print("Running single cell simulation for " + str(i) + "-th sample...")

            # current sample output directory
            SIMID = samplefolder + "/" + str(i)
            if not os.path.exists(SIMID):
                os.mkdir(SIMID)

            # folder with imp-combination text file(s)
            benchfolder = SIMID + "/BENCH"

            # initialise rejection flag
            reject = False

            # run single cell simulations for each combination of ionic model parameters
            for j in os.listdir(benchfolder):
                if not reject:
                    # AP stats output folder for current sample
                    APstatfolder = SIMID + "/Action_Potential"

                    # current imp paramfile
                    paramfile = (benchfolder + "/" + str(j)).replace(".txt", "")
                    if not os.path.exists(paramfile + ".txt"):
                        raise Exception("Cannot find {}.".format(paramfile + ".txt"))

                    carp_func.run_single_cell(
                        CARP_INSTALL_DIR,
                        APstatfolder,
                        paramfile,
                        ionicmodel,
                        numstim=numstim,
                        bcl=bcl,
                        dt_out=dt_out,
                        savetrace=savetrace,
                        savestate=savestate,
                    )

                    if savetrace:
                        extract_AP(
                            APstatfolder,
                            ionicmodel
                            + "_"
                            + os.path.basename(paramfile)
                            + "_AP_stats",
                            numstim=numstim,
                            tstart=tstart,
                            bcl=bcl,
                            cleandir=cleandir,
                            visualise=visualise,
                        )

                    reject = check_AP_stats(
                        APstatfolder,
                        range_dct=range_dct,
                    )

            if not reject:
                accepted_samples += str(i) + "\n"

    IO.write_string_text(readmefile, accepted_samples)

    return


def compute_qoi(samplefolder, outputfolder, pca_threshold=95, visualise=True):
    """Compute quantities of interest from action potential.

    Parameters
    ----------
    samplefolder: str
        Path to folder with all sample(s) details.
    outputfolder: str
        Path to the directory where the sample qoi will be saved.
    pca_threshold: int
        Percentage of output (i.e., Vm and dVm/dt) variance to be explained by the selected PCA modes. Default: 95 [%].
    visualise: bool
        If True, plot action potential trace. Default: False.
    """

    # check pca_threshold in (0,100]
    if pca_threshold <= 0 or pca_threshold > 100:
        raise Exception(
            "Select a viable percentage of output variance to be capture by selected PCA modes."
        )

    # IO
    if not os.path.exists(outputfolder):
        os.mkdir(outputfolder)

    Xfolder = outputfolder + "/X"
    if not os.path.exists(Xfolder):
        os.mkdir(Xfolder)

    yfolder = outputfolder + "/y"
    if not os.path.exists(yfolder):
        os.mkdir(yfolder)

    # save input labels and dataset
    xl_list = []
    X_list = []

    # find .txt file with accepted parameters
    readmefile = samplefolder + "/accepted_samples"
    if not os.path.exists(readmefile + ".txt"):
        readme = ["nan"] + [int(i) for i in os.listdir(samplefolder)]
    else:
        readme = IO.read_text(readmefile)

    # find ionic model parameters for each sample
    for i in natsorted(os.listdir(samplefolder)):
        if not ".txt" in i and str(i) in readme[1:]:
            # current simulation paramfile
            paramfile = samplefolder + "/" + str(i) + "/param"
            if not os.path.exists(paramfile + ".json"):
                raise Exception("Cannot find {}.".format(paramfile + ".json"))

            # load dictionary with all current simulation parameters
            par_imp_dct = IO.read_json(paramfile)

            # append current simulation scalar parameters
            xl_list.append([key for _, key in enumerate(par_imp_dct)])
            X_list.append([par_imp_dct[key] for _, key in enumerate(par_imp_dct)])

    xl_list = np.array(xl_list)
    X_list = np.array(X_list)

    # concatenate scalar and gradient input labels and dataset
    xlabels = np.array(xl_list)[0]
    X = np.array(X_list)
    np.savetxt(Xfolder + "/X.txt", X, fmt="%.6e")
    np.savetxt(Xfolder + "/xlabels.txt", xlabels, fmt="%s")

    # compute output(s) dataset(s): Action Potential (Vm), its derivative (dVm/dt), scalar outputs extracted from AP_stats.dat
    Vm_all = []
    dVmdt_all = []
    APD_all = []  # action potential duration
    DI_all = []  # diastolic interval
    Tri_all = []  # triangulation of the AP: APD90-APD40
    Vpeak_all = []  # Vm_max
    RMP_all = []  # Vm_min (~ resting membrane potential)

    for i in natsorted(os.listdir(samplefolder)):
        if not ".txt" in i and str(i) in readme[1:]:
            # extract last beat AP from current sample folder
            APfolder = samplefolder + "/" + str(i) + "/Action_Potential"
            if not os.path.exists(APfolder):
                print(
                    "Cannot compute AP-based quantities of interest. {} does not exist!".format(
                        APfolder
                    )
                )

            for filename in os.listdir(APfolder):
                if ".dat" in filename:
                    APstatfile = APfolder + "/" + filename.replace(".dat", "")

                    # load list with AP statistics for current sample
                    AP_stat = IO.read_data(APstatfile)

                    # extract statistics from last beat
                    if np.ndim(AP_stat) > 1:
                        AP_stat = AP_stat[-1, :]

                    # concatenate AP-based QoI (scalar) in a list
                    APD_all.append(AP_stat[3])
                    DI_all.append(AP_stat[4])
                    Tri_all.append(AP_stat[6])
                    Vpeak_all.append(AP_stat[7])
                    RMP_all.append(AP_stat[8])

                elif ".txt" in filename:
                    APfile = APfolder + "/" + filename.replace(".txt", "")

                    # load AP for current sample
                    AP = np.genfromtxt(APfile + ".txt", delimiter=" ")[:, 1]

                    # concatenate AP (array) in a list
                    Vm_all.append(AP)
                    dVmdt_all.append(np.gradient(AP))

    # save scalar outputs in text files
    IO.write_float_text(yfolder + "/y_APD", APD_all)
    IO.write_string_text(yfolder + "/ylabels_APD", ["APD"])

    IO.write_float_text(yfolder + "/y_DI", DI_all)
    IO.write_string_text(yfolder + "/ylabels_DI", ["DI"])

    IO.write_float_text(yfolder + "/y_Tri", Tri_all)
    IO.write_string_text(yfolder + "/ylabels_Tri", ["Tri"])

    IO.write_float_text(yfolder + "/y_Vpeak", Vpeak_all)
    IO.write_string_text(yfolder + "/ylabels_Vpeak", ["Vpeak"])

    IO.write_float_text(yfolder + "/y_RMP", RMP_all)
    IO.write_string_text(yfolder + "/ylabels_RMP", ["RMP"])

    # save time-dependent outputs dictionary
    np.savetxt(yfolder + "/Vm.txt", Vm_all, fmt="%.2f")
    np.savetxt(yfolder + "/dVmdt.txt", dVmdt_all, fmt="%.2f")

    # If 0 < n_components < 1 and svd_solver == "full", select the number of components such that the amount of variance that needs to be explained is greater than the percentage specified by n_components.
    qoi_dct = {"Vm": Vm_all, "dVmdt": dVmdt_all}
    for _, key in enumerate(qoi_dct):
        X_all = qoi_dct[key]

        if not visualise:
            # truncated matrix of principal component scores
            pca = PCA(n_components=pca_threshold / 100, svd_solver="full")
            Tl = pca.fit_transform(np.array(X_all))
            N_pca = pca.n_components_
            sv_l = pca.singular_values_

        else:
            # full matrix of principal component scores
            pca = PCA(n_components=0.99, svd_solver="full")
            T = pca.fit_transform(np.array(X_all))

            var_ratio = 100 * pca.explained_variance_ratio_
            cum_var_ratio = 100 * np.cumsum(pca.explained_variance_ratio_)
            N_pca = np.where(cum_var_ratio >= pca_threshold)[0][0] + 1
            Tl = T[:, :N_pca]
            sv_l = pca.singular_values_[:N_pca]

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
            plt.axvline(x=N_pca, color="blue", linestyle="--", label="95%")
            ax2.set_ylabel("Cumulative Explained Variance [%]", color="red")
            ax2.tick_params(axis="y", labelcolor="red", grid_alpha=0)
            ax2.set_ylim(
                ax1.set_ylim()[0] + 10 * int(int(var_ratio[0]) / 10),
                ax1.set_ylim()[1] + 10 * int(int(var_ratio[0]) / 10),
            )

            plt.title(str(pca_threshold) + " % explined variance", color="blue")
            plt.savefig(
                outputfolder
                + "/PCA_"
                + str(pca_threshold)
                + "_explained_variance_"
                + key
                + ".png"
            )

        print(
            "Number of modes capturing {}% of variation in {}: {}".format(
                pca_threshold, key, N_pca
            )
        )
        for l in range(np.size(Tl, 1)):
            np.savetxt(
                yfolder
                + "/y_PCA_"
                + key
                + "_"
                + str(l + 1)
                + "_of_"
                + str(np.size(Tl, 1))
                + ".txt",
                Tl[:, l],
                fmt="%.4f",
            )
            np.savetxt(
                yfolder
                + "/ylabels_PCA_"
                + key
                + "_"
                + str(l + 1)
                + "_of_"
                + str(np.size(Tl, 1))
                + ".txt",
                ["PCA_" + key + "_" + str(l + 1) + "_of_" + str(np.size(Tl, 1))],
                fmt="%s",
            )

            # save singular values in individual .txt files
            if not os.path.exists(yfolder + "/singular_values"):
                os.mkdir(yfolder + "/singular_values")
            np.savetxt(
                yfolder
                + "/singular_values/sv_PCA_"
                + key
                + "_"
                + str(l + 1)
                + "_of_"
                + str(np.size(Tl, 1))
                + ".txt",
                np.array(sv_l[l], ndmin=1, copy=False),
                fmt="%.4f",
            )
    return
