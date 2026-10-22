#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file is part of 'Cicci_2027_CommsMed'.

__author__ = Ludovica Cicci
__email__ = l.cicci@imperial.ac.uk
"""

# dependencies
from utils import IO
from utils import tools
from utils import plots

# extra imports
import os
import pandas as pd
import numpy as np
from natsort import natsorted


def extract_output_from_ecg(
    paramfolder,
    paramfilename,  # simulation_parameters
    labelfolder,
    labelfilename_scalar,
    labelfilename_gradient,
    samplefolder,
    simulationfolder,
    qoifolder,
    qoi_list=[
        "12leads",
        "12leads_truncated",
        "12leads_REP",
    ],
    tend=500,
    tend_truncated=150,
    threshold_PCA=0.95,
    PCAbasisfolder="",
    save_PCAbasis=True,
    visualise=True,
):
    """Extract output(s) of interest from 12-lead ECG.

    Parameters
    ----------
    paramfolder: str
        Path to the electrophysiology parameter datafile.
    paramfilename: str
        Name of .json file with names and reference values of all simulation parameters (without extension).
    labelfolder: str
        Path to the input parameter datafile.
    labelfilename_scalar: str
        Name of .txt file with lables of scalar parameters of interest (without extension).
    labelfilename_gradient: list
        Name of .txt file with lables of spatially varying parameters of interest (without extension), if any. Otherwise, use empty string "".
    samplefolder: str
        Path to folder with all sample(s) details.
    simulationfolder: str
        Path to folder with all simulation directories.
    qoifolder: str
        Path to the directory where the output(s) of interest will be saved.
    qoi_list: list
        List of outputs to extract. Default: ["12leads", "12leads_truncated", "12leads_REP"].
    tend: float
        Defines the point in time when the simulation stops. Default: 500 [ms].
    tend_truncated: float
        Defines the point in time to compute (truncated) ECG. Default: 150 [ms].
    threshold_PCA : float
        Ratio of variance to be explained by N principal components. Must be in [0,1].
    PCAbasisfolder: str
        Path to folder with existing PCA basis (and empirical ECG mean) on which to project the ECG data, if any. Otherwise, use empty string "" to use current ECG data to compute the basis. Default: "".
    save_PCAbasis : bool
        If True, save PCA basis. Default: True.
    visualise : bool
        If True, plot bar plot of explained varience against number of principal components. Default: True.
    """

    if labelfilename_gradient:
        labelfilename_list = [labelfilename_scalar, labelfilename_gradient]
    else:
        labelfilename_list = [labelfilename_scalar]

    if type(qoi_list) == str:
        qoi_list = [qoi_list]

    # IO
    if not os.path.exists(qoifolder):
        os.mkdir(qoifolder)

    # check qoi in predefine list
    for q in qoi_list:
        if not q in [
            "12leads",
            "12leads_truncated",
            "12leads_REP",
        ]:
            raise Exception(
                "{} cannot be computed. Choose QoI from {}".format(
                    q,
                    [
                        "12leads",
                        "12leads_truncated",
                        "12leads_REP",
                    ],
                )
            )

    # set path to parameters .json file and to input labels .txt file
    paramfile = paramfolder + "/" + paramfilename

    labelfile_list = [
        labelfolder + "/" + labelfilename for labelfilename in labelfilename_list
    ]

    # check folder structure
    print("\nChecking folder structure...")
    to_check = [paramfile + ".json"] + [
        labelfile + ".txt" for labelfile in labelfile_list
    ]
    for f in to_check:
        if not os.path.exists(f):
            raise Exception("Cannot find {}. Check path to file.".format(f))
        else:
            print(f + " found.")

    Xfolder = qoifolder + "/X"
    if not os.path.exists(Xfolder):
        os.mkdir(Xfolder)

    yfolder = qoifolder + "/y"
    if not os.path.exists(yfolder):
        os.mkdir(yfolder)

    # load dictionary with all simulation parameters (including those that are NOT considered as inputs)
    par_all_dct = IO.read_json(paramfile)

    # find parameters of interest as given in labelfile
    par_dct = {}
    for l in range(len(labelfile_list)):
        for key in IO.read_labels(labelfile_list[l]):
            if l == 0:
                par_dct[key] = par_all_dct[key]
            else:
                if key in list(par_dct.keys()):
                    raise Exception(
                        "{} can't be both scalar and spatially varying. Remove it from one of the label files.".format(
                            key
                        )
                    )
                par_dct[key + "_theta"] = par_all_dct[key]
                par_dct[key + "_grad_rho"] = 1
                par_dct[key + "_grad_z"] = 1

    xlabels = list(par_dct.keys())

    # keep only successful simulations
    successful_simulations = []
    successful_samples = ""
    for i in natsorted(os.listdir(simulationfolder)):
        if "EP" in i:
            bspdatfile = simulationfolder + "/" + i + "/lead_field_egms.dat"
            if os.path.exists(bspdatfile) and os.path.getsize(bspdatfile):
                successful_simulations.append(int(i.replace("SIM_EP_", "")))
                successful_samples += str(i) + "\n"
            else:
                print("{} was unsuccessful!".format(i))

    success_str = "Total number of successful simulations: {}\n".format(
        str(len(successful_simulations))
    )
    print("\n\n" + success_str + "\n")
    IO.write_string_text(qoifolder + "/README", success_str + successful_samples)

    # extract input(s) for successful simulations and compute output(s) dataset(s)
    X = []
    QoI_dct = {}
    for q in qoi_list:
        QoI_dct[q] = []

    successful_simulations = np.sort(successful_simulations)
    for i in successful_simulations:
        # current simulation paramfile
        paramfile = samplefolder + "/" + str(i) + "/param"
        if not os.path.exists(paramfile + ".json"):
            raise Exception("Cannot find {}.".format(paramfile + ".json"))

        # get all parameters from current sample
        par_dct = IO.read_json(paramfile)

        # append current simulation input parameters
        X.append([par_dct[k] for k in xlabels])

        # read BSP from datafile: array (tend/spacedt) x N_electrodes
        # NOTE: "lead_field_egms" is the name of the CARPentry output file
        bspfile = simulationfolder + "/SIM_EP_" + str(i) + "/lead_field_egms"
        BSP = IO.read_data(bspfile)

        # convert 10-lead BSP into 12-lead ECG signals
        ecgfile = qoifolder + "/ecg_data"
        tools.bsp_to_ecg(ecgfile, BSP)

        # process ECG and VCG
        ecg = IO.read_ecg_dat(ecgfile + ".dat")

        if "12leads" in qoi_list:
            # concatenate ECGs in a list
            ecg_signal = np.array(ecg).reshape(1, -1, order="F").tolist()[0]
            QoI_dct["12leads"].append(ecg_signal)

        if "12leads_truncated" in qoi_list:
            # concatenate truncated ECGs in a list
            tf = int((tend_truncated) / 2)
            ecg_signal = (
                np.array(ecg[:tf]).reshape(1, -1, order="F").tolist()[0]
            )
            QoI_dct["12leads_truncated"].append(ecg_signal)

        if "12leads_REP" in qoi_list:
            # get activation time
            LATfile = simulationfolder + "/SIM_EP_" + str(i) + "/init_acts_ACT-thresh"
            LAT = IO.read_data(LATfile)
            TAT = max(LAT) - min(LAT)

            # find time indices at the end of the activation
            time = np.linspace(0, tend, len(np.array(ecg.index)))
            tend_act = np.where(time > np.min(TAT))[0]

            # extract ECG, separate ACT (activation) and REP (repolarisation)
            # NOTE: REP/ACT points are set equal to 0, respectively.
            ECG_array = np.array(ecg)
            M = np.ones(np.shape(ECG_array))
            M[tend_act, :] = 0

            # concatenate ECGs in a list
            ecg_signal_rep = (
                (ECG_array * (1 - M)).reshape(1, -1, order="F").tolist()[0]
            )
            QoI_dct["12leads_REP"].append(ecg_signal_rep)

        print("\n")

    # save input matrix
    IO.write_array_text(Xfolder + "/X", np.array(X))
    IO.write_string_text(Xfolder + "/xlabels", [x + "\n" for x in xlabels])

    for q in qoi_list:
        if "leads" in q:
            # save ECG track in text file
            IO.write_array_text(qoifolder + "/ECG_" + q, QoI_dct[q])

            if PCAbasisfolder:
                # load PCA basis WL: Nt x L
                WL = np.genfromtxt(PCAbasisfolder + "/WL_" + q + ".txt").T

                # project new ECG dataset onto PCA basis --> TL: Nsamples x L
                Tl = np.matmul(np.array(QoI_dct[q]), WL)

            else:
                Tl, WL = tools.compute_PCA(
                    QoI_dct[q],
                    q,
                    threshold_PCA=threshold_PCA,
                    visualise=visualise,
                    plotfolder=qoifolder,
                )

                # add traning ECG dataset empirical mean
                emp_mean = np.mean(np.array(QoI_dct[q]), axis=0)
                Tl = Tl + np.matmul(emp_mean, WL.T)

                if save_PCAbasis:
                    # save PCA-basis matrix
                    IO.write_array_text(qoifolder + "/WL_" + q, np.array(WL))

            # save output matrix for each QoI
            for l in range(np.size(Tl, 1)):
                IO.write_float_text(
                    yfolder
                    + "/y_PCA_"
                    + q
                    + "_"
                    + str(l + 1)
                    + "_of_"
                    + str(np.size(Tl, 1)),
                    Tl[:, l],
                )
                IO.write_string_text(
                    yfolder
                    + "/ylabels_PCA_"
                    + q
                    + "_"
                    + str(l + 1)
                    + "_of_"
                    + str(np.size(Tl, 1)),
                    "PCA_" + q + "_" + str(l + 1) + "_of_" + str(np.size(Tl, 1)),
                )
        else:
            # save output matrix for each QoI
            QoI_dct[q] = (np.matrix(QoI_dct[q]).reshape(1, -1)).T
            IO.write_float_text(yfolder + "/y_" + q, QoI_dct[q])
            IO.write_string_text(yfolder + "/ylabels_" + q, q)

    os.remove(qoifolder + "/ecg_data.dat")

    return


###############################################################################


def extract_output_from_clinical_ecg(
    PCAbasisfolder,
    ecgdatafolder,  # /path/to/target
    ecgdatafilename,
    qoifolder,
    qoi_list=[
        "12leads",
        "12leads_truncated",
        "12leads_REP",
    ],
    tend_truncated=150
):
    """Extract output(s) of interest from 12-lead ECG."""

    if type(qoi_list) == str:
        qoi_list = [qoi_list]

    # IO
    if not os.path.exists(qoifolder):
        os.mkdir(qoifolder)

    # check qoi in predefine list
    for q in qoi_list:
        if not q in [
            "12leads",
            "12leads_truncated",
            "12leads_REP",
        ]:
            raise Exception(
                "{} cannot be computed. Choose QoI from {}".format(
                    q,
                    [
                        "12leads",
                        "12leads_truncated",
                        "12leads_REP",
                    ],
                )
            )

    yfolder = qoifolder + "/y"
    if not os.path.exists(yfolder):
        os.mkdir(yfolder)

    # extract input(s) for successful simulations and compute output(s) dataset(s)
    QoI_dct = {}
    for q in qoi_list:
        QoI_dct[q] = []

    # process ECG and VCG
    ecg = IO.read_ecg_dat(ecgdatafolder + "/" + ecgdatafilename + ".dat")

    if "12leads" in qoi_list:
        # concatenate ECGs in a list
        ecg_signal = np.array(ecg).reshape(1, -1, order="F").tolist()[0]
        QoI_dct["12leads"].append(ecg_signal)

    if "12leads_truncated" in qoi_list:
        # concatenate truncated ECGs in a list
        tf = int((tend_truncated) / 2)
        ecg_signal = np.array(ecg[:tf]).reshape(1, -1, order="F").tolist()[0]
        QoI_dct["12leads_truncated"].append(ecg_signal)

    if "12leads_REP" in qoi_list:
        tend_act = int((tend_truncated) / 2)

        # extract ECG, separate ACT (activation) and REP (repolarisation)
        # NOTE: REP/ACT points are set equal to 0, respectively.
        ECG_array = np.array(ecg)
        M = np.ones(np.shape(ECG_array))
        M[tend_act:, :] = 0

        # concatenate ECGs in a list
        ecg_signal_rep = (ECG_array * (1 - M)).reshape(1, -1, order="F").tolist()[0]
        QoI_dct["12leads_REP"].append(ecg_signal_rep)

    print("\n")

    for q in qoi_list:
        # save ECG track in text file
        IO.write_array_text(qoifolder + "/ECG_" + q, QoI_dct[q])

        # load PCA basis WL: Nt x L
        WL = np.genfromtxt(PCAbasisfolder + "/WL_" + q + ".txt").T

        # project new ECG dataset onto PCA basis --> TL: Nsamples x L
        Tl = np.matmul(np.array(QoI_dct[q]), WL)

        # save output matrix for each QoI
        for l in range(np.size(Tl, 1)):
            IO.write_float_text(
                yfolder
                + "/y_PCA_"
                + q
                + "_"
                + str(l + 1)
                + "_of_"
                + str(np.size(Tl, 1)),
                Tl[:, l],
            )
            IO.write_string_text(
                yfolder
                + "/ylabels_PCA_"
                + q
                + "_"
                + str(l + 1)
                + "_of_"
                + str(np.size(Tl, 1)),
                "PCA_" + q + "_" + str(l + 1) + "_of_" + str(np.size(Tl, 1)),
            )

    return


###############################################################################


def fit_linear_function(X, Y, xlabels):
    """Linear function interpolation of given points.

    Parameters
    ----------
    X: array
        Matrix of interpolation points. Must be defined accordingly to "linear_fitting_sampling.py" (see also "linear_fitting_screening.py").
    Y: array
        Matrix of interpolation values.
    xlabels: list
        Input labels.
    """

    # compute coefficients
    Ak_dct = {}
    for k, key in enumerate(xlabels):
        print("-------------------------------")
        print("-------------------------------")
        print("Assemble linear system Ha_l = v_l for {}".format(key))

        #      H       *  a_l  =     v_l
        #  ----------     ---     -----------
        # | min(X) 1 |   | a |   | Y(min(X)) |
        # |          | * |   | = |           |
        # | max(X) 1 |   | b |   | Y(max(X)) |

        # assemble H
        H = np.zeros((2, 2))
        H[0, :] = [X[2 * k, k], 1]
        H[1, :] = [X[2 * k + 1, k], 1]

        # assamble v_l
        # v_l 2 x 1 --> VL[:,l] = v_l, for l=1:L --> VL 2 x L
        VL = np.zeros((2, Y.shape[1]))
        VL[0, :] = Y[2 * k, :]
        VL[1, :] = Y[2 * k + 1, :]

        print("-------------------------------")
        print("Solve linear system Ha_l = v_l")

        # a_l 2 x 1 --> Ak[:,l] = a_l, for l=1:L --> Ak 2 x L
        # Ak(:,l) = [a_{k}^(l) b_{k}^(l)].T
        Ak = np.linalg.solve(H, VL)
        Ak_dct[xlabels[k]] = Ak

    return Ak_dct


def parameter_influence_on_PCAmodes(paramname, paramvar, threshold=0.1):
    """Auxiliary function to save boolean matrix of parameter influence on each PCA-modes.

    Parameters
    ----------
    paramname: str
        Name of current input parameter.
    paramvar: numpy array
        Variance of current parameter, for each mode.
    """

    N_pca = len(paramvar)

    # write summary
    headers = ["parameter"]
    data = [paramname]
    for i in range(N_pca):
        headers.append("mode{}".format(i + 1))
        if paramvar[i] > threshold:
            data.append("1")
        else:
            data.extend("0")

    return headers, data


def ranking_PCAmode_variance(
    samplefolder,
    PCAbasisfolder,
    qoifolder,
    plotfolder,
    WLname="WL_12leads",
    ECGname="ECG_12leads",
    cutoff=0.9,
    threshold=0.1,
    visualise=False,
):
    """Rank parameters according to the cumulative variance.

    Parameters
    ----------
    samplefolder: str
        Path to folder with all sample(s) details
    PCAbasisfolder: str
        Path to folder with existing PCA basis on which to project the ECG data.
    qoifolder: str
        Path to the directory with the output(s) of interest.
    plotfolder: str
        Name of the directory where the figures will be saved, if any. Default: "./".
    WLname: str
        Name of the PCAbasis. Default: "WL_12leads".
    ECGname: str
        Name of the ECG file. Default: "ECG_12leads".
    cutoff: float
        Cut-off value for most important parameters (i.e., those accounting for cutoff*100 % of PCAmode variance). Must be in [0,1]. Default: 0.9.
    threshold: float
        Threshold value to screen whether a specific input parameter has impact on each PCA-mode. Default: 0.1.
    visualise: bool
        If True, plot and save PCA-modes, for each parameter. Default: False.
    """

    ## compute PCA-based ranking criteria

    # load PCA basis WL: Nt x L
    WL = np.genfromtxt(PCAbasisfolder + "/" + WLname + ".txt").T
    L = WL.shape[1]

    print("-------------------------------")
    print("Extract interpolation points.")

    THETA = np.genfromtxt(samplefolder + "/ThetaArray.txt")

    print("-------------------------------")
    print("Interpolate PC-based function for each input.")

    # load ECG matrix (2*Nparam) x Nt --> ECG(h,:) = ecg(X(h,:)), h = 0:2*Nparam-1
    ECG = np.genfromtxt(qoifolder + "/" + ECGname + ".txt")

    # get input labels as used in linear_fitting_sampling
    xlabels = IO.read_text(samplefolder + "/ThetaArray_labels")

    # fit linear function
    Ak = fit_linear_function(THETA, np.matmul(ECG, WL), xlabels)

    # compute variance for each PC-component
    var_dct = {}
    var_array = np.zeros((len(xlabels), L))
    for k in range(len(xlabels)):
        for l in range(L):
            var_array[k, l] = (
                Ak[xlabels[k]][0, l] * (THETA[2 * k + 1, k] - THETA[2 * k, k])
            ) ** 2 / 12
        var_dct[xlabels[k]] = var_array[k, :]

    # save variances
    var_df = pd.DataFrame.from_dict(var_dct)
    var_df.to_csv(qoifolder + "/Var_uniform_variance.csv")

    ###########################################################################
    # normalise absolute variance for each input
    modevar_norm = var_array / np.sum(var_array, axis=0)
    
    # plot ranked normalised mode-based variance, with cut-off
    if visualise:
        for l in range(np.size(modevar_norm, 1)):
            if not os.path.exists(plotfolder + "/ModeVar"):
                os.mkdir(plotfolder + "/ModeVar")
            plots.plot_ranked_cutoff(
                xlabels,
                modevar_norm[:, l],
                cutoff=cutoff,
                figname=plotfolder + "/ModeVar/ranked_variance_mode_" + str(l + 1),
            )
    ###########################################################################

    # compute cumulative variance
    cumvar = np.sum(var_array, axis=1)

    # plot cumulative variance
    plots.plot_cumvar(xlabels, var_array, plotfolder, L_max=L)

    # normalise cumulative variance
    cumvar_norm = cumvar / np.sum(cumvar)

    # rank inputs according to cumulative variance (descending order)
    order_list = np.argsort(-cumvar_norm)
    xlabels_ranked = np.array(xlabels)[order_list].tolist()
    cumvar_norm_ranked = cumvar_norm[order_list]

    # save corresponding ranked inputs and normalised cumulative variance
    IO.write_array_text(
        qoifolder + "/ranked_inputs",
        xlabels_ranked,
        fmt="%s",
    )
    IO.write_float_text(
        qoifolder + "/ranked_normalised_cumulative_variance",
        cumvar_norm_ranked,
    )

    # save param impact on PC-modes
    data = []
    for j in range(len(xlabels)):
        headers, data_j = parameter_influence_on_PCAmodes(
            xlabels[j], var_array[j, :], threshold
        )
        data.append(data_j)
    df = pd.DataFrame(data, columns=headers)
    df.to_csv(qoifolder + "/modes_of_variation.csv", index=False)

    # plot PC-modes
    if visualise:
        if not os.path.exists(plotfolder + "/PCAmodes"):
            os.mkdir(plotfolder + "/PCAmodes")

        X = np.genfromtxt(PCAbasisfolder + "/" + ECGname + ".txt")
        for j in range(len(xlabels)):
            plots.plot_PCA_modes(
                X,
                xlabels[j],
                var_array[j, :],
                figname=plotfolder
                + "/PCAmodes/modes_of_variation_for_{}".format(xlabels[j]),
            )

    return cumvar
