#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file is part of 'Cicci_2027_CommsMed'.

__author__ = Ludovica Cicci
__email__ = l.cicci@imperial.ac.uk
"""

# extra imports
import numpy as np
import os
import sys

# make config.py (in the parent folder) importable
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# dependencies
from config import OUTPUT_DIR, check_paths
from utils import IO
from utils import global_sensitivity_analysis as gsa
from utils import plots

"""
Note: the following code runs a preliminary screening of the input parameters based on the assumption that inputs with negligible first-order effects have also negligible higher-order effects.
It uses a linear function for each input, built by interpolating the following points:
    - min value;
    - max value.
The matrix H of interpolation points is automatically built by running "linear_fitting_sampling.py" and has dimensions (2*Nparam) x Nparam.
"""


def main():
    check_paths("OUTPUT_DIR")

    meshname = "KCL1"  # NOTE: mesh basename
    basefolder = OUTPUT_DIR + "/" + meshname

    # current ECG segment
    ECGsegment = "QRScomplex"

    # directory where samples are saved
    samplefolder = basefolder + "/samples_fitting"

    # cut-off for model parameters [for HM]
    cutoff = 0.95

    # If True, plot ECGs and variance comparision
    visualise = False

    # threshold used to compute the PCA basis
    threshold_PCA = 0.99

    # ---> START

    thPCA = str(int(100 * threshold_PCA))

    if "QRScomplex" in ECGsegment:
        qoifolder = basefolder + "/qoi_fitting_150ms"

        # basename for 12-lead ECG signal and PCA-basis
        PCAbasisfolder = basefolder + "/PCAbasis_" + meshname + "_150ms_" + thPCA + "%"
        ECGname = "ECG_12leads_truncated"
        WLname = "WL_12leads_truncated"

    elif "Twave" in ECGsegment:
        qoifolder = basefolder + "/qoi_fitting_REP"

        # basename for 12-lead ECG signal and PCA-basis
        PCAbasisfolder = basefolder + "/PCAbasis_" + meshname + "_REP_" + thPCA + "%"
        ECGname = "ECG_12leads_REP"
        WLname = "WL_12leads_REP"

    else:
        qoifolder = basefolder + "/qoi_fitting"

        # basename for 12-lead ECG signal and PCA-basis
        PCAbasisfolder = basefolder + "/PCAbasis_" + meshname + "_" + thPCA + "%"
        ECGname = "ECG_12leads"
        WLname = "WL_12leads"

    # IO
    plotfolder = qoifolder + "/lsa"
    os.makedirs(plotfolder, exist_ok=True)

    # get input labels as used in linear_fitting_sampling
    xlabels = IO.read_text(samplefolder + "/ThetaArray_labels")

    # load ECG matrix (2*Nparam) x Nt --> ECG(h,:) = ecg(X(h,:)), h = 0:2*Nparam-1
    ECG = np.genfromtxt(qoifolder + "/" + ECGname + ".txt")

    # plot ECG
    if visualise:
        plots.plot_ecg_linear_fitting(xlabels, ECG)

    ## compute PCA-based ranking criterium
    cumvar = gsa.ranking_PCAmode_variance(
        samplefolder,
        PCAbasisfolder,
        qoifolder,
        plotfolder,
        WLname=WLname,
        ECGname=ECGname,
        cutoff=cutoff,
        visualise=visualise,
    )

    ## rank inputs in descending order of normalised cumulative variance
    cumvar_norm = cumvar / np.sum(cumvar)
    order_list = np.argsort(-cumvar_norm)
    xlabels_ranked = np.array(xlabels)[order_list].tolist()

    # save corresponding ranked inputs
    IO.write_array_text(
        qoifolder + "/final_ranking_cumvar",
        xlabels_ranked,
        fmt="%s",
    )

    # plot final ranking with cutoff
    plots.plot_final_ranking(
        cumvar_norm,
        xlabels,
        {"cumvar": cumvar_norm},
        cutoff=cutoff,
        barcol=["teal"],
        figname=plotfolder + "/final_ranking_cumvar",
    )
    # <--- END


if __name__ == "__main__":
    main()
