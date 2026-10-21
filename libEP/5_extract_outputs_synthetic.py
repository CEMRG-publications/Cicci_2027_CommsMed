#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file is part of 'Cicci_2027_CommsMed'.

__author__ = Ludovica Cicci
__email__ = l.cicci@imperial.ac.uk
"""

# extra imports
import os
import sys

# make config.py (in the parent folder) importable
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# dependencies
from config import libEP_DIR, OUTPUT_DIR, check_paths
from utils import global_sensitivity_analysis as gsa

def main():
    check_paths("OUTPUT_DIR")

    basefolder = OUTPUT_DIR + "/KCL"
    meshname = "KCL1"

    # current ECG segment
    ECGsegment = "QRScomplex"

    # name of .json file with default parameter values
    paramfilename = "simulation_parameters"

    # name of .txt file(s) with input parameters
    labelfilename_scalar = "xlabels_scalar"
    labelfilename_gradient = "xlabels_gradient"  # if any, otherwise use empty string "".

    # directory where samples are saved
    samplefolder = basefolder + "/samples"

    # directory with simulation output file(s)
    simulationfolder = basefolder + "/simulations"

    # directory where all qoi will be saved (usually, local HD)
    qoifolder = basefolder + "/qoi_" + ECGsegment

    # PCA threshold
    threshold_PCA = 0.99

    if "QRScomplex" in ECGsegment:
        # directory where PCA basis and training ECG dataset mean are saved, if any. Otherwise use empty string "".
        PCAbasisfolder = basefolder + "/PCAbasis_" + meshname + "_150ms_{}%".format(
            str(int(100 * threshold_PCA))
        )
        qoi_list = ["12leads_truncated"]
    elif "Twave" in ECGsegment:
        # directory where PCA basis and training ECG dataset mean are saved, if any. Otherwise use empty string "".
        PCAbasisfolder = basefolder + "/PCAbasis_" + meshname + "_REP_{}%".format(
            str(int(100 * threshold_PCA))
        )
        qoi_list = ["12leads_REP"]
    else:
        # directory where PCA basis and training ECG dataset mean are saved, if any. Otherwise use empty string "".
        PCAbasisfolder = basefolder + "/PCAbasis_" + meshname + "_{}%".format(
            str(int(100 * threshold_PCA))
        )
        qoi_list = ["12leads"]

    # ---> START
    gsa.extract_output_from_ecg(
        libEP_DIR + "/data",
        paramfilename,
        samplefolder,
        labelfilename_scalar,
        labelfilename_gradient,
        samplefolder,
        simulationfolder,
        qoifolder=qoifolder,
        qoi_list=qoi_list,
        tend=500,
        tend_truncated=150,
        threshold_PCA=threshold_PCA,
        PCAbasisfolder=PCAbasisfolder,
        save_PCAbasis=True,
        visualise=True,
    )
    # <--- END


if __name__ == "__main__":
    main()
