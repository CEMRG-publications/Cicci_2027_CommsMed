#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file is part of 'Cicci_2027_CommsMed'.

__author__ = Ludovica Cicci
__email__ = l.cicci@imperial.ac.uk
"""

# extra imports
import os
import shutil
import sys

# make config.py (in the parent folder) importable
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# dependencies
from config import libEP_DIR, OUTPUT_DIR, check_paths
from utils import sampling
from utils import IO


def main():
    check_paths("OUTPUT_DIR")

    basefolder = OUTPUT_DIR + "/KCL/History_Matching/QRScomplex"

    # name of .json file with default parameter values
    paramfilename = "simulation_parameters_fitting_biv" # "simulation_parameters_fitting_biv_best-fit_QRS"

    # name of .txt file(s) with input parameters
    labelfolder = libEP_DIR + "/data/xlabels"
    labelfilename = "xlabels_calibration_QRS" # "xlabels_calibration_Twave"

    # name of .json file(s) with parameter ranges
    rangefilename = "ranges_fitting_biv"

    # number of samples
    xlabels = IO.read_labels(labelfolder + "/" + labelfilename)
    N_samples = len(xlabels) * 20

    # directory where all samples and sampling data will be saved (usually, local HD)
    samplefolder = basefolder + "/samples_initial_space"

    # ---> START
    os.makedirs(basefolder, exist_ok=True)

    # check if output folder already exists to avoid overwriting
    if not os.path.exists(samplefolder):
        os.makedirs(samplefolder)
    elif os.listdir(samplefolder):
        raise Exception(
            "{} is non-empty. Change name or delete existing folder to avoid overwriting.".format(
                samplefolder
            )
        )

    # copy labels in samplefolder
    shutil.copy(libEP_DIR + "/data/xlabels/ToRORd_dynCl.txt", samplefolder)
    shutil.copy(libEP_DIR + "/data/xlabels/xlabels_empty.txt", samplefolder)
    shutil.copy(labelfolder + "/" + labelfilename + ".txt", samplefolder)

    # write fixed parameters
    sampling.save_fixed_parameters(
        libEP_DIR + "/data",
        paramfilename,
        samplefolder,
        labelfilename,
        samplefolder,
        N_samples=N_samples,
    )

    ## --------------- SCALAR SAMPLING --------------- ##
    sampling.sample_calibration_param(
        libEP_DIR + "/data",
        paramfilename,
        samplefolder,
        labelfilename,
        samplefolder,
        N_samples=N_samples,
        usr_range_json=libEP_DIR + "/data/ranges/" + rangefilename,
    )
    # <--- END


if __name__ == "__main__":
    main()
