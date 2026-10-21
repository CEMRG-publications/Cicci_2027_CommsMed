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

    basefolder = OUTPUT_DIR + "/KCL"

    # name of .json file with default parameter values
    paramfilename = "simulation_parameters"

    # name of .txt file(s) with input parameters
    labelfolder = libEP_DIR + "/data/xlabels"
    labelfilename_scalar = "xlabels_scalar"
    labelfilename_gradient = "xlabels_gradient"

    # name of .json file(s) with parameter ranges
    rangefilename_scalar = "ranges_scalar"
    rangefilename_gradient = "ranges_gradient"

    # number of samples
    xlabels_scalar = IO.read_labels(labelfolder + "/" + labelfilename_scalar)
    xlabels_gradient = IO.read_labels(labelfolder + "/" + labelfilename_gradient)
    N_samples = (len(xlabels_scalar) + 3*len(xlabels_gradient)) * 10

    # directory where all samples and sampling data will be saved (usually, local HD)
    samplefolder = basefolder + "/samples_LHS"

    # ---> START
    os.makedirs(basefolder, exist_ok=True)

    # check if output folder already exists to avoid overwriting
    if not os.path.exists(samplefolder):
        os.mkdir(samplefolder)
    elif os.listdir(samplefolder):
        raise Exception(
            "{} is non-empty. Change name or delete existing folder to avoid overwriting.".format(
                samplefolder
            )
        )

    # copy labels in samplefolder
    shutil.copy(libEP_DIR + "/data/xlabels/ToRORd_dynCl.txt", samplefolder)
    shutil.copy(
        libEP_DIR + "/data/xlabels/" + labelfilename_scalar + ".txt", samplefolder
    )
    if labelfilename_gradient:
        shutil.copy(
            libEP_DIR + "/data/xlabels/" + labelfilename_gradient + ".txt",
            samplefolder,
        )

    # list of labelfiles
    labelfilename_list = [labelfilename_scalar]
    if labelfilename_gradient:
        labelfilename_list.append(labelfilename_gradient)

    # write fixed parameters
    sampling.save_fixed_parameters(
        libEP_DIR + "/data",
        paramfilename,
        samplefolder,
        labelfilename_list,
        samplefolder,
        N_samples=N_samples,
    )

    ## --------------- SCALAR SAMPLING --------------- ##
    sampling.sample_scalar_param(
        libEP_DIR + "/data",
        paramfilename,
        samplefolder,
        labelfilename_scalar,
        samplefolder,
        N_samples=N_samples,
        usr_range_json=libEP_DIR + "/data/ranges/" + rangefilename_scalar,
    )

    ## --------------- GRADIENT SAMPLING --------------- ##
    if labelfilename_gradient:
        sampling.sample_gradient_param(
            libEP_DIR + "/data",
            paramfilename,
            samplefolder,
            labelfilename_gradient,
            "ToRORd_dynCl",
            samplefolder,
            N_samples=N_samples,
            usr_range_json=libEP_DIR + "/data/ranges/" + rangefilename_scalar,
            usr_grad_json=libEP_DIR + "/data/ranges/" + rangefilename_gradient,
        )
    # <--- END


if __name__ == "__main__":
    main()
