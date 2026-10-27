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
from config import CARP_INSTALL_DIR, libEP_DIR, OUTPUT_DIR, check_paths
from utils import sampling
from utils import IO
from utils import single_cell as cell


def main():
    check_paths("OUTPUT_DIR")

    meshname = "KCL1"  # NOTE: mesh basename
    basefolder = OUTPUT_DIR + "/" + meshname

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

    # If True, run a preliminary single-cell check of the samples, to discard those that are highly likely to fail (needs CARPentry).
    # If the sample folder already exists, the sampling is skipped and only the check is run (simulations already completed are not run again).
    check_single_cell = False

    # number of samples simulated at the same time in the single-cell check
    num_proc = 8

    # ---> START
    if check_single_cell:
        check_paths("CARP_INSTALL_DIR")

    os.makedirs(basefolder, exist_ok=True)

    # check if output folder already exists to avoid overwriting
    if not os.path.exists(samplefolder):
        os.mkdir(samplefolder)
    elif os.listdir(samplefolder):
        if not check_single_cell:
            raise Exception(
                "{} is non-empty. Change name or delete existing folder to avoid overwriting.".format(
                    samplefolder
                )
            )
        print(
            "{} is non-empty. Skipping samples generation and running the single-cell check only.".format(
                samplefolder
            )
        )
        cell.check_whole_organ_samples(
            CARP_INSTALL_DIR, samplefolder, num_proc=num_proc
        )
        return

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

    ## --------------- SINGLE-CELL CHECK --------------- ##
    if check_single_cell:
        cell.check_whole_organ_samples(
            CARP_INSTALL_DIR, samplefolder, num_proc=num_proc
        )
    # <--- END


if __name__ == "__main__":
    main()
