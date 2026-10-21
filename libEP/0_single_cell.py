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
from utils import IO
from utils import sampling
from utils import single_cell as cell


def main():
    check_paths("CARP_INSTALL_DIR", "OUTPUT_DIR")

    basefolder = OUTPUT_DIR + "/ToRORd_dynCl"

    # name of ionic model
    ionicmodelname = "ToRORd_dynCl"

    # path to and name of .json file with default parameter values
    paramfolder = libEP_DIR + "/data"
    paramfilename = "single_cell_simulation_parameters"

    # directory with label .txt file(s) for input parameters, and corresponding filenames
    labelfolder = libEP_DIR + "/data/xlabels"
    labelfilename = "ToRORd_dynCl"

    # name of .json file(s) with parameter ranges
    rangefilename = "ranges_scalar"

    # number of samples
    xlabels = IO.read_labels(labelfolder + "/" + labelfilename)
    N_samples = len(xlabels) * 20

    # directory where all samples and sampling data will be saved (usually, local HD)
    samplefolder = basefolder + "/samples"

    # number of samples simulated at the same time
    num_proc = 8

    # ---> START
    # IO
    if not os.path.exists(basefolder):
        os.mkdir(basefolder)

    # check if output folder already exists to avoid overwriting
    if not os.path.exists(samplefolder):
        os.mkdir(samplefolder)

        # copy labels in samplefolder
        shutil.copy(labelfolder + "/" + labelfilename + ".txt", samplefolder)

        # write fixed parameters
        sampling.save_fixed_parameters(
            paramfolder,
            paramfilename,
            samplefolder,
            labelfilename,
            samplefolder,
            N_samples=N_samples,
        )

        # sampling
        sampling.sample_scalar_param(
            paramfolder,
            paramfilename,
            samplefolder,
            labelfilename,
            samplefolder,
            N_samples=N_samples,
            usr_range_json=libEP_DIR + "/data/ranges/" + rangefilename,
        )
    else:
        print("{} is not empty! Skipping samples generation.".format(samplefolder))

    # extract ionic model parameters
    cell.extract_im_par(samplefolder, samplefolder, labelfilename)

    # run single cell simulations with checks
    cell.samples_screening(
        CARP_INSTALL_DIR,
        samplefolder,
        ionicmodelname,
        numstim=35,
        bcl=1000,
        dt_out=0.1,
        savetrace=True,
        savestate=False,
        cleandir=True,
        visualise=False,
        num_proc=num_proc,
    )

    # compute scalar outputs
    cell.compute_qoi(samplefolder, basefolder + "/qoi", pca_threshold=90)
    # <--- END


if __name__ == "__main__":
    main()
