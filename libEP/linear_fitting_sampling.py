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

    meshname = "KCL1"  # NOTE: mesh basename
    basefolder = OUTPUT_DIR + "/" + meshname

    # name of .json file with default parameter values
    paramfilename = "simulation_parameters_fitting_biv"

    # name of .txt file(s) with input parameters
    labelfilename = "xlabels_fitting_biv"

    # name of .json file(s) with parameter ranges
    rangefilename = "ranges_fitting_biv"

    # directory where all samples and sampling data will be saved (usually, local HD)
    samplefolder = basefolder + "/samples_fitting"

    # ---> START

    # IO
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
    shutil.copy(libEP_DIR + "/data/xlabels/" + labelfilename + ".txt", samplefolder)

    # get number of (scalar) parameters
    N_param = len(IO.read_labels(libEP_DIR + "/data/xlabels/" + labelfilename))

    # write fixed parameters
    sampling.save_fixed_parameters(
        libEP_DIR + "/data",
        paramfilename,
        samplefolder,
        labelfilename,
        samplefolder,
        N_samples=2 * N_param,
    )

    # scalar sampling
    sampling.sample_fitting(
        libEP_DIR + "/data",
        paramfilename,
        samplefolder,
        labelfilename,
        samplefolder,
        usr_range_json=libEP_DIR + "/data/ranges/" + rangefilename,
    )
    # <--- END


if __name__ == "__main__":
    main()
