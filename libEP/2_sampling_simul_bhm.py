#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file is part of 'Cicci_2027_CommsMed'.

__author__ = Ludovica Cicci
__email__ = l.cicci@imperial.ac.uk
"""

# other imports
import os
import shutil
import sys

# make config.py (in the parent folder) importable
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# dependencies
from config import libEP_DIR, OUTPUT_DIR, check_paths
from utils import sampling


def main():
    check_paths("OUTPUT_DIR")

    # current wave ID
    waveno = 1

    # folder of the current history matching run
    basefolder = OUTPUT_DIR + "/KCL/History_Matching/QRScomplex"

    # name of .json file with default parameter values
    paramfilename = "simulation_parameters_fitting_biv" # "simulation_parameters_fitting_biv_best-fit_QRS"

    # directory with label .txt file(s) for GPE/GSA/BHM input parameters
    labelfolder = libEP_DIR + "/data/xlabels"
    labelfilename = "xlabels_calibration_QRS" # "xlabels_calibration_Twave"

    # HM new samples for current wave (genereted in previous wave)
    simulfolder = basefolder + "/output"

    # directory where all samples and sampling data will be saved (usually, local HD)
    samplefolder = basefolder + "/samples_wave_" + str(waveno)

    # ---> START
    os.makedirs(samplefolder, exist_ok=True)

    sampling.assemble_simulation_param(
        libEP_DIR + "/data",
        paramfilename,
        labelfolder,
        labelfilename,
        simulfolder,
        "X_simul_" + str(waveno),
        samplefolder,
    )

    # copy labels in samplefolder
    shutil.copy(libEP_DIR + "/data/xlabels/ToRORd_dynCl.txt", samplefolder)
    shutil.copy(libEP_DIR + "/data/xlabels/xlabels_empty.txt", samplefolder)
    shutil.copy(labelfolder + "/" + labelfilename + ".txt", samplefolder)
    # <--- END


if __name__ == "__main__":
    main()
