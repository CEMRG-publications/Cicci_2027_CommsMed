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
from config import OUTPUT_DIR, check_paths
from utils import global_sensitivity_analysis as gsa


def main():
    check_paths("OUTPUT_DIR")

    basefolder = OUTPUT_DIR + "/KCL1"
    meshname = "KCL1"

    # current ECG segment
    ECGsegment = "QRScomplex"

    # directory where all qoi will be saved (usually, local HD)
    qoifolder = basefolder + "/target/qoi_" + ECGsegment

    # PCA threshold
    threshold_PCA = 0.99

    # ---> START
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

    # directory and name of the clinical 12-lead ECG .dat file (columns: time, LI, LII, LIII, aVR, aVL, aVF, V1-V6)
    ecgdatafolder = basefolder + "/target"
    ecgdatafilename = "ECG_from_BSP"

    gsa.extract_output_from_clinical_ecg(
        PCAbasisfolder,
        ecgdatafolder,
        ecgdatafilename,
        qoifolder,
        qoi_list=qoi_list,
        tend_truncated=150,
    )
    # <--- END


if __name__ == "__main__":
    main()
