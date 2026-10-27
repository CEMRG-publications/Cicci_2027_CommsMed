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

    meshname = "KCL1"  # NOTE: mesh basename
    basefolder = OUTPUT_DIR + "/" + meshname

    # target folder (input and outputs are/will be strored here)
    targetfolder = basefolder + "/target"

    # current ECG segment
    ECGsegment = "QRScomplex"

    # PCA threshold
    threshold_PCA = 0.99

    # ---> START
    qoifolder = targetfolder + "/qoi_" + ECGsegment

    thPCA = str(int(100 * threshold_PCA))

    if "QRScomplex" in ECGsegment:
        PCAbasisfolder = basefolder + "/PCAbasis_" + meshname + "_150ms_" + thPCA + "%"
        qoi_list = ["12leads_truncated"]
    elif "Twave" in ECGsegment:
        PCAbasisfolder = basefolder + "/PCAbasis_" + meshname + "_REP_" + thPCA + "%"
        qoi_list = ["12leads_REP"]
    else:
        PCAbasisfolder = basefolder + "/PCAbasis_" + meshname + "_" + thPCA + "%"
        qoi_list = ["12leads"]

    gsa.extract_output_from_clinical_ecg(
        PCAbasisfolder,
        targetfolder,
        "ECG_from_BSP", # clinical 12-lead ECG .dat file (columns: time, LI, LII, LIII, aVR, aVL, aVF, V1-V6)
        qoifolder,
        qoi_list=qoi_list,
        tend_truncated=150,
    )
    # <--- END


if __name__ == "__main__":
    main()
