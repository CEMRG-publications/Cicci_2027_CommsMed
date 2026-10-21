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
from config import CARP_INSTALL_DIR, libEP_DIR, OUTPUT_DIR, check_paths
from utils import processing_mesh as proc
from utils import carp_func


def main():
    check_paths("CARP_INSTALL_DIR", "OUTPUT_DIR")

    basefolder = OUTPUT_DIR + "/KCL"
    meshname = "KCL1"

    # directory with current mesh in carp_txt format
    meshfolder = OUTPUT_DIR + "/torsomesh"

    # ---> START
    os.makedirs(basefolder, exist_ok=True)

    ### MESH TAGS TO BE USED (mesh-dependent) ###
    _RVmyo="34"  # existing tag for the RV myocardium
    _LVmyo="35"  # existing tag for the LV myocardium
    _SE_RV="44"  # new tag for the RV sub-endocardial layer
    _SE_LV="45"  # new tag for the LV sub-endocardial layer
    _SE_septRV="46"  # new tag for the RV-septal sub-endocardial layer
    _RVbp="24"  # existing tag for the RV blood pool
    _LVbp="25"  # existing tag for the LV blood pool
    _upper_RH="11,12,17,19,22,23"  # existing tag for the right heart (excluding RV myocardium and blood pool)
    _upper_LH="13,14,16,18,20,21"  # existing tag for the left heart (excluding LV myocardium and blood pool)
    uvctagfile=libEP_DIR + "/data/tag_CARP_UVC.sh"
    ##############################################

    ## --------------- UVC --------------- ##
    # Compute UVC using CARPentry code
    carp_func.compute_uvc(
        CARP_INSTALL_DIR,
        meshfolder,
        meshname,
        uvctagfile,
        basefolder + "/UVC",
        _RVmyo=_RVmyo,
        _LVmyo=_LVmyo,
        _RVbp=_RVbp,
        _LVbp=_LVbp,
        _upper_RH=_upper_RH,
        _upper_LH=_upper_LH,
        cleandir=True,
    )

    ## --------------- SUBENDOCARDIAL LAYER --------------- ##
    proc.extract_SElayer(
        CARP_INSTALL_DIR,
        meshfolder,
        meshname,
        basefolder,
        meshname + "_subendolayer",
        meshname + "_SElayer",
        basefolder + "/UVC",
        _RVmyo=_RVmyo,
        _LVmyo=_LVmyo,
        _RVbp=_RVbp,
        _LVbp=_LVbp,
        _upper_RH=_upper_RH,
        _upper_LH=_upper_LH,
        _SE_RV=_SE_RV,
        _SE_LV=_SE_LV,
        _SE_septRV=_SE_septRV,
        cleandir=True,
    )

    ## --------------- FAST ENDOCARDIAL CONDUCTION LAYER(S) --------------- ##
    proc.compute_FEC_layers(
        CARP_INSTALL_DIR,
        basefolder,
        meshname + "_subendolayer",
        basefolder,
        meshname + "_FEC",
        basefolder + "/UVC",
        _RVmyo=_RVmyo,
        _LVmyo=_LVmyo,
        _SE_RV=_SE_RV,
        _SE_LV=_SE_LV,
        _SE_septRV=_SE_septRV,
        cleandir=True,
    )

    ## --------------- MYOCARDIUM SEGMENTATION --------------- ##
    proc.segment_biventricle(
        CARP_INSTALL_DIR,
        basefolder,
        meshname + "_FEC",
        basefolder,
        meshname + "_gregion",
        basefolder + "/UVC",
        _RVmyo=_RVmyo,
        _LVmyo=_LVmyo,
        _SE_RV=_SE_RV + "," + ",".join([str(i) for i in range(int(_SE_RV)*1000+50, int(_SE_RV)*1000+75, 2)]),
        _SE_LV=_SE_LV + "," + ",".join([str(i) for i in range(int(_SE_LV)*1000+50, int(_SE_LV)*1000+75, 2)]),
        _SE_septRV=_SE_septRV + "," + ",".join([str(i) for i in range(int(_SE_septRV)*1000+50, int(_SE_septRV)*1000+75, 2)]),
        Nz=10,
        Nrho=5,
        cleandir=True,
    )
    # <--- END


if __name__ == "__main__":
    main()
