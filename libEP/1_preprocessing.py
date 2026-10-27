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

    meshname = "KCL1"
    basefolder = OUTPUT_DIR + "/" + meshname

    # directory with current mesh in carp_txt format
    meshfolder = OUTPUT_DIR + "/mesh"

    # ---> START

    # IO
    os.makedirs(basefolder, exist_ok=True)

    # read mesh tags from json file
    tags, par = proc.get_heart_tags(libEP_DIR + "/data/tag_HEART", libEP_DIR + "/data/tag_BATH")

    ## --------------- UVC --------------- ##
    # Compute UVC using CARPentry code
    carp_func.compute_uvc(
        CARP_INSTALL_DIR,
        meshfolder,
        meshname,
        libEP_DIR + "/data/tag_CARP_UVC.sh",
        basefolder + "/UVC",
        _RVmyo=tags["RVmyo"],
        _LVmyo=tags["LVmyo"],
        _RVbp=tags["RVbp"],
        _LVbp=tags["LVbp"],
        _upper_RH=tags["upper_RH"],
        _upper_LH=tags["upper_LH"],
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
        _RVmyo=tags["RVmyo"],
        _LVmyo=tags["LVmyo"],
        _RVbp=tags["RVbp"],
        _LVbp=tags["LVbp"],
        _upper_RH=tags["upper_RH"],
        _upper_LH=tags["upper_LH"],
        _SE_RV=tags["SE_RV"],
        _SE_LV=tags["SE_LV"],
        _SE_septRV=tags["SE_septRV"],
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
        _RVmyo=tags["RVmyo"],
        _LVmyo=tags["LVmyo"],
        _SE_RV=tags["SE_RV"],
        _SE_LV=tags["SE_LV"],
        _SE_septRV=tags["SE_septRV"],
        factor_tag_FEC=par["factor_tag_FEC"],
        Z_FEC_min=par["Z_FEC_min"],
        Z_FEC_max=par["Z_FEC_max"],
        deltaZ_FEC=par["deltaZ_FEC"],
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
        _RVmyo=tags["RVmyo"],
        _LVmyo=tags["LVmyo"],
        # sub-endocardial layer tag, plus the tags of the fast endocardial conductive layers
        _SE_RV=tags["SE_RV_FEC"],
        _SE_LV=tags["SE_LV_FEC"],
        _SE_septRV=tags["SE_septRV_FEC"],
        factor_tag_myo=par["factor_tag_myo"],
        factor_tag_SE=par["factor_tag_SE"],
        Nz=par["Nz"],
        Nrho_LV=par["Nrho_LV"],
        Nrho_RV=par["Nrho_RV"],
        cleandir=True,
    )

    # number of regions of the segmented mesh (CARPentry allows at most 256)
    proc.count_heart_regions(basefolder, meshname + "_gregion", tags)
    # <--- END


if __name__ == "__main__":
    main()
