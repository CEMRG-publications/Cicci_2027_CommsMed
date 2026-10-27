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
from utils import setup


def main():
    check_paths("CARP_INSTALL_DIR", "OUTPUT_DIR")

    basefolder = OUTPUT_DIR + "/KCL1"

    # name of mesh to be used in EP simulation(s)
    meshname = "KCL1_gregion"

    # name of ionic model and name of labelfile
    ionicmodelname = "ToRORd_dynCl"
    labelfilename_im = "ToRORd_dynCl"

    # name of label .txt file with scalar input parameters
    labelfilename_scalar = "xlabels_scalar"

    # directory with ECG electrode coordinates
    electfolder = OUTPUT_DIR + "/KCL1_electrodes"

    # file with ECG electrode coordinates
    electfilename = "electrodes"

    # directory where samples are saved
    samplefolder = basefolder + "/samples_LHS"

    # directory where "mesh" and "data" are temporarily saved by this script.
    # Move or copy them into clusterfolder before running the simulations (see workflow.md)
    storagefolder = basefolder + "/basics"

    # directory where simulations will be run: absolute path to "mesh", "data" and samplefolder
    # This path is written in the parameter files and read by CARP.
    clusterfolder = basefolder

    # set propagation type
    propagation = "RE+"  # "RE-" or "RE+"

    # if True, use .vtx file(s) for EAS common among samples and stored in /path/to/ref_QRS/0/HPS
    hps_ref = False

    # folder with CVs .adj files, if common among samples (if not, use empty string "")
    cvfolder_ref = "" #basefolder + "/ref_QRS/0/CV"

    # folder with fibre .lon file, if different than default (if not, use empty string "")
    lonfolder_ref = "" #basefolder + "/ref_QRS/0"

    # ---> START
    # final mesh folder
    storage_meshfolder = storagefolder + "/mesh"
    os.makedirs(storage_meshfolder, exist_ok=True)

    # final data folder
    storage_datafolder = storagefolder + "/data"
    os.makedirs(storage_datafolder, exist_ok=True)

    # copy final mesh into storage folder
    setup.copy_mesh(
        CARP_INSTALL_DIR,
        basefolder,
        meshname,
        storage_meshfolder,
        meshname,
        extract_myo=False,
    )

    # copy datafile into storage folder
    setup.copy_par_file(libEP_DIR + "/data/par_files", storage_datafolder)

    # copy electrode file into storage folder
    setup.copy_electrode_file(electfolder, electfilename, storage_datafolder)

    # find electrode vertices on current torso mesh
    elect_vtx = setup.electrodes_pts2vtx(
        basefolder, meshname, electfolder, electfilename, storage_datafolder
    )

    # write lead_field basic commands
    setup.update_par_lead_field(
        storage_datafolder, "pfile_LF", elect_vtx[-2], clusterfolder + "/data"
    )

    # write ECG basic commands
    setup.update_par_ECG(
        storage_datafolder, "pfile_ECG", elect_vtx[-2], clusterfolder + "/data"
    )

    # write EP simulation basic commands
    setup.write_par_EPmodel(
        samplefolder,
        clusterfolder + "/" + os.path.basename(samplefolder),
        clusterfolder + "/data",
        propagation,
        "monodomain",
        hps_ref=hps_ref
    )

    # write sample-specific tissue paramfiles
    setup.write_par_tissue(
        CARP_INSTALL_DIR,
        basefolder,
        meshname,
        labelfilename_im,
        samplefolder,
        labelfilename_scalar,
        samplefolder,
        clusterfolder + "/" + os.path.basename(samplefolder),
        libEP_DIR + "/data/tag_BATH",
        libEP_DIR + "/data/tune_CV_ToRORd_dynCl",
        im=ionicmodelname,
        prepacing_lats="", # path to LATs.dat file
        prepacing_beats=100,
        prepacing_bcl=1000,
        lonfolder_ref=lonfolder_ref,
        cvfolder_ref=cvfolder_ref,
        tag_HEART_jsonfile=libEP_DIR + "/data/tag_HEART",
    )
    # <--- END

    return


if __name__ == "__main__":
    main()
