#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file is part of 'Cicci_2027_CommsMed'.

__author__ = Ludovica Cicci
__email__ = l.cicci@imperial.ac.uk
"""

# other imports
import os
import sys
from multiprocessing import Pool
from pytictoc import TicToc

# make config.py (in the parent folder) importable
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# dependencies
from config import CARP_INSTALL_DIR, libEP_DIR, OUTPUT_DIR, check_paths
from utils import pool_func as pf
from utils import his_purkinje_system as hps
from utils import heterogeneity as htr
from utils import carp_func
from utils import processing_mesh as proc
from utils import setup


if __name__ == "__main__":
    check_paths("CARP_INSTALL_DIR", "OUTPUT_DIR")

    basefolder = OUTPUT_DIR + "/KCL1"
    meshname = "KCL1"  # NOTE: mesh basename

    # name of .txt file(s) with input parameters
    labelfilename_scalar = "xlabels_scalar"
 
    # directory where samples are saved
    samplefolder = basefolder + "/samples_LHS"

    # name of .txt file with ionic model parameters for heterogeneity adjustment
    labelfilename_im = "ToRORd_dynCl" # if any, otherwise set to ""

    # If true, look for EAS in subendolayer mesh
    SElayer = True

    # If True, write heterogeneity adjustment-like file(s) for ionic model parameters and conduction velocities
    heterogeneity = True

    # number of processors
    num_proc = 8

    # ---> START

    time_total = TicToc()
    time_total.tic()

    # read mesh tags from json file
    tags, _ = proc.get_heart_tags(libEP_DIR + "/data/tag_HEART", libEP_DIR + "/data/tag_BATH")

    ## ------------------ find vertices of EAS roots ------------------ ##
    # EAS
    eas_list = ["LV_sf", "LV_pf", "LV_af", "RV_sf", "RV_mod"]

    # name of mesh currently processed
    if SElayer:
        current_mesh_name = meshname + "_subendolayer"
    else:
        current_mesh_name = meshname

    # pre-processing for EAS
    uvc_submsh_dct, pts_submsh, readme = hps.pool_vtx_stimuli_preprocessing(
        CARP_INSTALL_DIR,
        basefolder,
        current_mesh_name,
        basefolder + "/UVC",
        samplefolder,
        eas_list=eas_list,
        _RVmyo=tags["RVmyo_subendolayer"],
        _LVmyo=tags["LVmyo_subendolayer"],
        SElayer=SElayer,
        _SE_RV=tags["SE_RV_subendolayer"],
        _SE_LV=tags["SE_LV_subendolayer"],
        _SE_septRV=tags["SE_septRV_subendolayer"],
    )

    # define argument list for process
    argv = [
        (
            i,
            basefolder,
            current_mesh_name,
            basefolder + "/UVC",
            samplefolder,
            uvc_submsh_dct,
            pts_submsh,
            eas_list,
            SElayer,
        )
        for i in readme
    ]

    # start parallel process
    time = TicToc()
    time.tic()
    with Pool(processes=num_proc) as pool:
        pool.starmap(pf.process_vtx_stimuli, argv)
    time.toc("Compute vtx stimuli")

    ##  ------------------ HETEROGENEITY  ------------------ ##
    if heterogeneity:

        if labelfilename_im:

            ##  ------------------ compute adj-files for ion conductances  ------------------ ##
            # ADJ
            domain = "intra"  # either "intra" or "extra"

            # pre-processing for ADJ
            (
                compute_adj,
                RHO,
                Z,
                par_keys,
                ion_var_list,
                vtx_dct,
                readme,
            ) = htr.pool_adjustment_file_preprocessing(
                basefolder + "/UVC",
                samplefolder,
                samplefolder,
                labelfilename_im,
            )

            if compute_adj:
                # define argument list for process
                argv = [
                    (
                        i,
                        samplefolder,
                        vtx_dct,
                        RHO,
                        Z,
                        par_keys,
                        ion_var_list,
                        domain,
                    )
                    for i in readme
                ]

                # start parallel process
                time = TicToc()
                time.tic()
                with Pool(processes=num_proc) as pool:
                    pool.starmap(pf.process_adj_file, argv)
                time.toc("Compute adjustment files")

        ## ------------------ compute adj-files for conduction velocities ------------------ ##

        # name of mesh currently processed
        current_mesh_name = meshname + "_gregion"

        # pre-processing for CV adj
        (
            compute_CV,
            meshpath_i,
            vtx_dct,
            meanRHO_dct,
            meanZ_dct,
            meanRHO_SE_dct,
            meanZ_SE_dct,
            par_keys,
            readme,
        ) = htr.pool_CV_adjustment_file_preprocessing(
            CARP_INSTALL_DIR,
            basefolder,
            current_mesh_name,
            basefolder + "/UVC",
            samplefolder,
            _RVmyo=tags["RVmyo_gregion"],
            _LVmyo=tags["LVmyo_gregion"],
            _SE_RV=tags["SE_RV_gregion"],
            _SE_LV=tags["SE_LV_gregion"],
            _SE_septRV=tags["SE_septRV_gregion"],
        )

        if compute_CV:
            # define argument list for process
            argv = [
                (
                    i,
                    samplefolder,
                    meshpath_i,
                    vtx_dct,
                    meanRHO_dct,
                    meanZ_dct,
                    meanRHO_SE_dct,
                    meanZ_SE_dct,
                    par_keys,
                )
                for i in readme
            ]

            # start parallel process
            time = TicToc()
            time.tic()
            with Pool(processes=num_proc) as pool:
                pool.starmap(pf.process_cv_adj_file, argv)
            time.toc("Compute CV adj-like files")

    ## ------------------ compute fibres  ------------------ ##

    # name of mesh currently processed
    current_mesh_name = meshname  # original torso mesh, without subendocardial layer

    # pre-processing for fibres
    compute_fibres, readme = carp_func.pool_compute_fibre_preprocessing(
        CARP_INSTALL_DIR,
        OUTPUT_DIR + "/mesh",
        current_mesh_name,  # original torso mesh, without subendocardial layer
        samplefolder,
        samplefolder,
        labelfilename_scalar,
        _RVmyo=tags["RVmyo"],
        _LVmyo=tags["LVmyo"],
        _RVbp=tags["RVbp"],
        _LVbp=tags["LVbp"],
        _upper_RH=tags["upper_RH"],
        _upper_LH=tags["upper_LH"],
    )

    if compute_fibres:
        # define argument list for process
        argv = [
            (
                i,
                CARP_INSTALL_DIR,
                OUTPUT_DIR + "/mesh",
                current_mesh_name,  # original torso mesh, without subendocardial layer
                samplefolder,
                meshname + "_fibres",
            )
            for i in readme
        ]

        # start parallel process
        time = TicToc()
        time.tic()
        with Pool(processes=num_proc) as pool:
            pool.starmap(pf.process_fibre, argv)
        time.toc("Compute fibres")

    ## check sample folder(s)
    setup.check_sample_folder(samplefolder)

    time_total.toc()
    # <--- END
