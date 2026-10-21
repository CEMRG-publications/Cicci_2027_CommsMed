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
from config import CARP_INSTALL_DIR, OUTPUT_DIR, check_paths
from utils import pool_func as pf
from utils import his_purkinje_system as hps
from utils import heterogeneity as htr
from utils import carp_func
from utils import setup


if __name__ == "__main__":
    check_paths("CARP_INSTALL_DIR", "OUTPUT_DIR")

    basefolder = OUTPUT_DIR + "/KCL"
    meshname = "KCL1"  # NOTE: mesh basename

    # name of .txt file(s) with input parameters
    labelfilename_scalar = "xlabels_scalar"

    # name of .txt file(s) with input parameters for heterogeneity adjustment (ionic model parameters)
    labelfilename_im = "ToRORd_dynCl" # if any, otherwise set to ""
 
    # directory where samples are saved
    samplefolder = basefolder + "/samples"

    # If True, write heterogeneity adjustment-like file(s)
    heterogeneity = True

    # number of processors
    num_proc = 8

    # ---> START

    time_total = TicToc()
    time_total.tic()

    ## ------------------ find vertices of EAS roots ------------------ ##

    ### MESH TAGS TO BE USED (mesh-dependent) ###
    _RVmyo="34" + "," + ",".join([str(i) for i in range(3401, 3431)])  # existing tag for the RV myocardium
    _LVmyo="35" + "," + ",".join([str(i) for i in range(3501, 3551)])  # existing tag for the LV myocardium
    _SE_RV="44" + "," + ",".join([str(i) for i in range(4401, 4411)]) + "," + ",".join([str(i) for i in range(445001, 447511)])  # new tag for the RV sub-endocardial layer
    _SE_LV="45" + "," + ",".join([str(i) for i in range(4501, 4511)]) + "," + ",".join([str(i) for i in range(455001, 457511)])  # new tag for the LV sub-endocardial layer
    _SE_septRV="46" + "," + ",".join([str(i) for i in range(4601, 4611)]) + "," + ",".join([str(i) for i in range(465001, 467511)])  # new tag for the RV-septal sub-endocardial layer
    _RVbp="24"  # existing tag for the RV blood pool
    _LVbp="25"  # existing tag for the LV blood pool
    _upper_RH="11,12,17,19,22,23"  # existing tag for the right heart (excluding RV myocardium and blood pool)
    _upper_LH="13,14,16,18,20,21"  # existing tag for the left heart (excluding LV myocardium and blood pool)
    ##############################################

    # EAS
    eas_list = ["LV_sf", "LV_pf", "LV_af", "RV_sf", "RV_mod"]
    SElayer = True  # If true, look for EAS in subendolayer mesh

    # pre-processing for EAS
    uvc_submsh_dct, pts_submsh, readme = hps.pool_vtx_stimuli_preprocessing(
        CARP_INSTALL_DIR,
        basefolder,
        meshname + "_subendolayer",
        basefolder + "/UVC",
        samplefolder,
        eas_list=eas_list,
        _RVmyo=_RVmyo,
        _LVmyo=_LVmyo,
        SElayer=SElayer,
        _SE_RV=_SE_RV,
        _SE_LV=_SE_LV,
        _SE_septRV=_SE_septRV,
    )

    # define argument list for process
    argv = [
        (
            i,
            basefolder,
            meshname + "_subendolayer",
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
            meshname + "_gregion",
            basefolder + "/UVC",
            samplefolder,
            _RVmyo=_RVmyo,
            _LVmyo=_LVmyo,
            _SE_RV=_SE_RV,
            _SE_LV=_SE_LV,
            _SE_septRV=_SE_septRV,
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
    # pre-processing for fibres
    compute_fibres, readme = carp_func.pool_compute_fibre_preprocessing(
        CARP_INSTALL_DIR,
        OUTPUT_DIR + "/torsomesh",
        meshname,  # original torso mesh, without subendocardial layer
        samplefolder,
        samplefolder,
        labelfilename_scalar,
        _RVmyo=_RVmyo,
        _LVmyo=_LVmyo,
        _RVbp=_RVbp,
        _LVbp=_LVbp,
        _upper_RH=_upper_RH,
        _upper_LH=_upper_LH,
    )

    if compute_fibres:
        # define argument list for process
        argv = [
            (
                i,
                CARP_INSTALL_DIR,
                OUTPUT_DIR + "/torsomesh",
                meshname,  # original torso mesh, without subendocardial layer
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
