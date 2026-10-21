#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file is part of 'Cicci_2027_CommsMed'.

__author__ = Ludovica Cicci
__email__ = l.cicci@imperial.ac.uk
"""

# dependencies
from config import CARP_INSTALL_DIR
from utils import his_purkinje_system as hps
from utils import heterogeneity as htr
from utils import carp_func


def process_vtx_stimuli(
    id_sample,
    meshfolder,
    meshname,
    uvcfolder,
    samplefolder,
    uvc_submsh_dct,
    pts_submsh,
    eas_list,
    SElayer,
):
    hps.pool_vtx_stimuli_sample(
        id_sample,
        CARP_INSTALL_DIR,
        meshfolder,
        meshname,
        uvcfolder,
        samplefolder,
        uvc_submsh_dct,
        pts_submsh,
        eas_list=eas_list,
        SElayer=SElayer,
    )

    return


def process_adj_file(
    id_sample,
    samplefolder,
    vtx_dct,
    RHO,
    Z,
    par_keys,
    ion_var_list,
    domain,
):

    htr.pool_adjustment_file_sample(
        id_sample,
        samplefolder,
        vtx_dct,
        RHO,
        Z,
        par_keys,
        ion_var_list,
        domain,
    )

    return


def process_cv_adj_file(
    id_sample,
    samplefolder,
    meshpath_i,
    vtx_dct,
    meanRHO_dct,
    meanZ_dct,
    meanRHO_SE_dct,
    meanZ_SE_dct,
    par_keys,
):

    htr.pool_CV_adjustment_file_sample(
        id_sample,
        samplefolder,
        meshpath_i,
        vtx_dct,
        meanRHO_dct,
        meanZ_dct,
        meanRHO_SE_dct,
        meanZ_SE_dct,
        par_keys,
    )

    return


def process_fibre(
    id_sample,
    CARP_INSTALL_DIR,
    meshfolder,
    meshname,  # original torso mesh, without subendocardial layer
    samplefolder,
    finalmeshname,
):

    carp_func.pool_compute_fibre_sample(
        id_sample,
        CARP_INSTALL_DIR,
        meshfolder,
        meshname,  # original torso mesh, without subendocardial layer
        samplefolder,
        finalmeshname,
    )

    return
