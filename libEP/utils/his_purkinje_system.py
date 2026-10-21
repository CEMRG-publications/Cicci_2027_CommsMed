#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file is part of 'Cicci_2027_CommsMed'.

__author__ = Ludovica Cicci
__email__ = l.cicci@imperial.ac.uk
"""

# dependencies
from utils import IO
from utils import meshtool_func as meshtool

# extra imports
import os
import pandas as pd
import numpy as np
import shutil


def find_idx_from_UVC(uvc_target_dct, uvc_all_dct):
    """Find index of closest point on biventricle with given UVC

    Parameters
    ----------
    uvc_target_dct: dict
        Target UVC.
    uvc_all_dct: dct
        Biventricular UVC.
    """

    # Find coord of a point matching the target UVCs.
    idx1 = np.where(abs(uvc_all_dct["Z"] - uvc_target_dct["Z"]) < 0.02)[0]
    idx2 = np.where(abs(uvc_all_dct["RHO"] - uvc_target_dct["RHO"]) < 0.02)[0]
    idx3 = np.where(abs(uvc_all_dct["PHI"] - uvc_target_dct["PHI"]) < 0.02)[0]
    idx4 = np.where(uvc_all_dct["biv"] == uvc_target_dct["biv"])[0]
    idx = np.intersect1d(idx1, idx2)
    idx = np.intersect1d(idx, idx3)
    idx = np.intersect1d(idx, idx4)

    if len(idx) == 0:
        idx1 = np.where(abs(uvc_all_dct["Z"] - uvc_target_dct["Z"]) < 0.05)[0]
        idx2 = np.where(abs(uvc_all_dct["RHO"] - uvc_target_dct["RHO"]) < 0.05)[0]
        idx3 = np.where(abs(uvc_all_dct["PHI"] - uvc_target_dct["PHI"]) < 0.05)[0]
        idx4 = np.where(uvc_all_dct["biv"] == uvc_target_dct["biv"])[0]
        idx = np.intersect1d(idx1, idx2)
        idx = np.intersect1d(idx, idx3)
        idx = np.intersect1d(idx, idx4)

        if len(idx) == 0:
            idx1 = np.where(abs(uvc_all_dct["Z"] - uvc_target_dct["Z"]) < 0.1)[0]
            idx2 = np.where(abs(uvc_all_dct["RHO"] - uvc_target_dct["RHO"]) < 0.1)[0]
            idx3 = np.where(abs(uvc_all_dct["PHI"] - uvc_target_dct["PHI"]) < 0.1)[0]
            idx4 = np.where(uvc_all_dct["biv"] == uvc_target_dct["biv"])[0]
            idx = np.intersect1d(idx1, idx2)
            idx = np.intersect1d(idx, idx3)
            idx = np.intersect1d(idx, idx4)

            if len(idx) == 0:
                idx1 = np.where(abs(uvc_all_dct["Z"] - uvc_target_dct["Z"]) < 0.15)[0]
                idx2 = np.where(abs(uvc_all_dct["RHO"] - uvc_target_dct["RHO"]) < 0.15)[
                    0
                ]
                idx3 = np.where(abs(uvc_all_dct["PHI"] - uvc_target_dct["PHI"]) < 0.15)[
                    0
                ]
                idx4 = np.where(uvc_all_dct["biv"] == uvc_target_dct["biv"])[0]
                idx = np.intersect1d(idx1, idx2)
                idx = np.intersect1d(idx, idx3)
                idx = np.intersect1d(idx, idx4)

    return idx


def define_stimregion(
    pts_submsh, idx_target, idx_SElayer, radius=300, filename="early_activation_site"
):
    """Finds the closest points to a given target point

    Parameters
    ----------
    pts_submsh: dataframe
        Nodes of biventricular (sub)mesh, i.e., either whole myocardium or only subendocardial layer.
    idx_target: int
        Index of root point.
    idx_SElayer: list
        Indices of point with approximately the same transmural depth of the the root point.
    radius: int
        Number of points of the disc-like root. Default: 300.
    filename: str
        Path to output .vtx file with the stimulus indices on biv (sub)mesh. Default: early_activation_site
    """
    # Note: radius is defined by selecting the number of points aiming for
    #       radius = 300 is the optimal value on biv for a sphere diameter equal to 0.05 (z_uvc).

    pts_submsh = np.array(pts_submsh)
    d = ((pts_submsh - pts_submsh[idx_target]) ** 2).sum(axis=1)
    ndx = d.argsort()
    vtx = np.intersect1d(ndx, idx_SElayer)
    vtx = ndx[:radius]
    IO.write_vtx(filename, pd.DataFrame(vtx))

    return


def extract_submesh_EAS(
    CARP_INSTALL_DIR,
    meshfolder,
    meshname,
    finalmeshfolder,
    finalmeshname,
    uvcfolder,
    prestimfolder,
    _submsh,
):
    """Extract biv (sub)mesh from biventricle and map corresponding UVC.

    Parameters
    ----------
    CARP_INSTALL_DIR: str
        Path to CARPentry installation directory.
    meshfolder: str
        Path to the computational mesh.
    meshname: str
        Name of the biventricular computational mesh, eventually with subendocardial layer.
    finalmeshfolder: str
        Path to the directory where the final computational mesh will be saved.
    finalmeshname: str
        Name of the biventricular submesh.
    uvcfolder : str
        Path to folder with universal ventricular coordinates.
    prestimfolder: str
        Path to the directory where the pre-computed quantities for the stimuli are/will be saved.
    _submsh : str
        Pre-exising tag(s) for biventricular submesh.
    """

    # find executables
    meshtoolloc = CARP_INSTALL_DIR + "/meshtool/meshtool"

    # extract submesh (i.e., myo or SElayer) myocardium from biv
    if not meshtool.extract_mesh(
        meshtoolloc,
        meshfolder + "/" + meshname,
        _submsh,
        finalmeshfolder + "/" + finalmeshname,
    ):
        raise Exception("meshtool error. Abort!")

    # initialise dictionary containing UVC for all points in biv (sub)msh
    uvc_submsh_dct = {}

    # extract UVC for submesh (if submesh is biv, returns the same UVC)

    # Z
    if os.path.exists(uvcfolder + "/Z_subendocardium.dat"):
        if not meshtool.extract_data(
            meshtoolloc,
            finalmeshfolder + "/" + finalmeshname,
            uvcfolder + "/Z_subendocardium.dat",  # uvcfolder + "/Z.dat",
            prestimfolder + "/Z_submsh.dat",
            "0",
        ):
            raise Exception("meshtool error. Abort!")
    else:
        if not meshtool.extract_data(
            meshtoolloc,
            finalmeshfolder + "/" + finalmeshname,
            uvcfolder + "/Z.dat",
            prestimfolder + "/Z_submsh.dat",
            "0",
        ):
            raise Exception("meshtool error. Abort!")
    uvc_submsh_dct["Z"] = np.genfromtxt(prestimfolder + "/Z_submsh.dat")

    # RHO
    if not meshtool.extract_data(
        meshtoolloc,
        finalmeshfolder + "/" + finalmeshname,
        uvcfolder + "/RHO.dat",
        prestimfolder + "/RHO_submsh.dat",
        "0",
    ):
        raise Exception("meshtool error. Abort!")
    uvc_submsh_dct["RHO"] = np.genfromtxt(prestimfolder + "/RHO_submsh.dat")

    # PHI
    if not meshtool.extract_data(
        meshtoolloc,
        finalmeshfolder + "/" + finalmeshname,
        uvcfolder + "/PHI.dat",
        prestimfolder + "/PHI_submsh.dat",
        "0",
    ):
        raise Exception("meshtool error. Abort!")
    uvc_submsh_dct["PHI"] = np.genfromtxt(prestimfolder + "/PHI_submsh.dat")

    # biv
    if not meshtool.extract_data(
        meshtoolloc,
        finalmeshfolder + "/" + finalmeshname,
        uvcfolder + "/biv.dat",
        prestimfolder + "/biv_submsh.dat",
        "0",
    ):
        raise Exception("meshtool error. Abort!")
    uvc_submsh_dct["biv"] = np.genfromtxt(prestimfolder + "/biv_submsh.dat")

    return uvc_submsh_dct






##############################
# PARALLELISATION USING POOL #
##############################
def pool_vtx_stimuli_preprocessing(
    CARP_INSTALL_DIR,
    meshfolder,
    meshname,
    uvcfolder,
    samplefolder,
    eas_list=["LV_sf", "LV_pf", "LV_af", "RV_sf", "RV_mod"],
    _RVmyo="34",
    _LVmyo="35",
    SElayer=True,
    _SE_RV="44",
    _SE_LV="45",
    _SE_septRV="46",
):
    """Pre-processing step to compute early activation sites for the biventricle.

    Parameters
    ----------
    CARP_INSTALL_DIR: str
        Path to CARPentry installation directory.
    meshfolder: str
        Path to the computational mesh.
    meshname: str
        Name of the computational mesh, eventually with subendocardial layer.
    uvcfolder : str
        Path to folder with universal ventricular coordinates.
    samplefolder: str
        Path to folder with all sample(s) details.
    eas_list: list
        Early activation site name(s). Default: ["LV_sf", "LV_pf", "LV_af", "RV_sf", "RV_mod"].
    _RVmyo: str
        Pre-existing tag for RV myocardium. Default: "34".
    _LVmyo: str
        Pre-existing tag for LV myocardium. Default: "35".
    SElayer: bool
        If True, use subendocardial layer (from apex to base) to comput the stimulus locations. Default: True.
    _SE_RV: str
        Pre-existing tag for RV subendocardial layer. Default: "44".
    _SE_LV: str
        Pre-existing tag for LV subendocardial layer. Default: "45".
    _SE_septRV: str
        Pre-existing tag for septal subendocardial layer on the right side. Default: "46".
    """

    # find executables
    meshtoolloc = CARP_INSTALL_DIR + "/meshtool/meshtool"

    # IO directories
    stimfolder = samplefolder + "/../STIMULI"
    if not os.path.exists(stimfolder):
        os.mkdir(stimfolder)

    # volumetric meshes
    meshloc = stimfolder + "/volumes"
    if not os.path.exists(meshloc):
        os.mkdir(meshloc)

    # volumetric meshes extracted from biv
    mesh_iloc = meshloc + "/biv"
    if not os.path.exists(mesh_iloc):
        os.mkdir(mesh_iloc)

    # extract biv from torso
    if not meshtool.extract_mesh(
        meshtoolloc,
        meshfolder + "/" + meshname,
        _RVmyo + "," + _LVmyo + "," + _SE_RV + "," + _SE_LV + "," + _SE_septRV,
        meshloc + "/" + meshname + "_i",
    ):
        raise Exception("meshtool error. Abort!")

    # set tag for (sub)mesh on which EAS points are extracted
    _submsh = _SE_RV + "," + _SE_LV + "," + _SE_septRV
    if not SElayer:
        _submsh = _RVmyo + "," + _LVmyo + "," + _submsh

    # EAS pre-processing
    uvc_submsh_dct = extract_submesh_EAS(
        CARP_INSTALL_DIR,
        meshloc,
        meshname + "_i",
        mesh_iloc,
        meshname + "_i_submsh",
        uvcfolder,
        stimfolder,
        _submsh,
    )

    pts_submsh = IO.read_pts(mesh_iloc + "/" + meshname + "_i_submsh")

    # find accepted parameters
    readmefile = samplefolder + "/accepted_samples"
    if not os.path.exists(readmefile + ".txt"):
        readme = [str(i) for i in os.listdir(samplefolder) if not ".txt" in i]
    else:
        readme = IO.read_text(readmefile)
        readme = readme[1:]

    return uvc_submsh_dct, pts_submsh, readme


def pool_vtx_stimuli_sample(
    id_sample,
    CARP_INSTALL_DIR,
    meshfolder,
    meshname,
    uvcfolder,
    samplefolder,
    uvc_submsh_dct,
    pts_submsh,
    eas_list=["LV_sf", "LV_pf", "LV_af", "RV_sf", "RV_mod"],
    SElayer=True,
):
    """Compute early activation sites for the biventricle, for one sample.

    Parameters
    ----------
    id_sample: int
        ID of current sample.
    CARP_INSTALL_DIR: str
        Path to CARPentry installation directory.
    meshfolder: str
        Path to the computational mesh.
    meshname: str
        Name of the computational mesh, eventually with subendocardial layer.
    uvcfolder : str
        Path to folder with universal ventricular coordinates.
    samplefolder: str
        Path to folder with all sample(s) details.
    uvc_submsh_dct: dict
        Biv (sub)mesh from biventricle and corresponding UVC. Output of pool_vtx_stimuli_preprocessing.
    pts_submsh: dict
        Points of biv (sub)mesh. Output of pool_vtx_stimuli_preprocessing.
    eas_list: list
        Early activation site name(s). Default: ["LV_sf", "LV_pf", "LV_af", "RV_sf", "RV_mod"].
    SElayer: bool
        If True, use subendocardial layer (from apex to base) to comput the stimulus locations. Default: True.
    """

    # IO directories
    stimfolder = samplefolder + "/../STIMULI"
    meshloc = stimfolder + "/volumes"
    mesh_iloc = meshloc + "/biv"

    # current simulation paramfile
    paramfile = samplefolder + "/" + str(id_sample) + "/param"
    if not os.path.exists(paramfile + ".json"):
        raise Exception("Cannot find {}.".format(paramfile + ".json"))

    # current sample output directory
    SIMID = samplefolder + "/" + str(id_sample)
    if not os.path.exists(SIMID):
        os.mkdir(SIMID)

    # output directory for EASs
    hpsfolder = SIMID + "/HPS"
    if not os.path.exists(hpsfolder):
        os.mkdir(hpsfolder)

    if len(os.listdir(hpsfolder)) < 5:

        # load dictionary with all simulation parameters
        par_all_dct = IO.read_json(paramfile)

        # dict with early activation sites parameter values
        eas_all_dct = dict()
        for j in eas_list:
            eas_all_dct[j] = {
                "Z": par_all_dct[j + "_Z"],
                "RHO": par_all_dct[j + "_RHO"],
                "PHI": par_all_dct[j + "_PHI"],
                "biv": par_all_dct[j + "_biv"],
            }

        # compute early activation sites for HPS
        vtx_stimuli_mod(
            id_sample,
            CARP_INSTALL_DIR,
            meshloc,
            meshname + "_i",
            mesh_iloc,
            meshname + "_i_submsh",
            stimfolder,
            hpsfolder,
            pts_submsh,
            uvc_submsh_dct,
            eas_all_dct,
            SElayer,
        )
    else:
        print("EAS for sample {} already processed!".format(id_sample))

    return


def vtx_stimuli_mod(
    id_sample,
    CARP_INSTALL_DIR,
    meshfolder,
    meshname,
    finalmeshfolder,
    finalmeshname,
    prestimfolder,
    hpsfolder,
    pts_submsh,
    uvc_submsh_dct,
    eas_all_dct,
    SElayer=True,
):
    """Compute early activation sites for the biventricle.

    Parameters
    ----------
    id_samples: int
        ID of current sample.
    CARP_INSTALL_DIR: str
        Path to CARPentry installation directory.
    meshfolder: str
        Path to the computational mesh.
    meshname: str
        Name of the biventricular computational mesh, eventually with subendocardial layer.
    finalmeshfolder: str
        Path to the computational submesh.
    finalmeshname: str
        Name of the biventricular submesh.
    prestimfolder: str
        Path to the directory where the pre-computed quantities for the stimuli are/will be saved.
    hpsfolder: str
        Path to the directory where the final stimulus locations will be saved.
    pts_submsh: dataframe
        Nodes of biventricular (sub)mesh, i.e., either whole myocardium or only subendocardial layer.
    uvc_submsh_dct: dict
        Biventricular (sub)mesh UVC.
    eas_all_dct: dict
        Early activation site(s) UVC.
    SElayer: bool
        If True, use subendocardial layer (from apex to base) to comput the stimulus locations. Default: True.
    """

    # find executables
    meshtoolloc = CARP_INSTALL_DIR + "/meshtool/meshtool"

    # directory where the simulation results on biv will be stored
    loc = prestimfolder + "/biv_" + str(id_sample)  ##<--DIFFERENCE WRT "VTX_STIMULI"
    if not os.path.exists(loc):
        os.mkdir(loc)

    i_loc = loc + "/submsh"
    if not os.path.exists(i_loc):
        os.mkdir(i_loc)

    if SElayer:
        radius = 200
    else:
        radius = 300

    # define location of activation sites on biventricle (sub)mesh
    for _, eas in enumerate(eas_all_dct):
        uvc_target_dct = eas_all_dct[eas]
        idx_target = find_idx_from_UVC(uvc_target_dct, uvc_submsh_dct)

        if len(idx_target):
            define_stimregion(
                pts_submsh,
                idx_target[0],
                np.where(abs(uvc_submsh_dct["RHO"] - uvc_target_dct["RHO"]) < 0.05)[0],
                radius,
                i_loc + "/" + eas,
            )
            if not meshtool.map_file(
                meshtoolloc,
                finalmeshfolder + "/" + finalmeshname,
                i_loc + "/" + eas + ".vtx",
                loc,
                mode="s2m",
            ):
                raise Exception("meshtool error. Abort!")
            if not meshtool.map_file(
                meshtoolloc,
                meshfolder + "/" + meshname,
                loc + "/" + eas + ".vtx",
                hpsfolder,
                mode="s2m",
            ):
                raise Exception("meshtool error. Abort!")
        else:
            print("Skipping of generation {} for current sample!".format(eas))

    # delete unnecessary files
    shutil.rmtree(loc)  ##<--DIFFERENCE WRT "VTX_STIMULI"

    return
