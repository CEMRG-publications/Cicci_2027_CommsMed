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
import numpy as np


# ----------------------- IMP REGION ----------------------- #
def find_biv_vtx(uvcfolder):
    """Find vertices of biventricle in biventricular mesh.

    Parameters
    ----------
    uvcfolder : str
        Path to folder with universal ventricular coordinates.
    """

    # load UVC
    biv = IO.read_data(uvcfolder + "/biv")

    # dict of node indices
    vtx_dct = {}
    vtx_dct["biv"] = [i for i in range(len(biv))]

    return vtx_dct


def adjustment_file(
    adjfolder,
    adjfilename,  # same of ionic variable
    vtx_dct,
    RHO,
    Z,
    theta_dct,
    grad_rho_dct,
    grad_z_dct,
    domain="intra",
    visualise=False,
):
    """Generate adjustment file to impose heterogeneity in cell-scale dynamics (imp_region).

    Parameters
    ----------
    adjfolder: str
        Path to the directory where the vertex adjustment file will be saved.
    adjfilename: str
        Name of the vertex adjustment file. Must be the same name of the variable being adjusted.
    vtx_dct: dict
        List of vertices for each biventricular sub-region.
    RHO: list
        Rotational UVC coordinate.
    Z: list
        Apico-basal UVC coordinate.
    theta_dct: dct
        Value of theta (i.e., apex_epi) for current variable, in each sub-region.
    grad_rho_dct: dct
        Value of grad_rho (i.e., apex_endo / apex_epi) for current variable, in each sub-region.
    grad_z_dct: dct
        Value of grad_z (i.e., base_epi / apex_epi) for current variable, in each sub-region.
    domain: str
        Definition domain (either "intra" or "extra"). Default: "intra".
    visualise: bool
        If True, save vertex values in .dat datafile to be visualised in meshalyzer. Default: False.
    """

    # Define
    #   1- apex_endo/apex_epi = beta_endo/beta_epi =: grad_rho
    #   2- beta_epi/apex_epi = beta_endo/apex_endo =: grad_z
    # such that
    #   1- apex_epi =: theta
    #   2- apex_endo = theta * grad_rho
    #   3- beta_epi = theta * grad_z
    #   4- beta_endo = theta * grad_rho * grad_z

    # total number of vertices
    nvtx = len(RHO)

    # initialisation of adjustment files
    adjfile = adjfolder + "/" + adjfilename
    adjf = open(adjfile + ".adj", "w")

    # add header to .adj file
    adjf.write(str(nvtx) + "\n" + domain + "\n")

    # define adjusted ionic conductance wrt to each heterogeneity sub-region
    ion_var = np.zeros((nvtx,))
    for _, subreg in enumerate(vtx_dct):
        idx_subreg = list(vtx_dct[subreg])
        ion_var[idx_subreg] = (
            theta_dct[subreg] * np.multiply(1 - Z[idx_subreg], RHO[idx_subreg])
            + theta_dct[subreg]
            * grad_rho_dct[subreg]
            * np.multiply(1 - Z[idx_subreg], 1 - RHO[idx_subreg])
            + theta_dct[subreg]
            * grad_z_dct[subreg]
            * np.multiply(Z[idx_subreg], RHO[idx_subreg])
            + theta_dct[subreg]
            * grad_rho_dct[subreg]
            * grad_z_dct[subreg]
            * np.multiply(Z[idx_subreg], 1 - RHO[idx_subreg])
        )

    # output .adj file
    for i in range(nvtx):
        adjf.write("{} {}\n".format(i, ion_var[i]))
    adjf.close()

    # visualisation
    if visualise:
        datf = open(adjfile + ".dat", "w")
        for i in range(nvtx):
            datf.write("{}\n".format(ion_var[i]))
        datf.close()

    return




# ----------------------- GREGION ----------------------- #
def find_gregion_vtx(
    elem_biv, _RVmyo, _LVmyo, _SE_RV, _SE_LV, _SE_septRV
):
    """Find vertices of each gregion in biventricular mesh.

    Parameters
    ----------
    elem_biv: Dataframe
        Elements of biventricular mesh, including FEC and regional tags (if any).
    _RVmyo: str
        Pre-existing tag(s) for RV myocardium.
    _LVmyo: str
        Pre-existing tag(s) for LV myocardium.
    _SE_RV: str
        Pre-existing tag(s) for RV subendocardial layer.
    _SE_LV: str
        Pre-existing tag(s) for LV subendocardial layer.
    _SE_septRV: str
        Pre-existing tag(s) for septal subendocardial layer on the right side.
    """

    # extract tag(s) (int) for each region
    tag_dct = {}
    tag_dct["biv"] = [
        int(i)
        for i in (
            _RVmyo + "," + _LVmyo + "," + _SE_RV + "," + _SE_LV + "," + _SE_septRV
        ).split(",")
    ]
    tag_dct["SE_biv"] = [
        int(i) for i in (_SE_RV + "," + _SE_LV + "," + _SE_septRV).split(",")
    ]

    # find vertices for each tag in each subregion
    print("Getting vertices for conductivity regions...")
    vtx_dct_tmp = {}
    for _, subreg in enumerate(tag_dct):
        vtx_dct_tmp[subreg] = {}
        for _, subreg_tag in enumerate(tag_dct[subreg]):
            vtx_subreg_tag = elem_biv.loc[
                np.where(elem_biv[5] == subreg_tag)[0], [1, 2, 3, 4]
            ].to_numpy()
            vtx_dct_tmp[subreg][subreg_tag] = np.unique(vtx_subreg_tag).tolist()
        print("...{}. Done.".format(subreg))

    # get rid of empty fields
    vtx_dct = {}
    for _, subreg in enumerate(vtx_dct_tmp):
        vtx_dct[subreg] = {}
        for _, subreg_tag in enumerate(vtx_dct_tmp[subreg]):
            if vtx_dct_tmp[subreg][subreg_tag]:
                vtx_dct[subreg][subreg_tag] = vtx_dct_tmp[subreg][subreg_tag]

    return vtx_dct


def update_FEC_vtx(vtx_dct, Zmax_SE):
    """Update dictionary of gregion vertices based on FEC extension.

    Parameters
    ----------
    vtx_dct: dict
        Dictionary of vertices for each gregion.
    Zmax_SE: float
        Extension in the apico-basal direction of the fast subendocardial layer. Must be in (0,1].
    """

    vtx_red_dct = {}
    for _, subreg in enumerate(vtx_dct):
        if "SE_" in subreg:
            vtx_red_dct[subreg] = {}
            for _, subreg_tag in enumerate(vtx_dct[subreg]):
                if (
                    len(str(subreg_tag)) > 2
                    and int(str(subreg_tag)[2:4]) <= 1e2 * Zmax_SE
                ):
                    vtx_red_dct[subreg][subreg_tag] = vtx_dct[subreg][subreg_tag]
        else:
            vtx_red_dct[subreg] = vtx_dct[subreg]

    return vtx_red_dct


def adjustment_gregion(
    adjfolder,
    adjfilename,  # same of CV variable
    vtx_dct,
    RHO_dct,
    Z_dct,
    theta_dct,
    grad_rho_dct,
    grad_z_dct,
    meshpath_i,
    visualise=False,
):
    """Generate datafiles to impose heterogeneity in conductivity region (gregion).

    Parameters
    ----------
    adjfolder: str
        Path to the directory where the vertex adjustment-like file will be saved.
    adjfilename: str
        Name of the vertex adjustment-like file. Must be the same name of the variable being adjusted.
    vtx_dct: dict
        Dictionary of vertices for each gregion.
    RHO_dct: dict
        Rotational UVC coordinate for each gregion.
    Z_dct: dict
        Apico-basal UVC coordinate for each gregion.
    theta_dct: dct
        Value of theta (i.e., apex_epi) for current variable, in each heterogeneity sub-region.
    grad_rho_dct: dct
        Value of grad_rho (i.e., apex_endo / apex_epi) for current variable, in each heterogeneity sub-region.
    grad_z_dct: dct
        Value of grad_z (i.e., base_epi / apex_epi) for current variable, in each heterogeneity sub-region.
    meshpath_i: str
        Path to biventricular mesh, including FEC and regional tags (if any).
    visualise: bool
        If True, save regional element values in .elem file to be visualised in paraview. Default: False.
    """

    # Define
    #   1- apex_endo/apex_epi = beta_endo/beta_epi =: grad_rho
    #   2- beta_epi/apex_epi = beta_endo/apex_endo =: grad_z
    # such that
    #   1- apex_epi =: theta
    #   2- apex_endo = theta * grad_rho
    #   3- beta_epi = theta * grad_z
    #   4- beta_endo = theta * grad_rho * grad_z

    # initialisation of adjustment files
    adjfile = adjfolder + "/" + adjfilename
    adjf = open(adjfile + ".adj", "w")

    # add header to .adj file
    adjf.write("ID val\n")

    if visualise:
        elem_biv = IO.read_elem(meshpath_i)
        elem_tag = elem_biv[5]
        elem_biv[5] = 0

    # define adjusted conduction velocity wrt to each heterogeneity sub-region
    for _, subreg in enumerate(vtx_dct):
        for _, subreg_tag in enumerate(vtx_dct[subreg]):
            cv_var = (
                theta_dct[subreg]
                * np.multiply(
                    1 - Z_dct[subreg][subreg_tag], RHO_dct[subreg][subreg_tag]
                )
                + theta_dct[subreg]
                * grad_rho_dct[subreg]
                * np.multiply(
                    1 - Z_dct[subreg][subreg_tag], 1 - RHO_dct[subreg][subreg_tag]
                )
                + theta_dct[subreg]
                * grad_z_dct[subreg]
                * np.multiply(Z_dct[subreg][subreg_tag], RHO_dct[subreg][subreg_tag])
                + theta_dct[subreg]
                * grad_rho_dct[subreg]
                * grad_z_dct[subreg]
                * np.multiply(
                    Z_dct[subreg][subreg_tag], 1 - RHO_dct[subreg][subreg_tag]
                )
            )
            adjf.write("{} {}\n".format(subreg_tag, round(cv_var, 2)))

            if visualise:
                elem_biv.loc[np.where(elem_tag == subreg_tag)[0], 5] = int(1e3 * cv_var)

    adjf.close()

    if visualise:
        IO.write_elem(adjfile, elem_biv)

    return




##############################
# PARALLELISATION USING POOL #
##############################


def pool_adjustment_file_preprocessing(
    uvcfolder,
    samplefolder,
    labelfolder,
    imfilename,
):
    """Generate adjustment file to impose heterogeneity in cell-scale dynamics (imp_region), for each sample.

    Parameters
    ----------
    uvcfolder : str
        Path to folder with universal ventricular coordinates.
    samplefolder: str
        Path to folder with all sample(s) details.
    labelfolder: str
        Path to the ionic model parameters labelfile.
    imfilename: str
        Name of .txt file with labels of all ionic model parameters (without extension).
    """

    # check folder structure
    print("\nChecking folder structure...")
    to_check = [
        uvcfolder + "/RHO.dat",
        uvcfolder + "/Z.dat",
        labelfolder + "/" + imfilename + ".txt",
    ]
    for f in to_check:
        if not os.path.exists(f):
            raise Exception("Cannot find {}. Check path to file.".format(f))
        else:
            print(f + " found.")

    # load UVC
    RHO = IO.read_data(uvcfolder + "/RHO")
    Z = IO.read_data(uvcfolder + "/Z")

    # find spatially varying parameter(s) using first sample folder
    par_keys = []
    impar = list(IO.read_text(labelfolder + "/" + imfilename))
    folderlist = [
        d
        for d in os.listdir(samplefolder)
        if os.path.isdir(os.path.join(samplefolder, d))
    ]
    for _, key in enumerate(
        IO.read_json(samplefolder + "/{}/param".format(folderlist[-1]))
    ):
        if any(k in key for k in impar):
            if "_theta" in key or "_grad" in key:
                par_keys.append(key)

    if not par_keys:
        print("There are no spatially varying ionic model parameters! Skipping.\n")
        compute_adj = False
        return compute_adj, [], [], [], [], [], []

    else:
        compute_adj = True
        ion_var_list = [
            var.replace("_theta", "") for var in par_keys if "_theta" in var
        ]

        # extract node indices for each region (can be done only once)
        vtx_dct = find_biv_vtx(uvcfolder)

        # find accepted parameters
        readmefile = samplefolder + "/accepted_samples"
        if not os.path.exists(readmefile + ".txt"):
            readme = [str(i) for i in os.listdir(samplefolder) if not ".txt" in i]
        else:
            readme = IO.read_text(readmefile)
            readme = readme[1:]

        return (
            compute_adj,
            RHO,
            Z,
            par_keys,
            ion_var_list,
            vtx_dct,
            readme,
        )


def pool_adjustment_file_sample(
    id_sample,
    samplefolder,
    vtx_dct,
    RHO,
    Z,
    par_keys,
    ion_var_list,
    domain,
):
    """Generate adjustment file to impose heterogeneity in cell-scale dynamics (imp_region), for one sample.

    Parameters
    ----------
    id_sample: int
        ID of current sample.
    samplefolder: str
        Path to folder with all sample(s) details.
    vtx_dct: dict
        List of vertices for each biventricular sub-region.
    RHO: list
        Rotational UVC coordinate.
    Z: list
        Apico-basal UVC coordinate.
    par_keys: list
        Spatially varying parameter(s).
    ion_var_list: list
        Ionic variables in homogeneous biventricle.
    domain: str
        Definition domain (either "intra" or "extra").
    """

    # current simulation paramfile
    paramfile = samplefolder + "/" + str(id_sample) + "/param"
    if not os.path.exists(paramfile + ".json"):
        raise Exception("Cannot find {}.".format(paramfile + ".json"))

    # current sample output directory
    SIMID = samplefolder + "/" + str(id_sample)
    if not os.path.exists(SIMID):
        os.mkdir(SIMID)

    # output directory for adjustment file(s)
    adjfolder = SIMID + "/ADJ"
    if not os.path.exists(adjfolder) or len(os.listdir(adjfolder)) < len(ion_var_list):
        os.mkdir(adjfolder)

        print("Writing adjustment files for " + str(id_sample) + "-th paramfile...")

        # load dictionary with all current simulation parameters
        par_all_dct = IO.read_json(paramfile)

        # find values for spatially varying parameter(s)
        par_dct = {}
        for key in par_keys:
            par_dct[key] = par_all_dct[key]

        # write adjustment file for each spatially varying ionic variable
        for ion_var in ion_var_list:
            theta_dct = {"biv": par_dct[ion_var + "_theta"]}
            grad_rho_dct = {"biv": par_dct[ion_var + "_grad_rho"]}
            grad_z_dct = {"biv": par_dct[ion_var + "_grad_z"]}

            # adj-file generation
            adjustment_file(
                adjfolder,
                ion_var,
                vtx_dct,
                RHO,
                Z,
                theta_dct,
                grad_rho_dct,
                grad_z_dct,
                domain=domain,
                visualise=False,
            )
    else:
        print(
            "Adjustment files for "
            + str(id_sample)
            + "-th sample already computed! Skipping."
        )

    return


def pool_CV_adjustment_file_preprocessing(
    CARP_INSTALL_DIR,
    meshfolder,
    meshname,
    uvcfolder,
    samplefolder,
    _RVmyo="34",
    _LVmyo="35",
    _SE_RV="44",
    _SE_LV="45",
    _SE_septRV="46",
):
    """Generate datafiles to impose heterogeneity in conductivity region (gregion), for each sample.

    Parameters
    ----------
    CARP_INSTALL_DIR: str
        Path to CARPentry installation directory.
    meshfolder: str
        Path to the computational mesh.
    meshname: str
        Name of the torso computational mesh.
    uvcfolder : str
        Path to folder with universal ventricular coordinates.
    samplefolder: str
        Path to folder with all sample(s) details.
    _RVmyo: str
        Pre-existing tag(s) for RV myocardium. Default: "34".
    _LVmyo: str
        Pre-existing tag(s) for LV myocardium. Default: "35".
    _SE_RV: str
        Pre-existing tag(s) for RV subendocardial layer. Default: "44".
    _SE_LV: str
        Pre-existing tag(s) for LV subendocardial layer. Default: "45".
    _SE_septRV: str
        Pre-existing tag(s) for septal subendocardial layer on the right side. Default: "46".
    """

    # find executables
    meshtoolloc = CARP_INSTALL_DIR + "/meshtool/meshtool"

    # check if vel_f is a spatially varying parameter using first sample folder
    par_keys = []
    folderlist = [
        d
        for d in os.listdir(samplefolder)
        if os.path.isdir(os.path.join(samplefolder, d))
    ]
    for _, key in enumerate(
        IO.read_json(samplefolder + "/{}/param".format(folderlist[-1]))
    ):
        if "vel_f" in key:
            if "_theta" in key or "_grad" in key:
                par_keys.append(key)

    # folder with temporary outputs
    gregfolder = samplefolder + "/../GREGION"

    if not par_keys:
        print("vel_f is not a spatially varying parameter.")
        compute_CV = False
        return compute_CV, [], [], [], [], [], [], [], []

    else:
        compute_CV = True

        # generate temporary output folder
        if not os.path.exists(gregfolder):
            os.mkdir(gregfolder)

        # check folder structure
        print("\nChecking folder structure...")
        to_check = [
            uvcfolder + "/RHO.dat",
            uvcfolder + "/Z.dat",
            uvcfolder + "/Z_subendocardium.dat",
        ]
        for f in to_check:
            if not os.path.exists(f):
                raise Exception("Cannot find {}. Check path to file.".format(f))
            else:
                print(f + " found.")

        # find accepted parameters
        readmefile = samplefolder + "/accepted_samples"
        if not os.path.exists(readmefile + ".txt"):
            readme = [str(i) for i in os.listdir(samplefolder) if not ".txt" in i]
        else:
            readme = IO.read_text(readmefile)
            readme = readme[1:]

        # compute conductivity regions only if at least one CV datafile is missing
        if all(
            [os.path.exists(samplefolder + "/" + str(i) + "/CV") for i in readme]
        ) and all(
            [len(os.listdir(samplefolder + "/" + str(i) + "/CV")) == 3 for i in readme]
        ):
            print("CV datafiles already computed for all accepted samples! Skipping.\n")
            return compute_CV, [], [], [], [], [], [], [], []

        else:
            # extract biv from torso (including ALL subregions, either in the myocardium or in the subendocardial layer)
            meshpath_i = gregfolder + "/" + meshname + "_i"
            if not meshtool.extract_mesh(
                meshtoolloc,
                meshfolder + "/" + meshname,
                _RVmyo + "," + _LVmyo + "," + _SE_RV + "," + _SE_LV + "," + _SE_septRV,
                meshpath_i,
            ):
                raise Exception("meshtool error. Abort!")
            elem_biv = IO.read_elem(meshpath_i)

            # load UVC
            RHO = IO.read_data(uvcfolder + "/RHO")
            Z = IO.read_data(uvcfolder + "/Z")
            Z_SElayer = IO.read_data(uvcfolder + "/Z_subendocardium")

            # extract vertices for each region
            vtx_dct = find_gregion_vtx(
                elem_biv, _RVmyo, _LVmyo, _SE_RV, _SE_LV, _SE_septRV
            )

            # for each region, compute mean UVC
            meanZ_dct = {}
            meanRHO_dct = {}
            meanZ_SE_dct = {}
            meanRHO_SE_dct = {}
            for _, subreg in enumerate(vtx_dct):
                if "SE_" in subreg:
                    meanZ_SE_dct[subreg.replace("SE_", "")] = {}
                    meanRHO_SE_dct[subreg.replace("SE_", "")] = {}
                    for _, subreg_tag in enumerate(vtx_dct[subreg]):
                        meanZ_SE_dct[subreg.replace("SE_", "")][subreg_tag] = np.mean(
                            Z_SElayer[vtx_dct[subreg][subreg_tag]]
                        )
                        meanRHO_SE_dct[subreg.replace("SE_", "")][subreg_tag] = 0
                else:
                    meanZ_dct[subreg] = {}
                    meanRHO_dct[subreg] = {}
                    for _, subreg_tag in enumerate(vtx_dct[subreg]):
                        meanZ_dct[subreg][subreg_tag] = np.mean(
                            Z[vtx_dct[subreg][subreg_tag]]
                        )
                        meanRHO_dct[subreg][subreg_tag] = np.mean(
                            RHO[vtx_dct[subreg][subreg_tag]]
                        )

            return (
                compute_CV,
                meshpath_i,
                vtx_dct,
                meanRHO_dct,
                meanZ_dct,
                meanRHO_SE_dct,
                meanZ_SE_dct,
                par_keys,
                readme,
            )


def pool_CV_adjustment_file_sample(
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
    """Generate datafiles to impose heterogeneity in conductivity region (gregion), for one sample.

    Parameters
    ----------
    id_sample: int
        ID of current sample.
    samplefolder: str
        Path to folder with all sample(s) details.
    meshpath_i: str
        Path to biventricular mesh and sub-regions, if needed.
    vtx_dct: dict
        List of vertices for each biventricular sub-region.
    meanRHO_dct: dict
        Mean rotational UVC coordinate per element.
    meanZ_dct: dict
        Mean apico-basal UVC coordinate per element.
    meanRHO_SE_dct: dict
        Mean rotational UVC coordinate per element in subendocardial layer.
    meanZ_SE_dct: dict
        Mean apico-basal UVC coordinate per element in subendocardial layer.
    par_keys: list
        Spatially varying parameter(s).
    """

    # current simulation paramfile
    paramfile = samplefolder + "/" + str(id_sample) + "/param"
    if not os.path.exists(paramfile + ".json"):
        raise Exception("Cannot find {}.".format(paramfile + ".json"))

    # current sample output directory
    SIMID = samplefolder + "/" + str(id_sample)
    if not os.path.exists(SIMID):
        os.mkdir(SIMID)

    # output directory for adjustment-like datafile(s)
    cvfolder = SIMID + "/CV"
    if not os.path.exists(cvfolder) or len(os.listdir(cvfolder)) < 3:
        os.mkdir(cvfolder)

        print("Writing CV datafiles for " + str(id_sample) + "-th paramfile...")

        # load dictionary with all current simulation parameters
        par_all_dct = IO.read_json(paramfile)

        # find values for spatially varying parameter(s)
        par_dct = {}
        for key in par_keys:
            par_dct[key] = par_all_dct[key]

        # keep only indices for fast conducting regions in SElayer, depending on Zmax_SE
        vtx_red_dct = update_FEC_vtx(vtx_dct, par_all_dct["Zmax_SE"])

        idx_dct = {}
        idx_SE_dct = {}
        for _, subreg in enumerate(vtx_red_dct):
            if "SE_" in subreg:
                idx_SE_dct[subreg.replace("SE_", "")] = vtx_red_dct[subreg]
            else:
                idx_dct[subreg] = vtx_red_dct[subreg]

        # get anisotropy ratios
        ratio = {}
        ratio["vel_f"] = 1
        ratio["vel_s"] = par_all_dct["aniso_ratio"]
        ratio["SE_vel_f"] = par_all_dct["CV_SE_ratio"]

        # write adjustment file for vel_f, vel_s = aniso_ratio * vel_f and vel_f_SE = CV_SE_ratio * vel_f
        for cv_var in ["vel_f", "vel_s", "SE_vel_f"]:
            theta_dct = {"biv": ratio[cv_var] * par_dct["vel_f_theta"]}
            grad_rho_dct = {"biv": par_dct["vel_f_grad_rho"]}
            grad_z_dct = {"biv": par_dct["vel_f_grad_z"]}

            # dat-file generation
            if cv_var == "SE_vel_f":
                adjustment_gregion(
                    cvfolder,
                    cv_var,
                    idx_SE_dct,
                    meanRHO_SE_dct,
                    meanZ_SE_dct,
                    theta_dct,
                    grad_rho_dct,
                    grad_z_dct,
                    meshpath_i,
                    visualise=False,
                )

            else:
                adjustment_gregion(
                    cvfolder,
                    cv_var,
                    idx_dct,
                    meanRHO_dct,
                    meanZ_dct,
                    theta_dct,
                    grad_rho_dct,
                    grad_z_dct,
                    meshpath_i,
                    visualise=False,
                )
    else:
        print(
            "CV datafiles for "
            + str(id_sample)
            + "-th sample already computed! Skipping."
        )

    return
