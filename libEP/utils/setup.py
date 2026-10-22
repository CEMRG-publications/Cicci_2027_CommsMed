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
from utils import processing_mesh as proc

# extra imports
import os
import shutil
import numpy as np
import pandas as pd
import scipy


def copy_mesh(
    CARP_INSTALL_DIR,
    inputfolder,
    imeshname,
    outputfolder,
    omeshname,
    fmt_list=["elem", "lon", "pts"],
    extract_myo=False,
):
    """Copy mesh in specified format.

    Parameters
    ----------
    CARP_INSTALL_DIR: str
        Path to CARPentry installation directory.
    inputfolder: str
        Path to folder with input mesh to be copied.
    imeshname : str
        Input meshname.
    outputfolder: str
        Path to destination folder.
    omeshname : str
        Output meshname.
    fmt_list: list
        Mesh format. Default: ["elem","lon","pts"] (i.e., carp_txt).
    extract_myo: bool
        If True, extract myocardium from mesh. Default: False.
    """

    if type(fmt_list) == str:
        fmt_list = [fmt_list]

    for fmt in fmt_list:
        shutil.copyfile(
            inputfolder + "/" + imeshname + "." + fmt,
            outputfolder + "/" + omeshname + "." + fmt,
        )

    if extract_myo:
        # find executables
        meshtoolloc = CARP_INSTALL_DIR + "/meshtool/meshtool"

        # save tag info in a text file
        tag_file = outputfolder + "/info_mesh"
        meshtool.query_tags(meshtoolloc, outputfolder + "/" + omeshname, tag_file)

        # extract myocardium tag(s) only
        with open(tag_file + ".txt", "r") as file:
            lines = [line.rstrip() for line in file]
        init_myo_line = [l for l in range(len(lines)) if "Myocardium" in lines[l]][
            0
        ] + 1

        meshtool.extract_mesh(
            meshtoolloc,
            outputfolder + "/" + omeshname,
            lines[init_myo_line],
            outputfolder + "/" + omeshname + "_i",
        )

        # delete tag_file
        os.remove(tag_file + ".txt")
        os.remove(outputfolder + "/" + omeshname + "_i.eidx")
        os.remove(outputfolder + "/" + omeshname + "_i.nod")

    return


def copy_par_file(
    inputfolder, outputfolder, par_list=["pfile_BASE", "pfile_ECG", "pfile_LF"]
):
    """Copy simulation .par file(s).

    Parameters
    ----------
    inputfolder: str
        Path to folder with .par file(s) to be copied.
    outputfolder: str
        Path to destination folder.
    par_list : list
        Input simulation .par file(s). Default: ["pfile_BASE","pfile_ECG","pfile_LF"].
    """

    if type(par_list) == str:
        par_list = [par_list]

    for par in par_list:
        shutil.copyfile(
            inputfolder + "/" + par + ".par", outputfolder + "/" + par + ".par"
        )

    return


def copy_electrode_file(inputfolder, filename, outputfolder):
    """Copy electrode .pts file.

    Parameters
    ----------
    inputfolder: str
        Path to folder with .pts file to be copied.
    filename: str
        Name of electrode .pts file.
    outputfolder: str
        Path to destination folder.
    """

    shutil.copyfile(inputfolder + "/" + filename + ".pts", outputfolder + "/lf_src.pts")

    return


def copy_lats_file(inputfolder, filename, outputfolder):
    """Copy LATs .dat file.

    Parameters
    ----------
    inputfolder: str
        Path to folder with .dat file to be copied.
    filename: str
        Name of LATs .dat file.
    outputfolder: str
        Path to destination folder.
    """

    shutil.copyfile(
        inputfolder + "/" + filename + ".dat", outputfolder + "/" + filename + ".dat"
    )

    return


def check_sample_folder(
    samplefolder,
    to_check=[
        "HPS/LV_af.vtx",
        "HPS/LV_pf.vtx",
        "HPS/LV_sf.vtx",
        "HPS/RV_sf.vtx",
        "HPS/RV_mod.vtx",
    ],
):

    # find .txt file with accepted parameters. If already exists, rename it.
    readmefile = samplefolder + "/accepted_samples"
    if not os.path.exists(readmefile + ".txt"):
        readme = ["nan"] + [str(i) for i in os.listdir(samplefolder)]
    else:
        readme = IO.read_text(readmefile)

    # loop on each sample
    accepted_samples = "accepted_samples\n"
    print("\nChecking folder structure...")
    for i in os.listdir(samplefolder):
        if not ".txt" in i and str(i) in readme[1:]:
            incomplete = False
            for f in to_check:
                if not os.path.exists(samplefolder + "/" + str(i) + "/" + f):
                    print("Cannot find {} for sample {}.".format(f, i))
                    incomplete = True
            if not incomplete:
                accepted_samples += str(i) + "\n"
    IO.write_string_text(readmefile, accepted_samples)

    return


def electrodes_pts2vtx(meshfolder, meshname, electfolder, electfilename, outputfolder):
    """Find electrode vertices corresponding to given coordinates.

    Parameters
    ----------
    meshfolder: str
        Path to the computational mesh.
    meshname: str
        Name of the torso computational mesh.
    electfolder: str
        Path to file with electrode euclidean coordinates.
    electfilename:
        Name of .pts file with electrode euclidean coordinates (without extension).
    outputfolder: str
        Path to the directory where the sample datafile(s) will be saved.
    """

    # pts of torso mesh
    pts = IO.read_pts(meshfolder + "/" + meshname)

    # number of elem of torso mesh
    num_elem = len(IO.read_elem(meshfolder + "/" + meshname))

    # pts of electrodes
    subpts = IO.read_pts(electfolder + "/" + electfilename)

    # find electrode IDs in torso mesh
    vtx = np.zeros((np.size(subpts, 0), 1))
    for i in range(np.size(subpts, 0)):
        index = np.where(
            (abs(pts.loc[:, 0] - subpts.loc[i, 0]) < 0.1)
            & (abs(pts.loc[:, 1] - subpts.loc[i, 1]) < 0.1)
            & (abs(pts.loc[:, 2] - subpts.loc[i, 2]) < 0.1)
        )
        vtx[i] = int(index[0])

    np.savetxt(
        outputfolder + "/lf_src.vtx",
        vtx,
        fmt="%d",
        comments="",
        header=str(np.size(subpts, 0))
        + " # mesh extra "
        + str(int(num_elem))
        + " "
        + str(np.size(pts, 0)),
    )

    # save LF reference electrode in a separate file (needed when using bidomain source model)
    np.savetxt(
        outputfolder + "/lf_src_ref.vtx",
        np.array([vtx[-2]]),
        fmt="%d",
        comments="",
        header="1\nextra",
    )

    return [int(i) for i in vtx]


# TODO
# use template files rather than "updating". See write_slrm in utils/setup.py.
def update_par_lead_field(paramfolder, parfilename, lf_ref_vtx, clusterfolder):
    """Update .par file with carp command to compute lead field arrays.

    Parameters
    ----------
    paramfolder: str
        Path to the command parameter file.
    paramfilename: str
        Name of .par file with simulation commands (without extension).
    lf_ref_vtx: int
        Vertex of the reference ground (usually electrode RL).
    clusterfolder: str
        Path to folder where simulations will be run (usually on HPC facilities).
    """

    print("LEAD FIELD computation setting...")

    # append to par-file to compute lead field matrices
    with open(paramfolder + "/" + parfilename + ".par", "a") as file:
        file.write("lead_field_ref     = {}\n".format(str(lf_ref_vtx)))
        file.write("lead_field_sources = {}\n".format(clusterfolder + "/lf_src"))
        file.close()

    return


# TODO
# use template files rather than "updating". See write_slrm in utils/setup.py.
def update_par_ECG(paramfolder, parfilename, lf_ref_vtx, clusterfolder):
    """Update .par file with carp command to run electrophysiology simulations.

    Parameters
    ----------
    paramfolder: str
        Path to the command parameter file.
    paramfilename: str
        Name of .par file with simulation commands (without extension).
    lf_ref_vtx: int
        Vertex of the reference ground (usually electrode RL).
    clusterfolder: str
        Path to folder where simulations will be run (usually on HPC facilities).
    """

    print("ECG computing setting...")

    # append to par-file to compute ECG using the lead field approach
    with open(paramfolder + "/" + parfilename + ".par", "a") as file:
        file.write("lead_field_ref     = {}\n".format(str(lf_ref_vtx)))
        file.write("lead_field_sources = {}\n".format(clusterfolder + "/lf_src"))
        # file.write("phie_rec_ptf       = {}\n".format(clusterfolder + "/lf_src"))
        file.close()

    return


# --- ARCHER2 JOB SCRIPTS --- #




# --- EP MODEL --- #
def write_par_EPmodel(
    samplefolder,
    clusterfolder_samples,
    clusterfolder_data,
    propagation="RE+",
    sourcemodel="monodomain",
    parfile_ep="pfile_EP",
    eas_list=["LV_sf", "LV_pf", "LV_af", "RV_sf", "RV_mod"],
    hps_ref=False,
):
    """Write .par file to run electrophysiology simulations, for each sample.

    Parameters
    ----------
    samplefolder: str
        Path to folder with all sample(s) details (usually, local HD).
    clusterfolder_samples: str
        Path to folder will all sample(s) details (usually, HPC or 64-core machines)
    clusterfolder_data: str
        Path to folder will additional simulation data (usually, HPC or 64-core machines)
    propagation: str
        EP propagation method (e.g., "RE-", "RE+"). Default: "RE+".
    sourcemodel: str
        EP source model (e.g., "monodomain", "pseudobidomain", "bidomain"). Default: "monodomain".
    parfile_ep: str
        Basename of .par file with basic commands to run EP simulations (without extension). Default: "pfile_EP".
    eas_list: list
        Early activation site name(s). Default: ["LV_sf", "LV_pf", "LV_af", "RV_sf", "RV_mod"].
    hps_ref: bool
        If True, use .vtx file(s) with location of EAS(s) common among all samples (without extension), if any. Default: False.
    """

    # find .txt file with accepted parameters
    readmefile = samplefolder + "/accepted_samples"
    if not os.path.exists(readmefile + ".txt"):
        readme = ["nan"] + [str(i) for i in os.listdir(samplefolder)]
    else:
        readme = IO.read_text(readmefile)

    for i in os.listdir(samplefolder):
        if not ".txt" in i and str(i) in readme[1:]:
            print("Writing EP simulation setting for " + str(i) + "-th sample...")

            # write par-file to compute lead field matrices
            outputfile = samplefolder + "/" + str(i) + "/" + parfile_ep
            with open(outputfile + ".par", "w") as file:
                # current simulation paramfile
                par_file = samplefolder + "/" + str(i) + "/param"
                if not os.path.exists(par_file + ".json"):
                    raise Exception(
                        "Cannot find {}. Check path to file.".format(par_file + ".json")
                    )
                par_dct = IO.read_json(par_file)

                file.write(
                    "# ------------------------------- {} {} -------------------------------#\n".format(
                        propagation, sourcemodel
                    )
                )

                file.write(
                    "\n# -------------------------------\n# EP MODEL\n# -------------------------------\n"
                )

                # set source model type
                if sourcemodel == "pseudo_bidomain":
                    file.write("bidomain    = 2\n")
                elif sourcemodel == "bidomain":
                    file.write("bidomain    = 1\n")
                else:
                    file.write("bidomain    = 0\n")

                # set propagation model
                if propagation == "RE+":
                    file.write("diffusionOn = 1\n")
                else:
                    file.write("diffusionOn = 0\n")

                file.write(
                    "\n# -------------------------------\n# STIMULATION\n# -------------------------------\n"
                )

                # stimulus basename(s)
                stim_param = eas_list

                # total number of stimuli
                totNstims = len(stim_param) + 1
                if not sourcemodel == "monodomain":
                    totNstims += 1
                file.write("num_stim             = {}\n".format(str(totNstims)))

                # set up stimuli
                Nstims = len(stim_param)
                for j in range(Nstims):
                    file.write("stimulus[{}].name     = {}\n".format(j, stim_param[j]))
                    file.write("stimulus[{}].stimtype = 0\n".format(j))
                    file.write("stimulus[{}].duration = 1.0\n".format(j))
                    file.write("stimulus[{}].strength = 100\n".format(j))
                    file.write(
                        "stimulus[{}].start    = {}\n".format(
                            j, str(par_dct[stim_param[j] + "_t"])
                        )
                    )
                    
                    # write path to lon file
                    if hps_ref:
                        file.write(
                            "stimulus[{}].vtx_file = {}\n\n".format(
                                j,
                                clusterfolder_samples
                                + "/../ref_QRS/0/HPS/"
                                + stim_param[j],
                            )
                        )
                    else:
                        file.write(
                            "stimulus[{}].vtx_file = {}\n\n".format(
                                j,
                                clusterfolder_samples
                                + "/"
                                + str(i)
                                + "/HPS/"
                                + stim_param[j],
                            )
                        )

                # set additional stimulus (ground) for pseudo_bidomain and bidomain models (i.e., leadfield reference must be grounded)
                if not sourcemodel == "monodomain":
                    Nstims = Nstims + 1
                    file.write(
                        "stimulus[{}].name     = LF_REF_GND\n".format(str(Nstims - 1))
                    )
                    file.write("stimulus[{}].stimtype = 3\n".format(str(Nstims - 1)))
                    file.write("stimulus[{}].start    = 0\n".format(str(Nstims - 1)))
                    file.write(
                        "stimulus[{}].vtx_file = {}\n\n".format(
                            str(Nstims - 1), clusterfolder_data + "/lf_src_ref"
                        )
                    )

                # set additional stimulus (eikonal) for RE models
                Nstims = Nstims + 1
                file.write("stimulus[{}].stimtype = 8".format(str(Nstims - 1)))

                file.close()

    return


# --- FIBRE --- #
def write_par_fibre(
    samplefolder,
    clusterfolder_samples,
    lonfilename,
    parfile_tissue="pfile_TISSUE",
    LF_ref=False,
):
    """Write fibre section in .par file, for each sample.

    Parameters
    ----------
    samplefolder: str
        Path to folder with all sample(s) details (usually, local HD).
    clusterfolder_samples: str
        Path to folder will all sample(s) details (usually, HPC or 64-core machines)
    lonfilename: str
        Name of .lon file with fibre orientation.
    parfile_tissue: str
        Basename of .par file with basic commands to adjust tissue properties (without extension). Default: "pfile_TISSUE".
    LF_ref: bool
        If True, read fibres from common file. Default: False.
    """

    # find .txt file with accepted parameters
    readmefile = samplefolder + "/accepted_samples"
    if not os.path.exists(readmefile + ".txt"):
        readme = ["nan"] + [str(i) for i in os.listdir(samplefolder)]
    else:
        readme = IO.read_text(readmefile)

    # define fibres using UVC for each sample
    for i in os.listdir(samplefolder):
        if not ".txt" in i and str(i) in readme[1:]:
            print("FIBRE setting for " + str(i) + "-th sample...")

            # open existing parameter file
            outputfile = samplefolder + "/" + str(i) + "/" + parfile_tissue
            parf = open(outputfile + ".par", "a")
            parf.write(
                "\n# -------------------------------\n# FIBRES\n# -------------------------------\n"
            )

            if LF_ref:
                if ".lon" in lonfilename:
                    parf.write(
                        "orthoname = {}\n".format(
                            clusterfolder_samples
                            + "/0/"
                            + lonfilename.replace(".lon", "")
                        )
                    )
                # TODO: check if it works properly
                elif ".blon" in lonfilename:
                    parf.write(
                        "orthoname = {}\n".format(
                            clusterfolder_samples
                            + "/0/"
                            + lonfilename.replace(".blon", "")
                        )
                    )
            else:
                if ".lon" in lonfilename:
                    parf.write(
                        "orthoname = {}\n".format(
                            clusterfolder_samples
                            + "/"
                            + str(i)
                            + "/"
                            + lonfilename.replace(".lon", "")
                        )
                    )
                # TODO: check if it works properly
                elif ".blon" in lonfilename:
                    parf.write(
                        "orthoname = {}\n".format(
                            clusterfolder_samples
                            + "/"
                            + str(i)
                            + "/"
                            + lonfilename.replace(".blon", "")
                        )
                    )

            parf.close()

    return


# --- IMP REGION --- #
def write_par_imp_region(
    CARP_INSTALL_DIR,
    meshfolder,
    meshname,
    labelfilename_im,
    labelfolder,
    labelfilename_scalar,
    samplefolder,
    clusterfolder_samples,
    im="ToRORd_dynCl",
    im_sv_init="",
    prepacing_lats="",
    prepacing_beats=100,
    prepacing_bcl=1000,
    parfile_tissue="pfile_TISSUE",
):
    """Write ionic model section in .par file, for each sample.

    Parameters
    ----------
    CARP_INSTALL_DIR: str
        Path to CARPentry installation directory.
    meshfolder: str
        Path to the computational mesh.
    meshname: str
        Name of the torso computational mesh.
    labelfile_im: str
        Name of label file with ionic model parameters.
    labelfolder: str
        Path to the ionic model parameters labelfile.
    labelfilename_scalar: str
        Name of label file with scalar parameters.
    samplefolder: str
        Path to folder with all sample(s) details (usually, local HD).
    clusterfolder_samples: str
        Path to folder will all sample(s) details (usually, HPC or 64-core machines).
    im: str
        Name of ionic model. Default: "ToRORd_dynCl".
    im_sv_init: str
        Name of initial state file for single cell model used by the sv_init method (if any). Default: "".
    prepacing_lats: str
        Name of tissue activation file to guide state set-up used by the prepace method (if any). Default: "".
    prepacing_beats: int
        Number of pulses applied to single cell. Default: 100.
    prepacing_bcl: int
        Basic cycle length for repetitive single cell stimulation. Default: 1000 [ms].
    """

    # find executables
    meshtoolloc = CARP_INSTALL_DIR + "/meshtool/meshtool"

    # check folder structure
    print("\nChecking folder structure...")
    to_check = [
        labelfolder + "/" + labelfilename_im + ".txt",
        labelfolder + "/" + labelfilename_scalar + ".txt",
    ]
    for f in to_check:
        if not os.path.exists(f):
            raise Exception("Cannot find {}. Check path to file.".format(f))
        else:
            print(f + " found.")

    # read all ionic model parameters
    im_par_all = IO.read_text(labelfolder + "/" + labelfilename_im)

    # read all scalar input parameters
    scl_par_all = IO.read_text(labelfolder + "/" + labelfilename_scalar)

    # find .txt file with accepted parameters
    readmefile = samplefolder + "/accepted_samples"
    if not os.path.exists(readmefile + ".txt"):
        readme = ["nan"] + [str(i) for i in os.listdir(samplefolder) if not ".txt" in str(i)]
    else:
        readme = IO.read_text(readmefile)

    # get all parameters from first accepted sample
    par_file = samplefolder + "/" + str(readme[1]) + "/param"
    if not os.path.exists(par_file + ".json"):
        raise Exception(
            "Cannot find {}. Check path to file.".format(par_file + ".json")
        )
    par_dct = IO.read_json(par_file)
    par_keys = list(par_dct.keys())

    # find all scalar ionic model parameters and only spatially varying ones (if any)
    imp_scl_names = []
    imp_grad_names = []
    for imp in im_par_all:
        if imp in par_keys:
            imp_scl_names.append(imp)
        elif any(imp in key for key in par_keys):
            imp_grad_names.append(imp)

    # save tag info in a text file
    tag_file = meshfolder + "/info_mesh"
    meshtool.query_tags(meshtoolloc, meshfolder + "/" + meshname, tag_file)

    # extract myocardium tag(s) only
    with open(tag_file + ".txt", "r") as file:
        lines = [line.rstrip() for line in file]
    init_myo_line = [l for l in range(len(lines)) if "Myocardium" in lines[l]][0] + 1
    tag_myo = [int(tag) for tag in lines[init_myo_line].split(",")]

    # delete tag_file
    os.remove(tag_file + ".txt")

    for i in os.listdir(samplefolder):
        if not ".txt" in i and str(i) in readme[1:]:
            print("IMP_REGION setting for " + str(i) + "-th sample...")

            # current simulation paramfile
            par_file = samplefolder + "/" + str(i) + "/param"
            if not os.path.exists(par_file + ".json"):
                raise Exception(
                    "Cannot find {}. Check path to file.".format(par_file + ".json")
                )
            par_dct = IO.read_json(par_file)

            # find all scalar ionic model parameters (if any)
            imp_scl_dct = {}
            for imp in imp_scl_names:
                if imp in scl_par_all:
                    imp_scl_dct[imp] = par_dct[imp]

            # open existing parameter file
            outputfile = samplefolder + "/" + str(i) + "/" + parfile_tissue
            parf = open(outputfile + ".par", "a")
            parf.write(
                "\n# -------------------------------\n# IONIC MODEL\n# -------------------------------\n"
            )

            parf.write("num_imp_regions        = 1\n\n")
            parf.write(
                "imp_region[0].im       = {}\n".format(im)
            )  # we assume the same ionic model in each region
            parf.write("imp_region[0].num_IDs  = {}\n".format(len(tag_myo)))
            if im == "ToRORd_dynCl":
                im_param_string = 'imp_region[0].im_param = "flags=ENDO,{}"\n'.format(
                    ",".join(
                        [
                            par + "=" + "{:.4e}".format(imp_scl_dct[par])
                            for par in list(imp_scl_dct.keys())
                        ]
                    )
                )
            else:
                im_param_string = 'imp_region[0].im_param = "{}"\n'.format(
                    ",".join(
                        [
                            par + "=" + "{:.4e}".format(imp_scl_dct[par])
                            for par in list(imp_scl_dct.keys())
                        ]
                    )
                )

            # NOTE: CARPentry does NOT allow more than 256 characters in a string (i.e., imp_region[0].im_param = "...")
            if len(im_param_string) > 256:
                raise Exception(
                    "String for im_param in parameter file is {} characters long. Maximum number of characters allowed by CARP is 256. Abort!".format(
                        len(im_param_string)
                    )
                )
            parf.write(im_param_string)

            # single-cell: sv_init method
            if not imp_grad_names and not im_sv_init == "":
                parf.write(
                    "imp_region[0].im_sv_init = {}\n".format(
                        clusterfolder_samples + "/" + str(i) + "/" + im_sv_init + ".sv"
                    )
                )

            # single-cell: prepace method
            if not prepacing_lats == "" and im_sv_init == "":
                parf.write("prepacing_lats  = {}.dat\n".format(prepacing_lats))
                parf.write("prepacing_beats = {}\n".format(prepacing_beats))
                parf.write("prepacing_bcl   = {}\n".format(prepacing_bcl))

            for j in range(len(tag_myo)):
                parf.write("imp_region[0].ID[{}]    = {}\n".format(j, tag_myo[j]))

            # adjustment files to import heterogeneity nodewise
            if imp_grad_names:
                parf.write("\n\nnum_adjustments = {}\n\n".format(len(imp_grad_names)))
                for j in range(len(imp_grad_names)):
                    parf.write(
                        "adjustment[{}].file     = {}\n".format(
                            j,
                            clusterfolder_samples
                            + "/"
                            + str(i)
                            + "/ADJ/"
                            + imp_grad_names[j]
                            + ".adj",
                        )
                    )
                    parf.write(
                        "adjustment[{}].variable = {}\n".format(
                            j, im + "." + imp_grad_names[j]
                        )
                    )

            parf.close()

    return


# --- GREGION --- #
def gregion_adj_to_dct(adjfile):
    """Convert adjustment-like file for gregion to dictionary.

    Parameters
    ----------
    adjfile: str
        Path to adjustment-like file.
    """

    # read gregion adjustment-like file as dataframe
    adj_df = pd.read_csv(adjfile + ".adj", sep=" ")

    # convert dataframe into dictionary
    adj_dct = {}
    for i in range(len(adj_df["ID"])):
        adj_dct[adj_df["ID"][i]] = adj_df["val"][i]

    return adj_dct


def Bspline_tuneCV(csvfile):
    """Compute B-spline hyperparameters from CV-to-conductivity data pre-computed CARP tune_CV.py.

    Parameters
    ----------
    csvfile: str
        Path to .csv file with listed: CV, gi and ge (column-wise).
    """

    # import datafile with CV_to_g map
    df = pd.read_csv(csvfile + ".csv")

    # extract input (CV) and outputs (gi and ge)
    x = np.array(df["CV"])
    Y = np.vstack((np.array(df["gi"]), np.array(df["ge"]))).T

    # Bspline of degree 3 for interpolation
    b3 = scipy.interpolate.make_interp_spline(x, Y)

    # Bspline of degree 1 for extrapolation
    b1 = scipy.interpolate.make_interp_spline(x, Y, k=1)

    return x, Y, b3, b1


def Bspline(cv, xtrain, Ytrain, b_interp, b_extrap):
    """Map CV to intra/extracellular conductivities using B-spline of degree 1 (for extrapolation) or 3 (for interpolation).

    Parameters
    ----------
    cv: float
        Conduction velocity.
    xtrain: np.array
        Training 1D input matrix (i.e., CV).
    Ytrain: np.array
        Training 2D output matrix (i.e., gi and ge).
    b_interp: scipy.interpolate._bsplines.BSpline object
        Bspline of degree 3 for interpolation.
    b_extrap: scipy.interpolate._bsplines.BSpline object
        Bspline of degree 1 for extrapolation.
    """

    # find index of current cv in datafile
    idx = np.where(xtrain == cv)[0]

    # if current CV in not in datafile, interpolate or extrapolate using splines of order 1 or 3, respectively
    if not len(idx):
        if cv > max(xtrain):
            gi, ge = b_extrap(cv)
        elif cv < min(xtrain):
            gi, ge = abs(b_extrap(cv))
        else:
            gi, ge = b_interp(cv)
    else:
        gi, ge = Ytrain[idx][0]
    # print("CV = {}, gi = {}, ge = {}".format(cv,gi,ge))

    return round(gi, 3), round(ge, 3)


def write_par_gregion(
    CARP_INSTALL_DIR,
    meshfolder,
    meshname,
    samplefolder,
    tag_BATH_jsonfile,
    tune_CV_csvfile,
    _SE_biv_list=[44, 45, 46],
    parfile_tissue="pfile_TISSUE",
    cvfolder_ref="",
    tag_HEART_jsonfile="",
):
    """Write tissue conductivity section in .par file, for each sample.

    Parameters
    ----------
    CARP_INSTALL_DIR: str
        Path to CARPentry installation directory.
    meshfolder: str
        Path to the computational mesh.
    meshname: str
        Name of the torso computational mesh.
    samplefolder: str
        Path to folder with all sample(s) details (usually, local HD).
    tag_BATH_jsonfile: str
        Path to .json file with bath tag(s) for each bath region.
    tune_CV_csvfile: str
        Path to .csv file with CV, gi and ge training data.
    _SE_biv_list: list
        Pre-existing tag(s) for biventricular subendocardial layer. Default: [44,45,46].
        NOTE: Useful only when SElayer tags have the same number of digits as the myocardium (e.g., for "_subendocardial" meshes).
    parfile_tissue: str
        Basename of .par file with basic commands to adjust tissue properties (without extension). Default: "pfile_TISSUE".
    cvfolder_ref: str
        Path to .adj files with CVs common among all samples (without extension), if any. Default: "".
    tag_HEART_jsonfile: str
        Path to .json file with the tags of the heart regions (without extension). Needed when the conduction velocity is not spatially varying, to tell myocardium from fast endocardial layer. Default: "".
    """

    # find executables
    meshtoolloc = CARP_INSTALL_DIR + "/meshtool/meshtool"

    if type(_SE_biv_list) == int:
        _SE_biv_list = [_SE_biv_list]

    # check folder structure
    print("\nChecking folder structure...")
    to_check = [tag_BATH_jsonfile + ".json", tune_CV_csvfile + ".csv"]
    for f in to_check:
        if not os.path.exists(f):
            raise Exception("Cannot find {}. Check path to file.".format(f))
        else:
            print(f + " found.")

    # find .txt file with accepted parameters
    readmefile = samplefolder + "/accepted_samples"
    if not os.path.exists(readmefile + ".txt"):
        readme = ["nan"] + [str(i) for i in os.listdir(samplefolder) if not ".txt" in str(i)]
    else:
        readme = IO.read_text(readmefile)

    # get all parameters from first accepted sample
    par_file = samplefolder + "/" + str(readme[1]) + "/param"
    if not os.path.exists(par_file + ".json"):
        raise Exception(
            "Cannot find {}. Check path to file.".format(par_file + ".json")
        )
    par_dct = IO.read_json(par_file)

    # get different bath conductivities and conduction velocity parameters
    bath_names = []
    for _, key in enumerate(par_dct):
        if "_bath" in key:
            bath_names.append(key.replace("g_bath_", ""))

    # read traning CV and conductivities, and compute splines
    xtrain, Ytrain, b_interp, b_extrap = Bspline_tuneCV(tune_CV_csvfile)

    # get bath tag(s)
    bath_tag_dct = IO.read_json(tag_BATH_jsonfile)
    bath_tag_keys = list(bath_tag_dct.keys())
    num_bath_region = len(bath_tag_keys)

    # get myocardium tag(s) (including subendocardial layer)
    tag_file = meshfolder + "/info_mesh"
    meshtool.query_tags(meshtoolloc, meshfolder + "/" + meshname, tag_file)
    with open(tag_file + ".txt", "r") as file:
        lines = [line.rstrip() for line in file]
    os.remove(tag_file + ".txt")

    init_myo_line = [l for l in range(len(lines)) if "Myocardium" in lines[l]][0] + 1
    myo_SE_tag_keys = [int(tag) for tag in lines[init_myo_line].split(",")]
    num_myo_with_SE_region = len(myo_SE_tag_keys)

    # check if all bath regions are in the parameter files
    for key in bath_names:
        if not key in bath_tag_keys:
            raise Exception(
                "Missing tag for {}. Update {}.".format(key, tag_BATH_jsonfile)
            )

    # define fibres using UVC for each sample
    for i in os.listdir(samplefolder):
        if not ".txt" in i and str(i) in readme[1:]:
            print("GREGION and EKREGION setting for " + str(i) + "-th sample...")

            # current simulation paramfile
            par_file = samplefolder + "/" + str(i) + "/param"
            if not os.path.exists(par_file + ".json"):
                raise Exception(
                    "Cannot find {}. Check path to file.".format(par_file + ".json")
                )
            par_dct = IO.read_json(par_file)

            # get different bath conductivities and conduction velocity parameters
            bath_dct = {}
            cv_dct = {}
            for _, key in enumerate(par_dct):
                if "_bath" in key:
                    bath_dct[key.replace("g_bath_", "")] = par_dct[key]
                elif "vel_" in key or "_ratio" in key:
                    cv_dct[key] = par_dct[key]

            # ----- CONDUCTIVITY REGIONS ----- #
            # open existing parameter file
            outputfile = samplefolder + "/" + str(i) + "/" + parfile_tissue
            parf = open(outputfile + ".par", "a")
            parf.write(
                "\n# -------------------------------\n# CONDUCTIVITIES\n# -------------------------------\n"
            )

            # number of conductivity regions
            parf.write(
                "num_gregions       = {}\n\n".format(
                    num_bath_region + num_myo_with_SE_region
                )
            )

            if any("_theta" in key for key in list(cv_dct.keys())):
                
                if cvfolder_ref:
                    # get tag and corresponding value for vel_f, vel_s and SE_vel_f
                    vel_f_dct = gregion_adj_to_dct(
                        cvfolder_ref + "/vel_f"
                    )
                    vel_s_dct = gregion_adj_to_dct(
                        cvfolder_ref + "/vel_s"
                    )
                    SE_vel_f_dct = gregion_adj_to_dct(
                        cvfolder_ref + "/SE_vel_f"
                    )
                else:
                    # get tag and corresponding value for vel_f, vel_s and SE_vel_f
                    vel_f_dct = gregion_adj_to_dct(
                        samplefolder + "/" + str(i) + "/CV/vel_f"
                    )
                    vel_s_dct = gregion_adj_to_dct(
                        samplefolder + "/" + str(i) + "/CV/vel_s"
                    )
                    SE_vel_f_dct = gregion_adj_to_dct(
                        samplefolder + "/" + str(i) + "/CV/SE_vel_f"
                    )

                # get myo tag(s) without fast conductive endocardial layer
                myo_tag_values = list(
                    set(list(vel_f_dct.keys())) - set(list(SE_vel_f_dct.keys()))
                )

                # TODO: check if it works
                # remove SElayer tags from myo_tags manually
                myo_tag_values = list(set(myo_tag_values) - set(_SE_biv_list))

            else:
                # tags of ordinary myocardium (the other tags are fast endocardial layer), depending on current Zmax_SE
                if not tag_HEART_jsonfile:
                    raise Exception(
                        "tag_HEART_jsonfile is needed when the conduction velocity is not spatially varying."
                    )
                myo_tag_values = proc.get_myocardium_tags(
                    myo_SE_tag_keys, par_dct["Zmax_SE"], tag_HEART_jsonfile
                )

            for j, key in enumerate(bath_dct):
                parf.write("gregion[{}].name    = {}\n".format(j, key))
                parf.write(
                    "gregion[{}].num_IDs = {}\n".format(j, len(bath_tag_dct[key]))
                )
                parf.write(
                    "gregion[{}].ID      = {}\n".format(
                        j, ",".join([str(tag) for tag in bath_tag_dct[key]])
                    )
                )
                parf.write("gregion[{}].g_bath  = {}\n\n".format(j, bath_dct[key]))

            # number of bath regions
            # num_bath_region = len(bath_reg)

            # add ventricular gregions
            parf.write("# ----- Ventricles ----- #\n")
            if any("_theta" in key for key in list(cv_dct.keys())):
                for j, tag in enumerate(vel_f_dct):
                    if not tag in myo_tag_values:
                        parf.write(
                            "gregion[{}].name = fast_endocardial_layer\n".format(
                                j + num_bath_region
                            )
                        )
                        parf.write(
                            "gregion[{}].num_IDs = 1\n".format(j + num_bath_region)
                        )
                        parf.write(
                            "gregion[{}].ID      = {}\n".format(
                                j + num_bath_region, tag
                            )
                        )

                        gi, ge = Bspline(
                            SE_vel_f_dct[tag], xtrain, Ytrain, b_interp, b_extrap
                        )
                        parf.write(
                            "gregion[{}].g_il    = {}\n".format(j + num_bath_region, gi)
                        )
                        parf.write(
                            "gregion[{}].g_el    = {}\n".format(j + num_bath_region, ge)
                        )

                        gi, ge = Bspline(
                            0.105 * SE_vel_f_dct[tag],
                            xtrain,
                            Ytrain,
                            b_interp,
                            b_extrap,
                        )
                        parf.write(
                            "gregion[{}].g_it    = {}\n".format(j + num_bath_region, gi)
                        )
                        parf.write(
                            "gregion[{}].g_et    = {}\n".format(j + num_bath_region, ge)
                        )

                    else:
                        parf.write(
                            "gregion[{}].num_IDs = 1\n".format(j + num_bath_region)
                        )
                        parf.write(
                            "gregion[{}].ID      = {}\n".format(
                                j + num_bath_region, tag
                            )
                        )

                        gi, ge = Bspline(
                            vel_f_dct[tag], xtrain, Ytrain, b_interp, b_extrap
                        )
                        parf.write(
                            "gregion[{}].g_il    = {}\n".format(j + num_bath_region, gi)
                        )
                        parf.write(
                            "gregion[{}].g_el    = {}\n".format(j + num_bath_region, ge)
                        )

                        gi, ge = Bspline(
                            vel_s_dct[tag], xtrain, Ytrain, b_interp, b_extrap
                        )
                        parf.write(
                            "gregion[{}].g_it    = {}\n".format(j + num_bath_region, gi)
                        )
                        parf.write(
                            "gregion[{}].g_et    = {}\n".format(j + num_bath_region, ge)
                        )

            else:
                myo_cv = [
                    scaling * cv_dct["vel_f"]
                    for scaling in [1, cv_dct["aniso_ratio"], cv_dct["CV_SE_ratio"]]
                ]
                for j, tag in enumerate(myo_SE_tag_keys):
                    if not tag in myo_tag_values:
                        parf.write(
                            "gregion[{}].name = fast_endocardial_layer\n".format(
                                j + num_bath_region
                            )
                        )
                        parf.write(
                            "gregion[{}].num_IDs = 1\n".format(j + num_bath_region)
                        )
                        parf.write(
                            "gregion[{}].ID      = {}\n".format(
                                j + num_bath_region, tag
                            )
                        )

                        gi, ge = Bspline(myo_cv[2], xtrain, Ytrain, b_interp, b_extrap)
                        parf.write(
                            "gregion[{}].g_il    = {}\n".format(j + num_bath_region, gi)
                        )
                        parf.write(
                            "gregion[{}].g_el    = {}\n".format(j + num_bath_region, ge)
                        )

                        gi, ge = Bspline(
                            0.105 * myo_cv[2], xtrain, Ytrain, b_interp, b_extrap
                        )
                        parf.write(
                            "gregion[{}].g_it    = {}\n".format(j + num_bath_region, gi)
                        )
                        parf.write(
                            "gregion[{}].g_et    = {}\n".format(j + num_bath_region, ge)
                        )

                    else:
                        parf.write(
                            "gregion[{}].num_IDs = 1\n".format(j + num_bath_region)
                        )
                        parf.write(
                            "gregion[{}].ID      = {}\n".format(
                                j + num_bath_region, tag
                            )
                        )

                        gi, ge = Bspline(myo_cv[0], xtrain, Ytrain, b_interp, b_extrap)
                        parf.write(
                            "gregion[{}].g_il    = {}\n".format(j + num_bath_region, gi)
                        )
                        parf.write(
                            "gregion[{}].g_el    = {}\n".format(j + num_bath_region, ge)
                        )

                        gi, ge = Bspline(myo_cv[1], xtrain, Ytrain, b_interp, b_extrap)
                        parf.write(
                            "gregion[{}].g_it    = {}\n".format(j + num_bath_region, gi)
                        )
                        parf.write(
                            "gregion[{}].g_et    = {}\n".format(j + num_bath_region, ge)
                        )

            # ----- EIKONAL REGIONS ----- #
            parf.write(
                "\n# -------------------------------\n# CONDUCTION VELOCITIES\n# -------------------------------\n"
            )

            # number of eikonal regions (including bath)
            parf.write(
                "num_ekregions = {}\n\n".format(
                    num_bath_region + num_myo_with_SE_region
                )
            )

            if any("_theta" in key for key in list(cv_dct.keys())):
                # add bath ekregions
                for j in range(num_bath_region):
                    parf.write("ekregion[{}].ID      = {}\n".format(j, j))
                    parf.write("ekregion[{}].ignore  = 1\n".format(j))

                # "ekregion[*].ID" must be equal to "index" in "gregion[index]" (NOT the element tag).
                parf.write("\n# ----- Ventricles ----- #\n")

                for j, tag in enumerate(vel_f_dct):
                    if not tag in myo_tag_values:
                        parf.write(
                            "ekregion[{}].name  = fast_endocardial_layer\n".format(
                                j + num_bath_region
                            )
                        )
                        parf.write(
                            "ekregion[{}].ID    = {}\n".format(
                                j + num_bath_region, j + num_bath_region
                            )
                        )
                        parf.write(
                            "ekregion[{}].vel_f = {}\n".format(
                                j + num_bath_region, SE_vel_f_dct[tag]
                            )
                        )
                        parf.write(
                            "ekregion[{}].vel_s = {}\n".format(
                                j + num_bath_region, 0.105 * SE_vel_f_dct[tag]
                            )
                        )

                    else:
                        parf.write(
                            "ekregion[{}].ID    = {}\n".format(
                                j + num_bath_region, j + num_bath_region
                            )
                        )
                        parf.write(
                            "ekregion[{}].vel_f = {}\n".format(
                                j + num_bath_region, vel_f_dct[tag]
                            )
                        )
                        parf.write(
                            "ekregion[{}].vel_s = {}\n".format(
                                j + num_bath_region, vel_s_dct[tag]
                            )
                        )
            else:
                # add bath ekregions
                for j in range(num_bath_region):
                    parf.write("ekregion[{}].ID      = {}\n".format(j, j))
                    parf.write("ekregion[{}].ignore  = 1\n".format(j))

                # "ekregion[*].ID" must be equal to "index" in "gregion[index]" (NOT the element tag).
                parf.write("\n# ----- Ventricles ----- #\n")

                for j, tag in enumerate(myo_SE_tag_keys):
                    if not tag in myo_tag_values:
                        parf.write(
                            "ekregion[{}].name  = fast_endocardial_layer\n".format(
                                j + num_bath_region
                            )
                        )
                        parf.write(
                            "ekregion[{}].ID    = {}\n".format(
                                j + num_bath_region, j + num_bath_region
                            )
                        )
                        parf.write(
                            "ekregion[{}].vel_f = {}\n".format(
                                j + num_bath_region, myo_cv[2]
                            )
                        )
                        parf.write(
                            "ekregion[{}].vel_s = {}\n".format(
                                j + num_bath_region, 0.105 * myo_cv[2]
                            )
                        )

                    else:
                        parf.write(
                            "ekregion[{}].ID    = {}\n".format(
                                j + num_bath_region, j + num_bath_region
                            )
                        )
                        parf.write(
                            "ekregion[{}].vel_f = {}\n".format(
                                j + num_bath_region, myo_cv[0]
                            )
                        )
                        parf.write(
                            "ekregion[{}].vel_s = {}\n".format(
                                j + num_bath_region, myo_cv[1]
                            )
                        )
            parf.close()

    return


# --- ALL TISSUE (FIBRE, IMP_REGION, GREGION) --- #
def write_par_tissue(
    CARP_INSTALL_DIR,
    meshfolder,
    meshname,
    labelfilename_im,
    labelfolder,
    labelfilename_scalar,
    samplefolder,
    clusterfolder_samples,
    tag_BATH_jsonfile,
    tune_CV_csvfile,
    im="ToRORd_dynCl",
    im_sv_init="",
    prepacing_lats="",
    prepacing_beats=100,
    prepacing_bcl=1000,
    parfile_tissue="pfile_TISSUE",
    lonfolder_ref="",
    cvfolder_ref="",
    tag_HEART_jsonfile="",
):
    """Write tissue-related .par file, for each sample.

    Parameters
    ----------
    CARP_INSTALL_DIR: str
        Path to CARPentry installation directory.
    meshfolder: str
        Path to the computational mesh.
    meshname: str
        Name of the torso computational mesh.
    labelfile_imp: str
        Name of label file with ionic model parameters.
    labelfolder: str
        Path to the ionic model parameters labelfile.
    labelfilename_scalar: str
        Name of label file with scalar parameters.
    samplefolder: str
        Path to folder with all sample(s) details (usually, local HD).
    clusterfolder_samples: str
        Path to folder will all sample(s) details (usually, HPC or 64-core machines)
    tag_BATH_jsonfile: str
        Path to .json file with bath tag(s) for each bath region.
    tune_CV_csvfile: str
        Path to .csv file with CV, gi and ge training data.
    im: str
        Name of ionic model. Default: "ToRORd_dynCl"
    im_sv_init: str
        Name of initial state file for single cell model (if any). Default: "".
    prepacing_lats: str
        Name of tissue activation file to guide state set-up used by the prepace method (if any). Default: "".
    prepacing_beats: int
        Number of pulses applied to single cell. Default: 100.
    prepacing_bcl: int
        Basic cycle length for repetitive single cell stimulation. Default: 1000 [ms].
    parfile_tissue: str
        Basename of .par file with basic commands to adjust tissue properties (without extension). Default: "pfile_TISSUE".
    lonfolder_ref: str
        Path to .lon file with fibre orientation different than default (without extension), if any. Default: "".
    tag_HEART_jsonfile: str
        Path to .json file with the tags of the heart regions (without extension). Needed when the conduction velocity is not spatially varying. Default: "".
    """

    # find .txt file with accepted parameters
    readmefile = samplefolder + "/accepted_samples"
    if not os.path.exists(readmefile + ".txt"):
        readme = ["nan"] + [str(i) for i in os.listdir(samplefolder)]
    else:
        readme = IO.read_text(readmefile)

    for i in os.listdir(samplefolder):
        if not ".txt" in i and str(i) in readme[1:]:
            print("Initialise EP simulation .par file for " + str(i) + "-th sample...")

            # write par-file to compute lead field matrices
            outputfile = samplefolder + "/" + str(i) + "/" + parfile_tissue
            with open(outputfile + ".par", "w") as file:
                file.write(
                    "# ------------------------------- HETEROGENEITY -------------------------------#\n"
                )
                file.close()

            last_i = i

    # write path to lon file
    if lonfolder_ref:
        for lonfilename in os.listdir(lonfolder_ref):
            if ".lon" in lonfilename:
                write_par_fibre(
                    samplefolder,
                    clusterfolder_samples + "/../ref_QRS",
                    lonfilename,
                    LF_ref=True,
                )
    else:
        for lonfilename in os.listdir(samplefolder + "/" + str(last_i) + "/"):
            if ".lon" in lonfilename:
                write_par_fibre(samplefolder, clusterfolder_samples, lonfilename)
            elif ".blon" in lonfilename:
                write_par_fibre(samplefolder, clusterfolder_samples, lonfilename)

    # write imp section for param-file
    write_par_imp_region(
        CARP_INSTALL_DIR,
        meshfolder,
        meshname,
        labelfilename_im,
        labelfolder,
        labelfilename_scalar,
        samplefolder,
        clusterfolder_samples,
        im=im,
        im_sv_init=im_sv_init,
        prepacing_lats=prepacing_lats,
        prepacing_beats=prepacing_beats,
        prepacing_bcl=prepacing_bcl,
    )

    # write conductivity section of param-file
    write_par_gregion(
        CARP_INSTALL_DIR,
        meshfolder,
        meshname,
        samplefolder,
        tag_BATH_jsonfile,
        tune_CV_csvfile,
        cvfolder_ref=cvfolder_ref,
        tag_HEART_jsonfile=tag_HEART_jsonfile,
    )

    return
