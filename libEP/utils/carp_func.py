#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file is part of 'Cicci_2027_CommsMed'.

__author__ = Ludovica Cicci
__email__ = l.cicci@imperial.ac.uk
"""

# dependencies
from utils import meshtool_func as meshtool
from utils import processing_mesh as proc
from utils import IO

# extra imports
import os
import logging
import shutil
import glob
import numpy as np

logger = logging.getLogger(__name__)


# --- BENCH --- #
def run_single_cell(
    CARP_INSTALL_DIR,
    outputfolder,
    paramfile,
    ionicmodel,
    numstim=1,
    tstart=1.0,
    bcl=1000.0,
    dt_out=1.0,
    savetrace=False,
    savestate=False,
):
    """Single-cell action potential simulation.

    Parameters
    ----------
    CARP_INSTALL_DIR: str
        Path to CARPentry installation directory.
    outputfolder: str
        Path to the directory where bench simulation output will be saved.
    paramfile: str
        Path to text file with string of modified ionic model parameters.
    ionicmodel: str
        Name of the ionic model as given in "bench --list-imp".
    numstim: int
        Number of stimuli. Default: 1.
    tstart: double.
        Start stimulation. Default: 1.
    bcl: double
        Basic cycle length. Default: 1000.
    dt_out: double
        Temporal output granularity. Default: 1.
    savetrace: bool.
        If True, save action potential trace in text file. Default: False.
    savestate: bool.
        If True, save initial state at the end of the pacing protocol. Default: False.
    """

    # find executables
    BENCH = CARP_INSTALL_DIR + "/bin/bench.pt"

    # check folder structure
    print("Checking folder structure...")
    to_check = [paramfile + ".txt"]
    for f in to_check:
        if not os.path.exists(f):
            raise Exception("Cannot find {}. Check path to file.".format(f))
        else:
            print(f + " found.")

    if not os.path.exists(outputfolder):
        os.mkdir(outputfolder)

    APstatsfile = (
        outputfolder
        + "/"
        + ionicmodel
        + "_"
        + os.path.basename(paramfile)
        + "_AP_stats"
    )
    if not os.path.exists(APstatsfile + ".dat") and not os.path.exists(
        APstatsfile + "_AP_stats.dat"
    ):
        # initialise
        cmd = BENCH

        # setup stimulus
        tend = tstart + bcl * numstim

        # configure EP model
        cmd += " --imp=" + ionicmodel
        cmd += (
            " --stim-start "
            + str(tstart)
            + " --numstim "
            + str(numstim)
            + " --bcl "
            + str(bcl)
            + " --duration "
            + str(tend)
        )
        cmd += " --dt-out " + str(dt_out)
        cmd += " --fout=" + APstatsfile

        if savestate:
            # save initial state at the end of the pacing protocol
            cmd += (
                " --save-ini-time "
                + str(tend)
                + " --save-ini-file "
                + outputfolder
                + "/{}_numstim_{}_bcl_{}_dt_{}_baseline.sv".format(
                    ionicmodel, numstim, bcl, dt_out
                )
            )

        # change imp paramers according to param-file
        cmd += " --imp-par " + IO.read_text(paramfile)[0]

        # run bench simulation
        cmd += " --no-trace " + " --APstatistics on"
        print("Executing\n{}".format(cmd))
        os.system(cmd)

        if not savetrace:
            os.remove(APstatsfile + ".txt")

    else:
        print("Single-cell simulation already successful. Skipping bench run.")

    return


# --- UVC --- #
def compute_uvc(
    CARP_INSTALL_DIR,
    meshfolder,
    meshname,
    tags_file,
    uvcfolder,
    _RVmyo="34",
    _LVmyo="35",
    _RVbp="24",
    _LVbp="25",
    _upper_RH="11,12,17,19,22,23",
    _upper_LH="13,14,16,18,20,21",
    cleandir=False,
):
    """Extract biventricular meshes and surfaces from torso mesh, and compute UVC.

    Parameters
    ----------
    CARP_INSTALL_DIR: str
        Path to CARPentry installation directory.
    meshfolder: str
        Path to the computational mesh.
    meshname: str
        Name of the torso computational mesh.
    tags_file: str
        Path to sh file with mesh tags, as required by mguvc function.
    uvcfolder: str
        Path to the directory where universal ventricular coordinates will be saved.
    _RVmyo: str
        Pre-existing tag for RV myocardium. Default: "34".
    _LVmyo: str
        Pre-existing tag for LV myocardium. Default: "35".
    _RVbp: str
        Pre-existing tag for RV blood pool. Default: "24".
    _LVbp: str
        Pre-existing tag for LV blood pool. Default: "25".
    _upper_RH: str
        Pre-existing tag for "upper" right heart, including atrium, valves and vessels (if any). Default: "11,12,17,19,22,23".
    _upper_LH: str
        Pre-existing tag for "upper" left heart, including atrium, valves and vessels (if any). Default: "13,14,16,18,20,21".
    cleandir: bool
        If True, remove simulation output directory. Default: False.
    """

    # find executables
    meshtoolloc = CARP_INSTALL_DIR + "/meshtool/meshtool"

    # IO directories
    if not os.path.exists(uvcfolder):
        os.mkdir(uvcfolder)

    surfloc = uvcfolder + "/surfaces"
    if not os.path.exists(surfloc):
        os.mkdir(surfloc)

    # extract biventricular mesh from torso
    meshpath_i = uvcfolder + "/" + meshname + "_i"
    if not meshtool.extract_mesh(
        meshtoolloc,
        meshfolder + "/" + meshname,
        _RVmyo + "," + _LVmyo,
        meshpath_i,
    ):
        raise Exception("meshtool error. Abort!")

    # extract biventricular base
    if not meshtool.extract_surface(
        meshtoolloc,
        meshfolder + "/" + meshname,
        surfloc + "/" + meshname + "_i.base_left",
        op=_LVmyo + ":" + _upper_LH,
    ):
        raise Exception("meshtool error. Abort!")
    if not meshtool.extract_surface(
        meshtoolloc,
        meshfolder + "/" + meshname,
        surfloc + "/" + meshname + "_i.base_right",
        op=_RVmyo + ":" + _upper_RH,
    ):
        raise Exception("meshtool error. Abort!")
    if not meshtool.extract_surface(
        meshtoolloc,
        meshfolder + "/" + meshname,
        surfloc + "/" + meshname + "_i.base",
        surfloc
        + "/"
        + meshname
        + "_i.base_left"
        + ","
        + surfloc
        + "/"
        + meshname
        + "_i.base_right",
    ):
        raise Exception("meshtool error. Abort!")
    if not meshtool.map_file(
        meshtoolloc,
        meshpath_i,
        surfloc + "/" + meshname + "_i.base.surf",
        uvcfolder,
    ):
        raise Exception("meshtool error. Abort!")
    if not meshtool.map_file(
        meshtoolloc,
        meshpath_i,
        surfloc + "/" + meshname + "_i.base.surf.vtx",
        uvcfolder,
    ):
        raise Exception("meshtool error. Abort!")

    # extract septum surface
    if not meshtool.extract_surface(
        meshtoolloc,
        meshfolder + "/" + meshname,
        surfloc + "/" + meshname + "_i.rvsept",
        op=_RVbp + ":" + _LVmyo,
    ):
        raise Exception("meshtool error. Abort!")
    if not meshtool.map_file(
        meshtoolloc,
        meshpath_i,
        surfloc + "/" + meshname + "_i.rvsept.surf",
        uvcfolder,
    ):
        raise Exception("meshtool error. Abort!")
    if not meshtool.map_file(
        meshtoolloc,
        meshpath_i,
        surfloc + "/" + meshname + "_i.rvsept.surf.vtx",
        uvcfolder,
    ):
        raise Exception("meshtool error. Abort!")

    # extract endo surfaces (LV and RV)
    if not meshtool.extract_surface(
        meshtoolloc,
        meshfolder + "/" + meshname,
        surfloc + "/" + meshname + "_i.lvendo",
        op=_LVmyo + ":" + _LVbp,
    ):
        raise Exception("meshtool error. Abort!")
    if not meshtool.map_file(
        meshtoolloc,
        meshpath_i,
        surfloc + "/" + meshname + "_i.lvendo.surf",
        uvcfolder,
    ):
        raise Exception("meshtool error. Abort!")
    if not meshtool.map_file(
        meshtoolloc,
        meshpath_i,
        surfloc + "/" + meshname + "_i.lvendo.surf.vtx",
        uvcfolder,
    ):
        raise Exception("meshtool error. Abort!")

    if not meshtool.extract_surface(
        meshtoolloc,
        meshfolder + "/" + meshname,
        surfloc + "/" + meshname + "_i.rvendo",
        op=surfloc + "/" + meshname + "_i.rvsept" + "," + _RVmyo + ":" + _RVbp,
    ):
        raise Exception("meshtool error. Abort!")
    if not meshtool.map_file(
        meshtoolloc,
        meshpath_i,
        surfloc + "/" + meshname + "_i.rvendo.surf",
        uvcfolder,
    ):
        raise Exception("meshtool error. Abort!")
    if not meshtool.map_file(
        meshtoolloc,
        meshpath_i,
        surfloc + "/" + meshname + "_i.rvendo.surf.vtx",
        uvcfolder,
    ):
        raise Exception("meshtool error. Abort!")

    #######################################################################

    os.chdir(uvcfolder)

    # generate uvc
    cmd = "mguvc --input-model biv --output-model biv"
    cmd += (
        " --model-name "
        + uvcfolder
        + "/"
        + meshname
        + "_i"
        + " --tags-file "
        + tags_file
    )
    logger.debug(cmd)
    os.system(cmd)

    # rename output dat files to be consistent with the rest of the pipeline
    old_name = ["z", "rho", "phi", "ven"]
    new_name = ["Z", "RHO", "PHI", "biv"]
    for i in range(len(old_name)):
        cmd = (
            "mv "
            + uvcfolder
            + "/output/"
            + meshname
            + "_i.uvc_"
            + old_name[i]
            + ".dat "
            + uvcfolder
            + "/"
            + new_name[i]
            + ".dat"
        )
        logger.debug(cmd)
        os.system(cmd)

    # delete unnecessary files
    if cleandir:
        shutil.rmtree(uvcfolder + "/output")
        shutil.rmtree(uvcfolder + "/surfaces")
        for filename in glob.glob(uvcfolder + "/" + meshname + "_i.*"):
            os.remove(filename)
        os.remove(uvcfolder + "/map_info.txt")

    return


# --- FIBRES --- #
def compute_fibres_system(
    CARP_INSTALL_DIR,
    meshpath_i,
    outputfolder,
    _RVmyo="34",
    _LVmyo="35",
):
    """Compute fibres local reference system.

    Parameters
    ----------
    CARP_INSTALL_DIR: str
        Path to Path to CARPentry installation directory.
    meshpath_i: str
        Path to biventricular mesh.
    outputfolder: str
                Path to the directory where simulation outputs will be saved.
    _RVmyo: str
        Pre-existing tag(s) for RV myocardium.
    _LVmyo: str
        Pre-existing tag(s) for LV myocardium.
    """

    # load biv surfaces
    surf_lvendo = meshpath_i + ".lvendo"
    surf_rvendo = meshpath_i + ".rvendo"
    surf_epi = meshpath_i + ".epi"
    surf_apex = meshpath_i + ".apex"
    surf_base = meshpath_i + ".base"

    # check if outputfolder exists
    if not os.path.exists(outputfolder):
        os.mkdir(outputfolder)

    # define boundaries for each subproblem
    directory = ["EPI", "LV", "RV", "AB"]
    boundary0 = [
        [surf_lvendo, surf_rvendo],
        [surf_epi, surf_rvendo],
        [surf_epi, surf_lvendo],
        [surf_apex],
    ]
    boundary1 = [surf_epi, surf_lvendo, surf_rvendo, surf_base]

    # initialise command
    cmd_init = (
        " +F " + CARP_INSTALL_DIR + "/carputils/carputils/resources/options/pt_ell_amg"
    )
    cmd_init += (
        " +F " + CARP_INSTALL_DIR + "/carputils/carputils/resources/options/pt_para_amg"
    )
    cmd_init += (
        " +F " + CARP_INSTALL_DIR + "/carputils/carputils/resources/options/pt_mech_amg"
    )
    cmd_init += " -ellip_use_pt 1 -parab_use_pt 1 -purk_use_pt 0 -mech_use_pt 1"
    cmd_init += (
        " -purk_options_file "
        + CARP_INSTALL_DIR
        + "/carputils/carputils/resources/petsc_options/petsc_purk_gmres_opts"
    )
    cmd_init += " -mech_finite_element 0 -mech_activate_inertia 0 "

    # create dummy paramfile
    parfile = outputfolder + "/dummy_param.par"
    np.savetxt(parfile, [], header="")

    # solve Laplace problem for epi, LV and RV
    for i in range(3):
        SIMID = outputfolder + "/" + directory[i]
        args = " -experiment 2 -bidomain 1 -timedt 1 -spacedt 1"

        args += " -num_gregions 1 -gregion[0].num_IDs 2"
        args += " -gregion[0].ID[0] " + _RVmyo + " -gregion[0].ID[1] " + _LVmyo
        args += " -gregion[0].g_il 1 -gregion[0].g_it 1 -gregion[0].g_el 1 -gregion[0].g_et 1"

        args += " -num_stim 3"
        args += (
            " -stimulus[0].stimtype 3 -stimulus[0].vtx_file "
            + boundary0[i][0]
            + ".surf"
        )
        args += (
            " -stimulus[1].stimtype 3 -stimulus[1].vtx_file "
            + boundary0[i][1]
            + ".surf"
        )
        args += (
            " -stimulus[2].stimtype 2 -stimulus[2].start 0 -stimulus[2].strength 1 -stimulus[2].duration 1 -stimulus[2].vtx_file "
            + boundary1[i]
            + ".surf"
        )

        cmd = (
            "mpiexec "
            + CARP_INSTALL_DIR
            + "/bin/carp.pt"
            + " +F "
            + parfile
            + " -simID "
            + SIMID
            + " -meshname "
            + meshpath_i
            + args
        )
        cmd += cmd_init
        logger.debug(cmd)
        os.system(cmd)

        # make a datafile out of the igb-file
        cmd = "igbextract -o ascii_1pL -O " + SIMID + "/phie.dat " + SIMID + "/phie.igb"
        logger.debug(cmd)
        os.system(cmd)

    # solve Laplace problem for apex-base
    i = 3
    SIMID = outputfolder + "/" + directory[i]
    args = " -experiment 2 -bidomain 1 -timedt 1 -spacedt 1"

    args += " -num_gregions 1 -gregion[0].num_IDs 2"
    args += " -gregion[0].ID[0] " + _RVmyo + " -gregion[0].ID[1] " + _LVmyo
    args += (
        " -gregion[0].g_il 1 -gregion[0].g_it 1 -gregion[0].g_el 1 -gregion[0].g_et 1"
    )

    args += " -num_stim 2"
    args += " -stimulus[0].stimtype 3 -stimulus[0].vtx_file " + boundary0[i][0]
    args += (
        " -stimulus[1].stimtype 2 -stimulus[1].start 0 -stimulus[1].strength 1 -stimulus[1].duration 1 -stimulus[1].vtx_file "
        + boundary1[i]
        + ".surf"
    )

    cmd = (
        "mpiexec "
        + CARP_INSTALL_DIR
        + "/bin/carp.pt"
        + " +F "
        + parfile
        + " -simID "
        + SIMID
        + " -meshname "
        + meshpath_i
        + args
    )
    cmd += cmd_init
    logger.debug(cmd)
    os.system(cmd)

    # make a datafile out of the igb-file
    cmd = "igbextract -o ascii_1pL -O " + SIMID + "/phie.dat " + SIMID + "/phie.igb"
    logger.debug(cmd)
    os.system(cmd)

    # remove unnecessary files
    os.remove(parfile)
    for currentdir in directory:
        currentpath = outputfolder + "/" + currentdir
        for item in os.listdir(currentpath):
            if item != "phie.dat":
                os.remove(currentpath + "/" + item)
    return






##############################
# PARALLELISATION USING POOL #
##############################
def pool_compute_fibre_preprocessing(
    CARP_INSTALL_DIR,
    meshfolder,
    meshname,  # original torso mesh, without subendocardial layer
    samplefolder,
    labelfolder,
    labelfilename,
    _RVmyo="34",
    _LVmyo="35",
    _RVbp="24",
    _LVbp="25",
    _upper_RH="11,12,17,19,22,23",
    _upper_LH="13,14,16,18,20,21",
):
    """Preprocessing step to compute fibres.

    Parameters
    ----------
    CARP_INSTALL_DIR: str
        Path to Path to CARPentry installation directory.
    meshfolder: str
        Path to the computational mesh.
    meshname: str
        Name of the torso computational mesh. NOTE: Use original torso mesh, without subendocardial layer.
    samplefolder: str
        Path to folder with all sample(s) details.
    labelfolder: str
        Path to the input parameter datafile.
    labelfilename: str
        Name of .txt file with lables of parameters of interest, i.e., GSA inputs, (without extension).
    _RVmyo: str
        Pre-existing tag for RV myocardium. Default: "34".
    _LVmyo: str
        Pre-existing tag for LV myocardium. Default: "35".
    _RVbp: str
        Pre-existing tag for RV blood pool. Default: "24".
    _LVbp: str
        Pre-existing tag for LV blood pool. Default: "25".
    _upper_RH: str
        Pre-existing tag for "upper" right heart, including atrium, valves and vessels (if any). Default: "11,12,17,19,22,23".
    _upper_LH: str
        Pre-existing tag for "upper" left heart, including atrium, valves and vessels (if any). Default: "13,14,16,18,20,21".
    save_bin: bool
        If True, save fibres in binary format. Default: False.
    cleandir: bool
        If True, remove simulation output directory. Default: False.
    """

    # folder with Laplace solutions to compute fibres
    fibrefolder = samplefolder + "/../FIBRE"

    # check if at least one fibre angle is among the input parameters
    par_scl_keys = list(IO.read_text(labelfolder + "/" + labelfilename))
    if not any(
        key in par_scl_keys
        for key in ["alpha_epi", "alpha_endo", "beta_epi", "beta_endo"]
    ):
        print("Fibre angles are not input parameters! Skipping.\n")
        compute_fibres = False
        return compute_fibres, []

    else:
        compute_fibres = True

        print(
            "Fibre angle(s) as input parameter(s). Computing fibre orientation for each sample."
        )

        # pre-processing fibre reference coordinate system
        if not os.path.exists(fibrefolder):
            os.mkdir(fibrefolder)

            # extract biventricular surfaces
            proc.extract_biv_surf(
                CARP_INSTALL_DIR,
                meshfolder,
                meshname,
                fibrefolder + "/surfaces",
                _RVmyo,
                _LVmyo,
                _RVbp,
                _LVbp,
                _upper_RH,
                _upper_LH,
                cleandir=True,
            )

            print("Fibre pre-processing...")

            # compute fibres local reference system to be used with CARPentry RBM
            compute_fibres_system(
                CARP_INSTALL_DIR,
                fibrefolder + "/surfaces/" + meshname + "_i",
                fibrefolder,
                _RVmyo,
                _LVmyo,
            )
        else:
            print("Fibres already pre-processed! Skipping.")

        # find accepted parameters
        readmefile = samplefolder + "/accepted_samples"
        if not os.path.exists(readmefile + ".txt"):
            readme = [str(i) for i in os.listdir(samplefolder) if not ".txt" in i]
        else:
            readme = IO.read_text(readmefile)
            readme = readme[1:]

    return compute_fibres, readme


def pool_compute_fibre_sample(
    id_sample,
    CARP_INSTALL_DIR,
    meshfolder,
    meshname,  # original torso mesh, without subendocardial layer
    samplefolder,
    finalmeshname,
):
    """Compute fibres, for one sample.

    Parameters
    ----------
    id_sample: int
        ID of current sample.
    CARP_INSTALL_DIR: str
        Path to Path to CARPentry installation directory.
    meshfolder: str
        Path to the computational mesh.
    meshname: str
        Name of the torso computational mesh. NOTE: Use original torso mesh, without subendocardial layer.
    samplefolder: str
        Path to folder with all sample(s) details.
    fibrefolder: str
        Path to folder with Laplace solutions to compute fibres.
    finalmeshname: str
        Name of the torso computational mesh with new fibers.
    """

    # folder with Laplace solutions to compute fibres
    fibrefolder = samplefolder + "/../FIBRE"

    # current simulation paramfile
    paramfile = samplefolder + "/" + str(id_sample) + "/param"
    if not os.path.exists(paramfile + ".json"):
        raise Exception("Cannot find {}.".format(paramfile + ".json"))

    # current sample output directory
    SIMID = samplefolder + "/" + str(id_sample)
    if not os.path.exists(SIMID):
        os.mkdir(SIMID)

    # load dictionary with all current simulation parameters (scalar)
    par_all_dct = IO.read_json(paramfile)

    # find values for spatially varying parameter(s)
    angle_dct = {}
    for key in ["alpha_epi", "alpha_endo", "beta_epi", "beta_endo"]:
        angle_dct[key] = par_all_dct[key]

    # compute new fibre orientation file(s)
    if os.path.exists(SIMID + "/" + finalmeshname + ".lon"):
        print("Fibers for sample {} already computed! Skipping.".format(id_sample))
    else:
        rotate_fibre_orientation_mod(
            id_sample,
            CARP_INSTALL_DIR,
            meshfolder,
            meshname,
            fibrefolder + "/surfaces/" + meshname + "_i",
            SIMID,
            finalmeshname,
            fibrefolder,
            angle_dct,
        )

    return


def rotate_fibre_orientation_mod(
    id_sample,
    CARP_INSTALL_DIR,
    meshfolder,
    meshname,
    meshpath_i,
    finalmeshfolder,
    finalmeshname,
    fibrefolder,
    angle_dct,
):
    """Rotate fiber directions.

    Parameters
    ----------
    id_sample: int
        ID of current sample.
    CARP_INSTALL_DIR: str
        Path to CARPentry installation directory.
    meshfolder: str
        Path to the computational mesh.
    meshname: str
        Name of the torso computational mesh.
    meshpath_i: str
        Path to biventricular mesh.
    finalmeshfolder: str
        Path to the directory where the final computational mesh will be saved.
    finalmeshname: str
        Name of the torso computational mesh with new fibers.
    fibrefolder: str
        Path to folder with Laplace solutions to computr fibres.
    angle_dct: dict
        Values for the fibre angles at the epi/endocardium.
    """

    # find executables
    meshtoolloc = CARP_INSTALL_DIR + "/meshtool/meshtool"

    # generate temporary output folder
    outputfolder = finalmeshfolder + "/lon_" + str(id_sample)
    if not os.path.exists(outputfolder):
        os.mkdir(outputfolder)

    # define output meshes
    finalmeshpath = finalmeshfolder + "/" + finalmeshname
    finalmeshpath_i = outputfolder + "/" + finalmeshname + "_i"

    # run fibres generation. NOTE: original lon-file is overwritten
    cmd = (
        "GlRuleFibers -m "
        + meshpath_i
        + " -o "
        + finalmeshpath_i
        + ".lon"
        + " -t biv"
        + " --alpha_endo "
        + str(angle_dct["alpha_endo"])
        + " --alpha_epi "
        + str(angle_dct["alpha_epi"])
        + " --beta_endo "
        + str(angle_dct["beta_endo"])
        + " --beta_epi "
        + str(angle_dct["beta_epi"])
    )
    cmd += " -e " + fibrefolder + "/EPI/phie.dat"
    cmd += " -l " + fibrefolder + "/LV/phie.dat"
    cmd += " -r " + fibrefolder + "/RV/phie.dat"
    cmd += " -a " + fibrefolder + "/AB/phie.dat"
    logger.debug(cmd)
    os.system(cmd)

    # remove sheet from lon file
    lon = IO.read_lon(finalmeshpath_i)
    IO.write_lon(finalmeshpath_i, lon.loc[:, :2])

    # copy file to insert final biventricular mesh back into torso mesh
    shutil.copyfile(meshpath_i + ".elem", finalmeshpath_i + ".elem")
    shutil.copyfile(meshpath_i + ".pts", finalmeshpath_i + ".pts")
    shutil.copyfile(meshpath_i + ".eidx", finalmeshpath_i + ".eidx")
    shutil.copyfile(meshpath_i + ".nod", finalmeshpath_i + ".nod")

    # insert fibers on torso mesh. A new mesh is created
    if not meshtool.insert_submesh(
        meshtoolloc, meshfolder + "/" + meshname, finalmeshpath_i, finalmeshpath
    ):
        raise Exception("meshtool error. Abort!")

    # delete unnecessary files
    shutil.rmtree(outputfolder)
    os.remove(finalmeshpath + ".elem")
    os.remove(finalmeshpath + ".pts")

    return
