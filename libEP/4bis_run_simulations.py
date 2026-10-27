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
from natsort import natsorted
from pytictoc import TicToc

# make config.py (in the parent folder) importable
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# dependencies
from config import CARP_INSTALL_DIR, OUTPUT_DIR, check_paths
from utils import IO


def main():
    check_paths("CARP_INSTALL_DIR", "OUTPUT_DIR")

    meshname = "KCL1" # NOTE: mesh basename

    # folder with "mesh", "data" and the samples (same as clusterfolder in 4_simulation_assembling.py)
    clusterfolder = OUTPUT_DIR + "/" + meshname

    # name of mesh to be used in EP simulation(s)
    meshname_EP = meshname + "_gregion"

    # folder with the samples (same as samplefolder in 4_simulation_assembling.py)
    samplefolder = clusterfolder + "/samples_LHS"

    # output folder
    simulationfolder = clusterfolder + "/simulations_LHS"
    os.makedirs(simulationfolder, exist_ok=True)

    # parameter file with the EP model parameters
    parfile_base = clusterfolder + "/data/pfile_BASE.par"
    parfile_lf = clusterfolder + "/data/pfile_LF.par"
    parfile_ecg = clusterfolder + "/data/pfile_ECG.par"

    # number of processors
    num_proc = 64

    # ---> START

    # find executables
    CARP = CARP_INSTALL_DIR + "/bin/carp.pt"
    Parab = CARP_INSTALL_DIR + "/carputils/carputils/resources/petsc_options/ilu_cg_opts"
    OPTS = CARP_INSTALL_DIR + "/carputils/carputils/resources/petsc_options/gamg_cg_opts"

    # find .txt file with accepted parameters
    readmefile = samplefolder + "/accepted_samples"
    if not os.path.exists(readmefile + ".txt"):
        readme = ["nan"] + [str(i) for i in os.listdir(samplefolder) if not ".txt" in str(i)]
    else:
        readme = IO.read_text(readmefile)

    # loop on samples
    for i in natsorted(os.listdir(samplefolder)):
        if not ".txt" in i and str(i) in readme[1:]:

            SIMID = samplefolder + "/" + str(i)

            # parameter file with the EP model parameters
            parfile_ep = SIMID + "/pfile_EP.par"
            parfile_tissue = SIMID + "/pfile_TISSUE.par"

            time = TicToc()
            time.tic()

            # path to the folder with the lead field arrays (if used, otherwise set equal to None)
            lffolder = simulationfolder + "/SIM_LF_" + str(i)
            if not os.path.exists(lffolder):
                os.mkdir(lffolder)

                # compute lead field arrays
                cmd = (
                    "mpiexec -np "
                    + str(num_proc)
                    + " "
                    + CARP
                    + " +F "
                    + parfile_base
                    + " +F "
                    + parfile_lf
                    + " +F "
                    + parfile_tissue
                    + " -meshname "
                    + clusterfolder
                    + "/mesh/"
                    + meshname_EP
                    + " -simID "
                    + lffolder
                    + " -parab_options_file "
                    + Parab
                    + " -ellip_options_file "
                    + OPTS
                )
                print(cmd)
                os.system(cmd)
            else:
                print(
                    "WARNING: {} already exists. Skipping the lead field simulation. If it failed, delete the folder to run it again.".format(
                        lffolder
                    )
                )

            # directories where the simulation (i.e., voltages) results are stored
            epfolder = simulationfolder + "/SIM_EP_" + str(i)
            if not os.path.exists(epfolder):
                os.mkdir(epfolder)

                # run EP simulation
                cmd = (
                    "mpiexec -np "
                    + str(num_proc)
                    + " "
                    + CARP
                    + " +F "
                    + parfile_base
                    + " +F "
                    + parfile_ecg
                    + " +F "
                    + parfile_tissue
                    + " +F "
                    + parfile_ep
                    + " -meshname "
                    + clusterfolder
                    + "/mesh/"
                    + meshname_EP
                    + " -simID "
                    + epfolder
                    + " -lead_field_dir "
                    + lffolder
                    + " -parab_options_file "
                    + Parab
                    + " -ellip_options_file "
                    + OPTS
                )
                print(cmd)
                os.system(cmd)
            else:
                print(
                    "WARNING: {} already exists. Skipping the EP simulation. If it failed, delete the folder to run it again.".format(
                        epfolder
                    )
                )

            time.toc()
    # <--- END


if __name__ == "__main__":
    main()
