#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file is part of 'Cicci_2027_CommsMed'.

__author__ = Ludovica Cicci
__email__ = l.cicci@imperial.ac.uk
"""

# imports
import logging
import os

logger = logging.getLogger(__name__)


def extract_mesh(meshtoolloc, msh, tags, submsh, ifmt="carp_txt", ofmt="carp_txt"):
    cmd = (
        meshtoolloc
        + " extract mesh -msh="
        + msh
        + " -tags="
        + tags
        + " -submsh="
        + submsh
        + " -ifmt="
        + ifmt
        + " -ofmt="
        + ofmt
    )
    logger.debug(cmd)

    return os.system(cmd) == 0


def extract_surface(meshtoolloc, msh, surf, op=None, ifmt="carp_txt", ofmt="carp_txt", *argv):
    cmd = meshtoolloc + " extract surface -msh=" + msh + " -surf=" + surf + " -ifmt=" + ifmt + " -ofmt=" + ofmt
    if not op == None:
        cmd = cmd + " -op=" + op
    if not argv == None:
        for arg in argv:
            cmd = cmd + " " + arg
    logger.debug(cmd)

    return os.system(cmd) == 0


def extract_data(meshtoolloc, submsh, msh_data, submsh_data, mode=0):
    cmd = (
        meshtoolloc
        + " extract data -submsh="
        + submsh
        + " -submsh_data="
        + submsh_data
        + " -msh_data="
        + msh_data
        + " -mode="
        + str(mode)
    )
    logger.debug(cmd)

    return os.system(cmd) == 0


    

def insert_submesh(meshtoolloc, msh, submsh, outmsh):
    cmd = (
        meshtoolloc
        + " insert submesh -msh="
        + msh
        + " -submsh="
        + submsh
        + " -outmsh="
        + outmsh
    )

    logger.debug(cmd)
    return os.system(cmd) == 0


def insert_data(meshtoolloc, submsh, submsh_data, msh, odat, mode):
    cmd = (
        meshtoolloc
        + " insert data -submsh="
        + submsh
        + " -submsh_data="
        + submsh_data
        + " -msh="
        + msh
        + " -odat="
        + odat
        + " -mode="
        + mode
    )
    logger.debug(cmd)
    return os.system(cmd) == 0


def map_file(meshtoolloc, submsh, files, outdir, mode="m2s"):
    cmd = (
        meshtoolloc
        + " map -submsh="
        + submsh
        + " -files="
        + files
        + "  -outdir="
        + outdir
        + " -mode="
        + mode
    )

    logger.debug(cmd)
    return os.system(cmd) == 0


def generate_distancefield(meshtoolloc, msh, ssurf, odat, *argv):
    cmd = (
        meshtoolloc
        + " generate distancefield -msh="
        + msh
        + " -ssurf="
        + ssurf
        + " -odat="
        + odat
        + " "
    )
    if not argv == None:
        for arg in argv:
            cmd = cmd + " " + arg

    logger.debug(cmd)
    return os.system(cmd) == 0




def query_tags(meshtoolloc, msh, outputfile):
    cmd = meshtoolloc + " query tags -msh=" + msh
    logger.debug(cmd)

    pipe = os.popen(cmd)
    with open(outputfile + ".txt", "w") as file:
        file.write(pipe.read())
        file.close()

    return pipe.close() == None


