#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file is part of 'Cicci_2027_CommsMed'.

__author__ = Ludovica Cicci
__email__ = l.cicci@imperial.ac.uk
"""

import os

# Edit the two paths below before running any script.
CARP_INSTALL_DIR = ""  # folder where CARPentry is installed
OUTPUT_DIR = ""  # folder where results are written

# Folder that contains this file, i.e. the root of 'Cicci_2027_CommsMed'
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
libEP_DIR = os.path.join(ROOT_DIR, "libEP")


def check_paths(*names):
    for name in names:
        if not globals()[name]:
            raise SystemExit(
                "Please set {} in {}".format(name, os.path.join(ROOT_DIR, "config.py"))
            )
