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
import copy


# --- HEART TAGS --- #
def fec_layers(Z_FEC_min=0.50, Z_FEC_max=0.75, deltaZ_FEC=0.02):
    """Apico-basal extension (in %) at which each fast endocardial conductive layer starts, as tagged in compute_FEC_layers.

    Parameters
    ----------
    Z_FEC_min: float
        Lower bound for fast endocardial conductive layer extension. Must be a multiple of 0.01. Default: 0.50.
    Z_FEC_max: float
        Upper bound for fast endocardial conductive layer extension. Must be a multiple of 0.01. Default: 0.75.
    deltaZ_FEC: float
        Extension of each fast endocardial conductive layer. Must be a multiple of 0.01. Default: 0.02.
    """

    values = [100 * Z_FEC_min, 100 * Z_FEC_max, 100 * deltaZ_FEC]
    if any(abs(v - round(v)) > 1e-6 for v in values):
        raise Exception(
            "Z_FEC_min, Z_FEC_max and deltaZ_FEC must be multiples of 0.01 (given: {}, {}, {}).".format(
                Z_FEC_min, Z_FEC_max, deltaZ_FEC
            )
        )
    if round(values[2]) < 1:
        raise Exception("deltaZ_FEC must be at least 0.01.")

    return list(range(round(values[0]), round(values[1]) + 1, round(values[2])))


def _derive_tag_sets(base, par, layers):
    """Tags of the segments and of the fast endocardial conductive layers, derived from the base tags and the construction parameters (see segment_submesh and compute_FEC_layers)."""

    factor_myo = par["factor_tag_myo"]
    factor_SE = par["factor_tag_SE"]
    factor_FEC = par["factor_tag_FEC"]
    Nz = par["Nz"]

    sets = {"myo_segments": {}, "SE_segments": {}, "FEC_tags": {}, "FEC_segments": {}}
    for name, Nrho in (("RVmyo", par["Nrho_RV"]), ("LVmyo", par["Nrho_LV"])):
        sets["myo_segments"][name] = [
            factor_myo * int(base[name]) + k for k in range(1, Nz * Nrho + 1)
        ]
    for name in ("SE_RV", "SE_LV", "SE_septRV"):
        b = int(base[name])
        sets["SE_segments"][name] = [factor_SE * b + z for z in range(1, Nz + 1)]
        sets["FEC_tags"][name] = [factor_FEC * b + layer for layer in layers]
        sets["FEC_segments"][name] = [
            factor_SE * t + z for t in sets["FEC_tags"][name] for z in range(1, Nz + 1)
        ]

    return sets


def get_heart_tags(tagfile, tagbathfile=""):
    """Read the mesh tags and the construction parameters of the regions from a .json file, derive the tags of all the regions, and check that they do not clash.

    Parameters
    ----------
    tagfile: str
        Path to .json file with the tags of the input mesh ("base") and the parameters used to segment the biventricle ("construction"), without extension.
    tagbathfile: str
        Path to .json file with the tags of the bath regions, without extension. If given, the tags of the regions are checked against them too. Default: "".

    Returns
    -------
    tags: dict
        Comma-separated tags (as strings) with keys
        - "RVmyo", "LVmyo", "SE_RV", "SE_LV", "SE_septRV", "RVbp", "LVbp", "upper_RH", "upper_LH": tags in the input mesh;
        - "RVmyo_subendolayer", "LVmyo_subendolayer", "SE_RV_subendolayer", "SE_LV_subendolayer", "SE_septRV_subendolayer": tags in the sub-endocardial layer mesh (<meshname>_subendolayer);
        - "RVmyo_FEC", "LVmyo_FEC": myocardium tags in the mesh with the fast endocardial conductive layers (<meshname>_FEC);
        - "SE_RV_FEC", "SE_LV_FEC", "SE_septRV_FEC": sub-endocardial layer tag, plus the tags of the fast endocardial conductive layers (as given by compute_FEC_layers; tags in the mesh <meshname>_FEC, input of segment_biventricle);
        - "RVmyo_gregion", "LVmyo_gregion": myocardium tag, plus the tags of its segments in the segmented mesh;
        - "SE_RV_gregion", "SE_LV_gregion", "SE_septRV_gregion": sub-endocardial layer tag, plus the tags of its segments and the tags of the segments of the fast endocardial conductive layers in the segmented mesh.
    par: dict
        Construction parameters, as in the .json file.
    """

    data = IO.read_json(tagfile)
    base = data["base"]
    par = data["construction"]

    # factors and numbers of regions must be integers (otherwise tags would be written as e.g. 3401.0). Floats with an integer value (e.g., 1e2, 100.0) are converted
    for key in ("factor_tag_myo", "factor_tag_SE", "factor_tag_FEC", "Nz", "Nrho_LV", "Nrho_RV"):
        value = par[key]
        if isinstance(value, bool) or not isinstance(value, (int, float)) or value != int(value):
            raise Exception(
                "'{}' in {} must be an integer (given: {}).".format(key, tagfile, value)
            )
        par[key] = int(value)

    factor_myo = par["factor_tag_myo"]
    factor_SE = par["factor_tag_SE"]
    factor_FEC = par["factor_tag_FEC"]
    Nz = par["Nz"]
    layers = fec_layers(par["Z_FEC_min"], par["Z_FEC_max"], par["deltaZ_FEC"])

    # setup.write_par_gregion reads the layer number from the digits 3 and 4 of the tag: this works only with the factor 100 and layers with two digits
    if (factor_myo, factor_SE, factor_FEC) != (100, 100, 100) or layers[0] < 10 or layers[-1] > 99:
        raise Exception(
            "Tags in {} are not supported: factor_tag_myo, factor_tag_SE and factor_tag_FEC must be 100 and the fast layers (100 * Z_FEC_min to 100 * Z_FEC_max) must have two digits (10 to 99), "
            "because setup.write_par_gregion reads the layer number from the digits 3 and 4 of the tag (given: factors {}, {}, {}; layers {} to {}).".format(
                tagfile, factor_myo, factor_SE, factor_FEC, layers[0], layers[-1]
            )
        )

    # tags are integers
    for key in base:
        for t in base[key].split(","):
            if not t.strip().isdigit():
                raise Exception("Tag '{}' of '{}' in {} is not a positive integer.".format(t, key, tagfile))

    # tags of the segments: see segment_submesh
    sets = _derive_tag_sets(base, par, layers)

    def to_str(lst):
        return ",".join([str(t) for t in lst])

    tags = dict(base)
    families = {}
    for name in ("RVmyo", "LVmyo"):
        segments = sets["myo_segments"][name]
        tags[name + "_subendolayer"] = base[name]
        tags[name + "_FEC"] = base[name]
        tags[name + "_gregion"] = to_str([base[name]] + segments)
        families[name] = [int(base[name])] + segments
    for name in ("SE_RV", "SE_LV", "SE_septRV"):
        b = base[name]
        tags[name + "_subendolayer"] = b
        tags[name + "_FEC"] = to_str([b] + sets["FEC_tags"][name])
        tags[name + "_gregion"] = to_str(
            [b] + sets["SE_segments"][name] + sets["FEC_segments"][name]
        )
        families[name] = (
            [int(b)]
            + sets["SE_segments"][name]
            + sets["FEC_tags"][name]
            + sets["FEC_segments"][name]
        )

    # checks
    if par["Nrho_RV"] > 5 or par["Nrho_LV"] > 5:
        raise Exception("Nrho_RV and Nrho_LV must be at most 5.")
    if par["Nz"] > 10:
        raise Exception("Nz must be at most 10.")
    problems = []
    if Nz * max(par["Nrho_RV"], par["Nrho_LV"]) >= factor_myo:
        problems.append(
            "Nz * Nrho ({}) must be smaller than factor_tag_myo ({}).".format(
                Nz * max(par["Nrho_RV"], par["Nrho_LV"]), factor_myo
            )
        )
    if Nz >= factor_SE:
        problems.append("Nz ({}) must be smaller than factor_tag_SE ({}).".format(Nz, factor_SE))
    if layers[-1] >= factor_FEC:
        problems.append(
            "The last fast layer ({}) must be smaller than factor_tag_FEC ({}).".format(layers[-1], factor_FEC)
        )

    # tags of the regions must be different from each other, and from the tags of the other regions of the input mesh
    owner = {}
    for name, lst in families.items():
        for t in lst:
            owner.setdefault(t, []).append(name)
    others = {name: base[name] for name in ("RVbp", "LVbp", "upper_RH", "upper_LH")}
    if tagbathfile:
        bath = IO.read_json(tagbathfile)
        for name in bath:
            others["bath/" + name] = ",".join(bath[name])
    for name, lst in others.items():
        for t in lst.split(","):
            if t.strip().isdigit() and int(t) in owner:
                owner[int(t)].append(name)
    for t in sorted(owner):
        if len(owner[t]) > 1:
            problems.append("Tag {} is used by {}.".format(t, " and ".join(owner[t])))
    if problems:
        raise Exception(
            "Tags in {} are not valid:\n".format(tagfile) + "\n".join(problems[:20])
        )

    return tags, par


def get_myocardium_tags(tags_in_mesh, Zmax_SE, tagfile):
    """Tags (among the tags of a mesh) of ordinary myocardium, i.e., not of the fast endocardial conductive layer, for a given extension of the fast layer.

    Myocardium tags (as in the input mesh, and their segments) are myocardium. Sub-endocardial layer tags (and their segments) are fast layer.
    Tags of the fast endocardial conductive layers (and their segments) are fast layer if their layer (in % of the apico-basal extension) is not above 100 * Zmax_SE, otherwise they are myocardium.

    Parameters
    ----------
    tags_in_mesh: list
        Tags (integers) present in the mesh.
    Zmax_SE: float
        Extension in the apico-basal direction of the fast subendocardial layer. Must be in (0,1].
    tagfile: str
        Path to .json file with the tags of the input mesh ("base") and the parameters used to segment the biventricle ("construction"), without extension.
    """

    tags, par = get_heart_tags(tagfile)
    layers = fec_layers(par["Z_FEC_min"], par["Z_FEC_max"], par["deltaZ_FEC"])
    sets = _derive_tag_sets(tags, par, layers)

    myo = set()
    for name in ("RVmyo", "LVmyo"):
        myo.update([int(tags[name])] + sets["myo_segments"][name])
    SE = set()
    layer_of = {}
    for name in ("SE_RV", "SE_LV", "SE_septRV"):
        SE.update([int(tags[name])] + sets["SE_segments"][name])
        for layer, t in zip(layers, sets["FEC_tags"][name]):
            layer_of[t] = layer
            for z in range(1, par["Nz"] + 1):
                layer_of[par["factor_tag_SE"] * t + z] = layer

    unknown = [t for t in tags_in_mesh if t not in myo and t not in SE and t not in layer_of]
    if unknown:
        raise Exception(
            "Tags {} of the mesh are not in {}. Update the file.".format(unknown[:20], tagfile)
        )

    return [t for t in tags_in_mesh if t in myo or (t in layer_of and layer_of[t] > 100 * Zmax_SE)]


def count_heart_regions(meshfolder, meshname, tags, max_regions=256):
    """Print the number of different myocardial region tags of a segmented mesh, and warn if it is larger than the maximum number of regions allowed by CARPentry.

    Parameters
    ----------
    meshfolder: str
        Path to the computational mesh.
    meshname: str
        Name of the segmented mesh.
    tags: dict
        Tags of the regions, as returned by get_heart_tags.
    max_regions: int
        Maximum number of regions allowed by CARPentry. Default: 256.
    """

    elem = IO.read_elem(meshfolder + "/" + meshname)
    heart = set()
    for name in ("RVmyo", "LVmyo", "SE_RV", "SE_LV", "SE_septRV"):
        heart.update([int(t) for t in tags[name + "_gregion"].split(",")])
    n_regions = len(heart.intersection(set(elem[5].unique())))
    print("Mesh {} has {} different myocardial region tags.".format(meshname, n_regions))
    if n_regions > max_regions:
        print(
            "WARNING: CARPentry does not allow more than {} different regions.".format(max_regions)
        )

    return n_regions


def extract_biv_surf(
    CARP_INSTALL_DIR,
    meshfolder,
    meshname,
    outputfolder,
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
    outputfolder: str
        Path to the directory where universal ventricular coordinates will be saved.
    _RVmyo: str
        Pre-existing tag for RV myocardium. Default: "34".
    _LVmyo: str
        Pre-existing tag for LV myocardium. Default: "35".
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
    if not os.path.exists(outputfolder):
        os.mkdir(outputfolder)

    mesh_iloc = outputfolder + "/volumes"
    if not os.path.exists(mesh_iloc):
        os.mkdir(mesh_iloc)

    surfloc = outputfolder + "/surfaces"
    if not os.path.exists(surfloc):
        os.mkdir(surfloc)

    apexloc = outputfolder + "/apex"
    if not os.path.exists(apexloc):
        os.mkdir(apexloc)

    # extract biventricular mesh from torso
    meshpath_i = outputfolder + "/" + meshname + "_i"
    if not meshtool.extract_mesh(
        meshtoolloc,
        meshfolder + "/" + meshname,
        _RVmyo + "," + _LVmyo,
        meshpath_i,
    ):
        raise Exception("meshtool error. Abort!")

    # extract LV from biv
    meshpath_i_LV = mesh_iloc + "/" + meshname + "_i_LV"
    if not meshtool.extract_mesh(meshtoolloc, meshpath_i, _LVmyo, meshpath_i_LV):
        raise Exception("meshtool error. Abort!")
    if not meshtool.extract_surface(meshtoolloc, meshpath_i, meshpath_i_LV, op=_LVmyo):
        raise Exception("meshtool error. Abort!")

    # extract RV from biv
    meshpath_i_RV = mesh_iloc + "/" + meshname + "_i_RV"
    if not meshtool.extract_mesh(
        meshtoolloc,
        meshpath_i,
        _RVmyo,
        meshpath_i_RV,
    ):
        raise Exception("meshtool error. Abort!")
    if not meshtool.extract_surface(
        meshtoolloc,
        meshpath_i,
        meshpath_i_RV,
        op=_RVmyo,
    ):
        raise Exception("meshtool error. Abort!")

    # extract septum surface
    meshpath_i_rvsept = surfloc + "/" + meshname + "_i.rvsept"
    if not meshtool.extract_surface(
        meshtoolloc,
        meshfolder + "/" + meshname,
        meshpath_i_rvsept,
        op=_RVbp + ":" + _LVmyo,
    ):
        raise Exception("meshtool error. Abort!")
    if not meshtool.map_file(
        meshtoolloc,
        meshpath_i,
        meshpath_i_rvsept + ".surf.vtx",
        outputfolder,
    ):
        raise Exception("meshtool error. Abort!")

    # extract endo surfaces (LV and RV)
    meshpath_i_lvendo = surfloc + "/" + meshname + "_i.lvendo"
    if not meshtool.extract_surface(
        meshtoolloc,
        meshfolder + "/" + meshname,
        meshpath_i_lvendo,
        op=_LVmyo + ":" + _LVbp,
    ):
        raise Exception("meshtool error. Abort!")
    if not meshtool.map_file(
        meshtoolloc,
        meshpath_i,
        meshpath_i_lvendo + ".surf.vtx",
        outputfolder,
    ):
        raise Exception("meshtool error. Abort!")

    meshpath_i_rvendo = surfloc + "/" + meshname + "_i.rvendo"
    if not meshtool.extract_surface(
        meshtoolloc,
        meshfolder + "/" + meshname,
        meshpath_i_rvendo,
        op=meshpath_i_rvsept + "," + _RVmyo + ":" + _RVbp,
    ):
        raise Exception("meshtool error. Abort!")
    if not meshtool.map_file(
        meshtoolloc,
        meshpath_i,
        meshpath_i_rvendo + ".surf.vtx",
        outputfolder,
    ):
        raise Exception("meshtool error. Abort!")

    # extract epi surface
    meshpath_i_lvepi = surfloc + "/" + meshname + "_i.lvepi"
    if not meshtool.extract_surface(
        meshtoolloc,
        meshfolder + "/" + meshname,
        meshpath_i_lvepi,
        op=_LVmyo + "-" + _LVbp,
    ):
        raise Exception("meshtool error. Abort!")
    if not meshtool.map_file(
        meshtoolloc,
        meshpath_i,
        meshpath_i_lvepi + ".surf.vtx",
        outputfolder,
    ):
        raise Exception("meshtool error. Abort!")

    meshpath_i_rvepi = surfloc + "/" + meshname + "_i.rvepi"
    if not meshtool.extract_surface(
        meshtoolloc,
        meshfolder + "/" + meshname,
        meshpath_i_rvepi,
        op=_RVmyo + "-" + _RVbp,
    ):
        raise Exception("meshtool error. Abort!")
    if not meshtool.map_file(
        meshtoolloc,
        meshpath_i,
        meshpath_i_rvepi + ".surf.vtx",
        outputfolder,
    ):
        raise Exception("meshtool error. Abort!")

    meshpath_i_epi = surfloc + "/" + meshname + "_i.epi"
    if not meshtool.extract_surface(
        meshtoolloc,
        meshfolder + "/" + meshname,
        meshpath_i_epi,
        op=meshpath_i_lvepi + "," + meshpath_i_rvepi,
    ):
        raise Exception("meshtool error. Abort!")
    if not meshtool.map_file(
        meshtoolloc,
        meshpath_i,
        meshpath_i_epi + ".surf.vtx",
        outputfolder,
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
        surfloc + "/" + meshname + "_i.base.surf.vtx",
        outputfolder,
    ):
        raise Exception("meshtool error. Abort!")

    # extract biventricular apex

    # extract upper LH
    meshpath_LH = apexloc + "/LH"
    if not meshtool.extract_mesh(
        meshtoolloc, meshfolder + "/" + meshname, _LVmyo + "," + _upper_LH, meshpath_LH
    ):
        raise Exception("meshtool error. Abort!")

    # extract upper LH and compute distance from LV (in upperLH)
    meshpath_upperLH = apexloc + "/upperLH"
    if not meshtool.extract_surface(
        meshtoolloc, meshpath_LH, meshpath_upperLH, op=_upper_LH
    ):
        raise Exception("meshtool error. Abort!")
    if not meshtool.generate_distancefield(
        meshtoolloc,
        meshpath_LH,
        meshpath_upperLH,
        apexloc + "/dist_LHlayerfromLV.dat",
    ):
        raise Exception("meshtool error. Abort!")

    # find apex as furthest point in LV from everything else in LH
    dist_from_LV = np.genfromtxt(apexloc + "/dist_LHlayerfromLV.dat")
    idx_apex_LH = np.where(dist_from_LV == np.amax(dist_from_LV))[0]
    pts_LH = np.array(IO.read_pts(meshpath_LH))

    # find the apex point' index in biv mesh
    pts_biv = np.array(IO.read_pts(meshpath_i))
    idx_apex_biv = np.where(np.all(pts_biv == pts_LH[idx_apex_LH], axis=1))[0]
    IO.write_vtx(
        outputfolder + "/" + meshname + "_i.apex",
        pd.DataFrame(idx_apex_biv),
    )

    # delete unnecessary files
    if cleandir:
        shutil.rmtree(mesh_iloc)
        shutil.rmtree(surfloc)
        shutil.rmtree(apexloc)

    return


def extract_SElayer(
    CARP_INSTALL_DIR,
    meshfolder,
    meshname,
    finalmeshfolder,
    SElayermeshname,
    finalmeshname,
    uvcfolder,
    _RVmyo="34",
    _LVmyo="35",
    _RVbp="24",
    _LVbp="25",
    _upper_RH="11,12,17,19,22,23",
    _upper_LH="13,14,16,18,20,21",
    _SE_RV="44",
    _SE_LV="45",
    _SE_septRV="46",
    Zmax_SE=0.75,
    cleandir=False,
):
    """Compute subendocardial layer (from apex to base, 1 element thick in the transmural direction) and fast conducting endocardial layer (from apex to Zmax_SE).

    Parameters
    ----------
    CARP_INSTALL_DIR: str
        Path to CARPentry installation directory.
    meshfolder: str
        Path to the computational mesh.
    meshname: str
        Name of the torso computational mesh.
    finalmeshfolder: str
        Path to the directory where the final computational mesh will be saved.
    SElayermeshname: str
        Name of the torso computational mesh with new tags for the subendocardial layer.
    finalmeshname: str
        Name of the torso computational mesh with new tags for the fast conducting endocardial layer.
    uvcfolder: str
        Path to the directory where universal ventricular coordinates are saved.
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
    _SE_RV: str
        Tag for RV subendocardial layer.
    _SE_LV: str
        Tag for LV subendocardial layer.
    _SE_septRV: str
        Tag for septal subendocardial layer on the right side.
    Zmax_SE: float
        Extension in the apico-basal direction of the fast subendocardial layer. Must be in (0,1]. Default: 0.75.
    cleandir: bool
        If True, remove simulation output directory. Default: False.
    """

    # find executables
    meshtoolloc = CARP_INSTALL_DIR + "/meshtool/meshtool"

    # check folder structure
    if not os.path.exists(uvcfolder):
        raise Exception(
            "UVC folder does not exist. Must be same folder where the universal ventricular coordinates are saved."
        )

    # check value for Zmax_SE
    if Zmax_SE < 0 or Zmax_SE > 1:
        raise Exception("Maximal extension for subendocardial layer must be in [0,1].")

    # IO directories
    simresult = finalmeshfolder + "/SElayer"
    if not os.path.exists(simresult):
        os.mkdir(simresult)

    # volumetric meshes
    meshloc = simresult + "/volumes"
    if not os.path.exists(meshloc):
        os.mkdir(meshloc)

    # volumetric meshes extracted from biv
    mesh_iloc = meshloc + "/biv"
    if not os.path.exists(mesh_iloc):
        os.mkdir(mesh_iloc)

    # surfaces extracted from torso
    surfloc = simresult + "/surfaces"
    if not os.path.exists(surfloc):
        os.mkdir(surfloc)

    # surfaces extracted from biv
    surf_iloc = surfloc + "/biv"
    if not os.path.exists(surf_iloc):
        os.mkdir(surf_iloc)

    # surfaces extracted from biv-LV
    surf_left_loc = surf_iloc + "/LV"
    if not os.path.exists(surf_left_loc):
        os.mkdir(surf_left_loc)

    # surfaces extracted from biv-RV
    surf_right_loc = surf_iloc + "/RV"
    if not os.path.exists(surf_right_loc):
        os.mkdir(surf_right_loc)

    # surfaces extracted from LV/RV
    surf_iloc_layer = surf_iloc + "/layer"
    if not os.path.exists(surf_iloc_layer):
        os.mkdir(surf_iloc_layer)

    apexloc = simresult + "/apex"
    if not os.path.exists(apexloc):
        os.mkdir(apexloc)

    SIMID = simresult + "/layer"
    if not os.path.exists(SIMID):
        os.mkdir(SIMID)

    # extract biv from torso
    meshpath_i = meshloc + "/" + meshname + "_i"
    if not meshtool.extract_mesh(
        meshtoolloc, meshfolder + "/" + meshname, _RVmyo + "," + _LVmyo, meshpath_i
    ):
        raise Exception("meshtool error. Abort!")

    # extract LV from biv
    meshpath_i_LV = mesh_iloc + "/" + meshname + "_i_LV"
    if not meshtool.extract_mesh(meshtoolloc, meshpath_i, _LVmyo, meshpath_i_LV):
        raise Exception("meshtool error. Abort!")

    # extract RV from biv
    meshpath_i_RV = mesh_iloc + "/" + meshname + "_i_RV"
    if not meshtool.extract_mesh(meshtoolloc, meshpath_i, _RVmyo, meshpath_i_RV):
        raise Exception("meshtool error. Abort!")

    # extract RV septum surface from torso and map it on biv
    septumsurf = meshname + "_i.rvsept"
    if not meshtool.extract_surface(
        meshtoolloc,
        meshfolder + "/" + meshname,
        surfloc + "/" + septumsurf,
        op=_RVbp + ":" + _LVmyo,
    ):
        raise Exception("meshtool error. Abort!")
    if not meshtool.map_file(
        meshtoolloc, meshpath_i, surfloc + "/" + septumsurf + ".surf", surf_iloc
    ):
        raise Exception("meshtool error. Abort!")

    # extract LV endo surface from torso and map it on biv and on LV
    endosurf_LV = meshname + "_i.lvendo"
    if not meshtool.extract_surface(
        meshtoolloc,
        meshfolder + "/" + meshname,
        surfloc + "/" + endosurf_LV,
        op=_LVmyo + ":" + _LVbp,
    ):
        raise Exception("meshtool error. Abort!")
    if not meshtool.map_file(
        meshtoolloc, meshpath_i, surfloc + "/" + endosurf_LV + ".surf", surf_iloc
    ):
        raise Exception("meshtool error. Abort!")
    if not meshtool.map_file(
        meshtoolloc,
        meshpath_i_LV,
        surf_iloc + "/" + endosurf_LV + ".surf",
        surf_left_loc,
    ):
        raise Exception("meshtool error. Abort!")

    # extract RV endo surface from torso and map it on biv and on RV
    endosurf_RV = meshname + "_i.rvendo"
    if not meshtool.extract_surface(
        meshtoolloc,
        meshfolder + "/" + meshname,
        surfloc + "/" + endosurf_RV,
        op=surfloc + "/" + septumsurf + "," + _RVmyo + ":" + _RVbp,
    ):
        raise Exception("meshtool error. Abort!")
    if not meshtool.map_file(
        meshtoolloc, meshpath_i, surfloc + "/" + endosurf_RV + ".surf", surf_iloc
    ):
        raise Exception("meshtool error. Abort!")
    if not meshtool.map_file(
        meshtoolloc,
        meshpath_i_RV,
        surf_iloc + "/" + endosurf_RV + ".surf",
        surf_right_loc,
    ):
        raise Exception("meshtool error. Abort!")

    # ==== SElayer FULL EXTENSION ==== #
    # set name and path for mesh with SElayer (extented from apex to base)
    SElayermeshpath = finalmeshfolder + "/" + SElayermeshname

    # extract biv from torso (new mesh to be overwritten)
    SElayermeshpath_i = meshloc + "/" + SElayermeshname + "_i"
    if not meshtool.extract_mesh(
        meshtoolloc,
        meshfolder + "/" + meshname,
        _RVmyo + "," + _LVmyo,
        SElayermeshpath_i,
    ):
        raise Exception("meshtool error. Abort!")
    elem_SE = IO.read_elem(SElayermeshpath_i)

    # extract LV subendocardial layer (1 element thick)
    index_lvendo = IO.read_surf_to_nodeList(surf_iloc + "/" + endosurf_LV)
    SElayer_LV = elem_SE[[1, 2, 3, 4]].map(lambda x: x in index_lvendo)
    SElayer_LV_idx = np.where(SElayer_LV.any(axis=1).map(int) == 1)[0]
    elem_SE.loc[SElayer_LV_idx, 5] = int(_SE_LV) + (
        int(_LVmyo) - elem_SE.loc[SElayer_LV_idx, 5]
    )

    # extract RV subendocardial layer (1 element thick)
    index_rvendo = IO.read_surf_to_nodeList(surf_iloc + "/" + endosurf_RV)
    SElayer_RV = elem_SE[[1, 2, 3, 4]].map(lambda x: x in index_rvendo)
    SElayer_RV_idx = np.where(SElayer_RV.any(axis=1).map(int) == 1)[0]
    elem_SE.loc[SElayer_RV_idx, 5] = int(_SE_RV) + (
        int(_RVmyo) - elem_SE.loc[SElayer_RV_idx, 5]
    )
    elem_SE.loc[
        np.where(elem_SE[5] == int(_SE_RV) - (int(_LVmyo) - int(_RVmyo)))[0], 5
    ] = int(_SE_septRV)

    # overwrite new biv mesh
    IO.write_elem(SElayermeshpath_i, elem_SE)

    # insert biv mesh with SElayer into torso mesh
    if not meshtool.insert_submesh(
        meshtoolloc, meshfolder + "/" + meshname, SElayermeshpath_i, SElayermeshpath
    ):
        raise Exception("meshtool error. Abort!")

    # extract SElayer from biv
    SElayermeshpath_i_layer = mesh_iloc + "/" + SElayermeshname + "_i_layer"
    if not meshtool.extract_mesh(
        meshtoolloc,
        SElayermeshpath_i,
        _SE_septRV + "," + _SE_RV + "," + _SE_LV,
        SElayermeshpath_i_layer,
    ):
        raise Exception("meshtool error. Abort!")

    # extract left layer from SElayer
    SElayermeshpath_i_layer_LV = mesh_iloc + "/" + SElayermeshname + "_i_layer_LV"
    if not meshtool.extract_mesh(
        meshtoolloc, SElayermeshpath_i_layer, _SE_LV, SElayermeshpath_i_layer_LV
    ):
        raise Exception("meshtool error. Abort!")
    if not meshtool.extract_surface(
        meshtoolloc, SElayermeshpath_i_layer, SElayermeshpath_i_layer_LV, op=_SE_LV
    ):
        raise Exception("meshtool error. Abort!")

    # extract right layer from SElayer
    SElayermeshpath_i_layer_RV = mesh_iloc + "/" + SElayermeshname + "_i_layer_RV"
    if not meshtool.extract_mesh(
        meshtoolloc,
        SElayermeshpath_i_layer,
        _SE_septRV + "," + _SE_RV,
        SElayermeshpath_i_layer_RV,
    ):
        raise Exception("meshtool error. Abort!")
    if not meshtool.extract_surface(
        meshtoolloc,
        SElayermeshpath_i_layer,
        SElayermeshpath_i_layer_RV,
        op=_SE_septRV + "," + _SE_RV,
    ):
        raise Exception("meshtool error. Abort!")

    # extract LV SElayer base from torso mesh and map it on biv and on LVlayer
    basesurf_LVlayer = SElayermeshname + "_i.SElayer_base_left"
    if not meshtool.extract_surface(
        meshtoolloc,
        SElayermeshpath,
        surfloc + "/" + basesurf_LVlayer,
        op=_SE_LV + ":" + _upper_LH,
    ):
        raise Exception("meshtool error. Abort!")
    if not meshtool.map_file(
        meshtoolloc,
        SElayermeshpath_i,
        surfloc + "/" + basesurf_LVlayer + ".surf",
        surf_iloc,
    ):
        raise Exception("meshtool error. Abort!")
    if not meshtool.map_file(
        meshtoolloc,
        SElayermeshpath_i_layer,
        surf_iloc + "/" + basesurf_LVlayer + ".surf",
        surf_iloc_layer,
    ):
        raise Exception("meshtool error. Abort!")
    if not meshtool.map_file(
        meshtoolloc,
        SElayermeshpath_i_layer_LV,
        surf_iloc_layer + "/" + basesurf_LVlayer + ".surf",
        surf_left_loc,
    ):
        raise Exception("meshtool error. Abort!")

    # extract RV SElayer base from torso mesh and map it on biv and on RVlayer
    basesurf_RVlayer = SElayermeshname + "_i.SElayer_base_right"
    if not meshtool.extract_surface(
        meshtoolloc,
        SElayermeshpath,
        surfloc + "/" + basesurf_RVlayer,
        op=_SE_RV + ":" + _upper_RH,
    ):
        raise Exception("meshtool error. Abort!")
    if not meshtool.map_file(
        meshtoolloc,
        SElayermeshpath_i,
        surfloc + "/" + basesurf_RVlayer + ".surf",
        surf_iloc,
    ):
        raise Exception("meshtool error. Abort!")
    if not meshtool.map_file(
        meshtoolloc,
        SElayermeshpath_i_layer,
        surf_iloc + "/" + basesurf_RVlayer + ".surf",
        surf_iloc_layer,
    ):
        raise Exception("meshtool error. Abort!")
    if not meshtool.map_file(
        meshtoolloc,
        SElayermeshpath_i_layer_RV,
        surf_iloc_layer + "/" + basesurf_RVlayer + ".surf",
        surf_right_loc,
    ):
        raise Exception("meshtool error. Abort!")

    # extract left heart (LH)
    meshpath_LHlayer = apexloc + "/SElayer_LH"
    if not meshtool.extract_mesh(
        meshtoolloc, SElayermeshpath, _SE_LV + "," + _upper_LH, meshpath_LHlayer
    ):
        raise Exception("meshtool error. Abort!")

    # extract upper LH and compute distance from left SElayer (in upperLH)
    meshpath_upperLHlayer = apexloc + "/SElayer_upperLH"
    if not meshtool.extract_surface(
        meshtoolloc, meshpath_LHlayer, meshpath_upperLHlayer, op=_upper_LH
    ):
        raise Exception("meshtool error. Abort!")
    if not meshtool.generate_distancefield(
        meshtoolloc,
        meshpath_LHlayer,
        meshpath_upperLHlayer,
        apexloc + "/dist_LHlayerfromLVlayer.dat",
    ):
        raise Exception("meshtool error. Abort!")

    # find apex as furthest point in LV from everything else in LH
    dist_from_LVlayer = np.genfromtxt(apexloc + "/dist_LHlayerfromLVlayer.dat")
    idx_apex_LHlayer = np.where(dist_from_LVlayer == np.amax(dist_from_LVlayer))[0]
    pts_LVlayer = np.array(IO.read_pts(SElayermeshpath_i_layer_LV))
    pts_LHlayer = np.array(IO.read_pts(meshpath_LHlayer))
    idx_apex_LVlayer = np.where(
        np.all(pts_LVlayer == pts_LHlayer[idx_apex_LHlayer], axis=1)
    )[0]

    apexsurf_LVlayer = SElayermeshname + "_i.SElayer_lvapex"
    IO.write_vtx(surf_left_loc + "/" + apexsurf_LVlayer, pd.DataFrame(idx_apex_LVlayer))

    # find surface around apex point in LVlayer
    if not meshtool.extract_surface(
        meshtoolloc, SElayermeshpath_i_layer_LV, apexloc + "/LVlayermyo", op=_SE_LV
    ):
        raise Exception("meshtool error. Abort!")

    surf_LVlayer = IO.read_surf(apexloc + "/LVlayermyo")
    surf_apex_LVlayer = surf_LVlayer.loc[
        (surf_LVlayer[1] == int(idx_apex_LVlayer))
        | (surf_LVlayer[2] == int(idx_apex_LVlayer))
        | (surf_LVlayer[3] == int(idx_apex_LVlayer)),
        :,
    ]
    surf_apex_LVlayer_vtx = np.unique(surf_apex_LVlayer.loc[:, 1:3])
    IO.write_vtx(
        surf_left_loc + "/" + apexsurf_LVlayer + ".surf",
        pd.DataFrame(surf_apex_LVlayer_vtx),
    )
    IO.write_surf(
        surf_left_loc + "/" + apexsurf_LVlayer, pd.DataFrame(surf_apex_LVlayer)
    )

    # compute distance to apex in LV SElayer, and then map it to biv SElayer
    data_LV = SIMID + "/SElayer_Z_LV"
    if not meshtool.generate_distancefield(
        meshtoolloc,
        SElayermeshpath_i_layer_LV,
        surf_left_loc + "/" + apexsurf_LVlayer,
        data_LV + ".dat",
        "-esurf=" + surf_left_loc + "/" + basesurf_LVlayer,
    ):
        raise Exception("meshtool error. Abort!")
    if not meshtool.insert_data(
        meshtoolloc,
        SElayermeshpath_i_layer_LV,
        data_LV + ".dat",
        SElayermeshpath_i_layer,
        data_LV + "_biv.dat",
        "0",
    ):
        raise Exception("meshtool error. Abort!")

    # find apex as nearest point in RV to LV-apex
    pts_RVlayer = np.array(IO.read_pts(SElayermeshpath_i_layer_RV))
    dist_RVlayerapex = np.sqrt(
        (pts_RVlayer[:, 0] - pts_LVlayer[idx_apex_LVlayer, 0]) ** 2
        + (pts_RVlayer[:, 1] - pts_LVlayer[idx_apex_LVlayer, 1]) ** 2
        + (pts_RVlayer[:, 2] - pts_LVlayer[idx_apex_LVlayer, 2]) ** 2
    )
    idx_apex_RVlayer = np.argmin(dist_RVlayerapex)

    apexsurf_RVlayer = SElayermeshname + "_i.SElayer_rvapex"
    IO.write_vtx(
        surf_right_loc + "/" + apexsurf_RVlayer, pd.DataFrame([idx_apex_RVlayer])
    )

    # find surface around apex point in RVlayer
    if not meshtool.extract_surface(
        meshtoolloc,
        SElayermeshpath_i_layer_RV,
        apexloc + "/RVlayermyo",
        op=_SE_septRV + "," + _SE_RV,
    ):
        raise Exception("meshtool error. Abort!")

    surf_RVlayer = IO.read_surf(apexloc + "/RVlayermyo")
    surf_apex_RVlayer = surf_RVlayer.loc[
        (surf_RVlayer[1] == int(idx_apex_RVlayer))
        | (surf_RVlayer[2] == int(idx_apex_RVlayer))
        | (surf_RVlayer[3] == int(idx_apex_RVlayer)),
        :,
    ]
    surf_apex_RVlayer_vtx = np.unique(surf_apex_RVlayer.loc[:, 1:3])
    IO.write_vtx(
        surf_right_loc + "/" + apexsurf_RVlayer + ".surf",
        pd.DataFrame(surf_apex_RVlayer_vtx),
    )
    IO.write_surf(
        surf_right_loc + "/" + apexsurf_RVlayer, pd.DataFrame(surf_apex_RVlayer)
    )

    # compute distance to apex in RV SElayer, and then map it to biv SElayer
    data_RV = SIMID + "/SElayer_Z_RV"
    if not meshtool.generate_distancefield(
        meshtoolloc,
        SElayermeshpath_i_layer_RV,
        surf_right_loc + "/" + apexsurf_RVlayer,
        data_RV + ".dat",
        "-esurf=" + surf_right_loc + "/" + basesurf_RVlayer,
    ):
        raise Exception("meshtool error. Abort!")

    if not meshtool.insert_data(
        meshtoolloc,
        SElayermeshpath_i_layer_RV,
        data_RV + ".dat",
        SElayermeshpath_i_layer,
        data_RV + "_biv.dat",
        "0",
    ):
        raise Exception("meshtool error. Abort!")

    # merge Z data on biv
    uvc_LV = np.genfromtxt(data_LV + "_biv.dat")
    uvc_RV = np.genfromtxt(data_RV + "_biv.dat")
    Z_SElayer = uvc_LV + uvc_RV

    # save results
    IO.write_data(uvcfolder + "/Z_subendocardium_SE", Z_SElayer)

    # insert
    if not meshtool.insert_data(
        meshtoolloc,
        SElayermeshpath_i_layer,
        uvcfolder + "/Z_subendocardium_SE.dat",
        SElayermeshpath_i,
        uvcfolder + "/Z_subendocardium.dat",
        "0",
    ):
        raise Exception("meshtool error. Abort!")

    # ==== SElayer REDUCED EXTENSION ==== #
    # extract biventricle from SElayer torso mesh
    finalmesh_i = meshloc + "/" + finalmeshname + "_i"
    if not meshtool.extract_mesh(
        meshtoolloc,
        SElayermeshpath,
        _RVmyo + "," + _LVmyo + "," + _SE_septRV + "," + _SE_RV + "," + _SE_LV,
        finalmesh_i,
    ):
        raise Exception("meshtool error. Abort!")

    # extract SElayer from biventricle mesh
    finalmesh_i_SE = meshloc + "/" + finalmeshname + "_i_SE"
    if not meshtool.extract_mesh(
        meshtoolloc,
        finalmesh_i,
        _SE_septRV + "," + _SE_RV + "," + _SE_LV,
        finalmesh_i_SE,
    ):
        raise Exception("meshtool error. Abort!")
    elem_SE = IO.read_elem(finalmesh_i_SE)

    # find RV and LV indices in biventricle
    idxRV = np.where(elem_SE[5] == int(_SE_RV))[0]
    idxLV = np.where(elem_SE[5] == int(_SE_LV))[0]
    idxseptRV = np.where(elem_SE[5] == int(_SE_septRV))[0]

    # set tags of layers such that Z > Zmax_SE back to the corresponding myocardial tags
    subindexall = np.where(100 * Z_SElayer > round(100 * Zmax_SE))[0]
    subSElayer = elem_SE[[1, 2, 3, 4]].map(lambda x: x in subindexall)
    subSE_RVidx = np.intersect1d(np.where(subSElayer.any(axis=1).map(int) == 1), idxRV)
    subSE_LVidx = np.intersect1d(np.where(subSElayer.any(axis=1).map(int) == 1), idxLV)
    subSE_septRVidx = np.intersect1d(
        np.where(subSElayer.any(axis=1).map(int) == 1), idxseptRV
    )
    elem_SE.loc[subSE_RVidx, 5] = _RVmyo
    elem_SE.loc[subSE_LVidx, 5] = _LVmyo
    elem_SE.loc[subSE_septRVidx, 5] = _LVmyo

    # overwrite finalmesh
    IO.write_elem(finalmesh_i_SE, elem_SE)

    # insert submeshes
    if not meshtool.insert_submesh(
        meshtoolloc, finalmesh_i, finalmesh_i_SE, finalmesh_i
    ):
        raise Exception("meshtool error. Abort!")
    if not meshtool.insert_submesh(
        meshtoolloc, SElayermeshpath, finalmesh_i, finalmeshfolder + "/" + finalmeshname
    ):
        raise Exception("meshtool error. Abort!")

    # delete unnecessary files
    if cleandir:
        shutil.rmtree(simresult)

    return


def compute_FEC_layers(
    CARP_INSTALL_DIR,
    meshfolder,
    meshname,
    finalmeshfolder,
    finalmeshname,
    uvcfolder,
    _RVmyo="34",
    _LVmyo="35",
    _SE_RV="44",
    _SE_LV="45",
    _SE_septRV="46",
    factor_tag_FEC=100,
    Z_FEC_min=0.50,
    Z_FEC_max=0.75,
    deltaZ_FEC=0.02,
    cleandir=False,
):
    """Compute multiple fast endocardial layers for biventricle. NOTE: requires running extract_SElayer first!

    Parameters
    ----------
    CARP_INSTALL_DIR: str
        Path to CARPentry installation directory.
    meshfolder: str
        Path to the computational mesh.
    meshname: str
        Name of the torso computational mesh with subendocardial layers (from apex to base).
    finalmeshfolder: str
        Path to the directory where the final computational mesh will be saved.
    finalmeshname: str
        Name of the torso computational mesh with new tags for the subendocardial layer.
    uvcfolder : str
        Path to folder with universal ventricular coordinates.
    Z: str
        Path to .dat file with apico-basal UVC on subendocardial layer only.
    _RVmyo: str
        Pre-existing tag for RV myocardium. Default: "34".
    _LVmyo: str
        Pre-existing tag for LV myocardium. Default: "35".
    _SE_RV: str
        Pre-existing tag for RV subendocardial layer. Default: "44".
    _SE_LV: str
        Pre-existing tag for LV subendocardial layer. Default: "45".
    _SE_septRV: str
        Pre-existing tag for septal subendocardial layer on the right side. Default: "46".
    factor_tag_FEC: int
        Multiplicative factor for fast endocardial conductive layers with respect corresponding subendocardial layer tag. Default: 100.
    Z_FEC_min: float
        Lower bound for fast endocardial conductive layer extension. Must be in [0,1]. Default: 0.50.
    Z_FEC_max: float
        Upper bound for fast endocardial conductive layer extension. Must be in [0,1]. Default: 0.75.
    deltaZ_FEC: float
        Apico-basal extension of each fast endocardial conducting layer. Default: 0.02.
    cleandir: bool
        If True, remove simulation output directory. Default: False.
    """

    # find executables
    meshtoolloc = CARP_INSTALL_DIR + "/meshtool/meshtool"

    # check value for Z_FEC_min and Z_FEC_max
    if Z_FEC_min < 0 or Z_FEC_min > 1:
        raise Exception(
            "Minimal extension for fast endocardial conductive layer must be in [0,1]"
        )
    if Z_FEC_max < 0 or Z_FEC_max > 1:
        raise Exception(
            "Maximal extension for fast endocardial conductive layer must be in [0,1]"
        )
    if Z_FEC_min > Z_FEC_max:
        raise Exception(
            "Maximal extension for fast endocardial conductive layer must be bigger than minimal extension"
        )

    # IO directories
    simresult = finalmeshfolder + "/FEC"
    if not os.path.exists(simresult):
        os.mkdir(simresult)

    # volumetric meshes
    meshloc = simresult + "/volumes"
    if not os.path.exists(meshloc):
        os.mkdir(meshloc)

    # volumetric meshes extracted from biv
    mesh_iloc = meshloc + "/biv"
    if not os.path.exists(mesh_iloc):
        os.mkdir(mesh_iloc)

    # extract biventricle (including SElayer) from torso mesh
    finalmesh_i = meshloc + "/" + finalmeshname + "_i"
    if not meshtool.extract_mesh(
        meshtoolloc,
        meshfolder + "/" + meshname,
        _RVmyo + "," + _LVmyo + "," + _SE_septRV + "," + _SE_RV + "," + _SE_LV,
        finalmesh_i,
    ):
        raise Exception("meshtool error. Abort!")

    # extract SElayer from biventicular mesh
    finalmesh_i_SE = mesh_iloc + "/" + finalmeshname + "_i_SE"
    if not meshtool.extract_mesh(
        meshtoolloc,
        finalmesh_i,
        _SE_septRV + "," + _SE_RV + "," + _SE_LV,
        finalmesh_i_SE,
    ):
        raise Exception("meshtool error. Abort!")
    elem_SE = IO.read_elem(finalmesh_i_SE)

    # find RV and LV indices in biventricle
    idxRV = np.where(elem_SE[5] == int(_SE_RV))[0]
    idxLV = np.where(elem_SE[5] == int(_SE_LV))[0]
    idxseptRV = np.where(elem_SE[5] == int(_SE_septRV))[0]

    # load apico-basal coordinate on SElayer
    Z_SElayer = np.genfromtxt(uvcfolder + "/Z_subendocardium_SE.dat")

    # set tags of layers such that Z > Zmax_SE back to the corresponding myocardial tags
    subindexall = np.where(100 * Z_SElayer > round(100 * Z_FEC_max))[0]
    subSElayer = elem_SE[[1, 2, 3, 4]].map(lambda x: x in subindexall)
    subSE_RVidx = np.intersect1d(np.where(subSElayer.any(axis=1).map(int) == 1), idxRV)
    subSE_LVidx = np.intersect1d(np.where(subSElayer.any(axis=1).map(int) == 1), idxLV)
    subSE_septRVidx = np.intersect1d(
        np.where(subSElayer.any(axis=1).map(int) == 1), idxseptRV
    )
    elem_SE.loc[subSE_RVidx, 5] = _RVmyo
    elem_SE.loc[subSE_LVidx, 5] = _LVmyo
    elem_SE.loc[subSE_septRVidx, 5] = _LVmyo

    # define multiple SE layers
    print("Tagging multiple FEC layers...")
    for layer in fec_layers(Z_FEC_min, Z_FEC_max, deltaZ_FEC):
        print("... {}% apico-basal extension".format(layer))
        subindexall = np.intersect1d(
            np.where(100 * Z_SElayer >= layer),
            np.where(100 * Z_SElayer < layer + 100 * deltaZ_FEC),
        )
        subSElayer = elem_SE[[1, 2, 3, 4]].map(lambda x: x in subindexall)
        subSE_RVidx = np.intersect1d(
            np.where(subSElayer.any(axis=1).map(int) == 1), idxRV
        )
        subSE_LVidx = np.intersect1d(
            np.where(subSElayer.any(axis=1).map(int) == 1), idxLV
        )
        subSE_septRVidx = np.intersect1d(
            np.where(subSElayer.any(axis=1).map(int) == 1), idxseptRV
        )
        elem_SE.loc[subSE_RVidx, 5] = str(factor_tag_FEC * int(_SE_RV) + layer)
        elem_SE.loc[subSE_LVidx, 5] = str(factor_tag_FEC * int(_SE_LV) + layer)
        elem_SE.loc[subSE_septRVidx, 5] = str(factor_tag_FEC * int(_SE_septRV) + layer)
    IO.write_elem(finalmesh_i_SE, elem_SE)

    # insert submeshes
    if not meshtool.insert_submesh(
        meshtoolloc, finalmesh_i, finalmesh_i_SE, finalmesh_i
    ):
        raise Exception("meshtool error. Abort!")
    if not meshtool.insert_submesh(
        meshtoolloc,
        meshfolder + "/" + meshname,
        finalmesh_i,
        finalmeshfolder + "/" + finalmeshname,
    ):
        raise Exception("meshtool error. Abort!")

    # delete unnecessary files
    if cleandir:
        shutil.rmtree(simresult)

    return


def segment_submesh(
    CARP_INSTALL_DIR,
    meshfolder,
    meshname,
    finalmeshfolder,
    finalmeshname,
    uvcfolder,
    _submsh_myo,
    _submsh_SE,
    _submsh_septSE,
    factor_tag_myo,
    factor_tag_SE,
    Nz,
    Nrho,
):
    """Segment biventricular submesh (e.g., RV, LV) into different regions along apico-basal and tranmural direction.

    Parameters
    ----------
    CARP_INSTALL_DIR: str
        Path to CARPentry installation directory.
    meshfolder: str
        Path to the computational mesh.
    meshname: str
        Name of the biventricular computational mesh, eventually with fast endocardial conductive layer(s).
    finalmeshfolder: str
        Path to the directory where the final computational mesh will be saved.
    finalmeshname: str
        Name of the biventricular submesh with new tags for each region.
    uvcfolder : str
        Path to folder with universal ventricular coordinates.
    _submsh_myo : str
        Pre-exising tag(s) for submesh myocardium.
    _submsh_SE : str
        Pre-exising tag(s) for submesh fast conducting endocardial layer(s). Use meshtool convention for multiple tags.
    _submsh_septSE : str
        Pre-exising tag(s) for submesh fast conducting endocardial layer(s) located on the septum. Use meshtool convention for multiple tags. If any, set "".
    factor_tag_myo : int
        Multiplicative factor for segmented regions with respect corresponding myocardium tag.
    factor_tag_septSE : int
        Multiplicative factor for segmented regions with respect corresponding fast endocardial conductive layer tag.
    Nz : int
        Number of regions along the apico-basal direction.
    Nrho : int
        Number of regions along the transmural direction.
    """

    # find executables
    meshtoolloc = CARP_INSTALL_DIR + "/meshtool/meshtool"

    # set submsh tag as the union of myocardium and subendocardial layer tags
    _submsh = _submsh_myo + "," + _submsh_SE + "," + _submsh_septSE

    # extract submesh (i.e., RV or LV) myocardium from biv
    if not meshtool.extract_mesh(
        meshtoolloc,
        meshfolder + "/" + meshname,
        _submsh,
        finalmeshfolder + "/" + finalmeshname,
    ):
        raise Exception("meshtool error. Abort!")
    elem_submsh = IO.read_elem(finalmeshfolder + "/" + finalmeshname)

    # extract Z and RHO of current submsh
    if not meshtool.extract_data(
        meshtoolloc,
        finalmeshfolder + "/" + finalmeshname,
        uvcfolder + "/Z.dat",
        finalmeshfolder + "/" + finalmeshname + "_Z.dat",
        "0",
    ):
        raise Exception("meshtool error. Abort!")
    if not meshtool.extract_data(
        meshtoolloc,
        finalmeshfolder + "/" + finalmeshname,
        uvcfolder + "/RHO.dat",
        finalmeshfolder + "/" + finalmeshname + "_RHO.dat",
        "0",
    ):
        raise Exception("meshtool error. Abort!")

    # load UVCs
    Zsubmsh = np.genfromtxt(finalmeshfolder + "/" + finalmeshname + "_Z.dat")
    RHOsubmsh = np.genfromtxt(finalmeshfolder + "/" + finalmeshname + "_RHO.dat")

    # get subintervals in each direction
    Zlevel = np.linspace(min(Zsubmsh), max(Zsubmsh), Nz + 1)
    RHOlevel = np.linspace(min(RHOsubmsh), max(RHOsubmsh), Nrho + 1)

    # tag and idx submsh
    tagsubmsh = elem_submsh[5]
    idxsubmsh = np.array(elem_submsh[[1, 2, 3, 4]])

    # new tags
    tags_new = copy.deepcopy(tagsubmsh)
    for el in range(len(elem_submsh)):
        tt = idxsubmsh[el, :]
        z = 0
        while z < Nz and np.min(Zsubmsh[tt]) >= Zlevel[z]:
            z = z + 1
        rho = 0
        while rho < Nrho and np.min(RHOsubmsh[tt]) >= RHOlevel[rho]:
            rho = rho + 1
        if str(tagsubmsh[el]) in _submsh_myo.split(","):
            tags_new[el] = factor_tag_myo * tagsubmsh[el] + z + (rho - 1) * Nz
        elif str(tagsubmsh[el]) in _submsh_SE.split(","):
            tags_new[el] = factor_tag_SE * tagsubmsh[el] + z
        elif str(tagsubmsh[el]) in _submsh_septSE.split(","):
            tags_new[el] = factor_tag_SE * tagsubmsh[el] + z
    elem_submsh[5] = tags_new

    # Save submesh with different segments for gregion
    IO.write_elem(finalmeshfolder + "/" + finalmeshname, elem_submsh)

    # insert submesh back into biv
    if not meshtool.insert_submesh(
        meshtoolloc,
        meshfolder + "/" + meshname,
        finalmeshfolder + "/" + finalmeshname,
        meshfolder + "/" + meshname,
    ):
        raise Exception("meshtool error. Abort!")

    return


def segment_biventricle(
    CARP_INSTALL_DIR,
    meshfolder,
    meshname,
    finalmeshfolder,
    finalmeshname,
    uvcfolder,
    _RVmyo="34",
    _LVmyo="35",
    _SE_RV="44",
    _SE_LV="45",
    _SE_septRV="46",
    factor_tag_myo=100,
    factor_tag_SE=100,
    Nz=10,
    Nrho_LV=5,
    Nrho_RV=3,
    cleandir=False,
):
    """Segment biventricular myocardium (including FEC layer(s), if any) into different regions along apico-basal and tranmural direction.

    Parameters
    ----------
    CARP_INSTALL_DIR: str
        Path to CARPentry installation directory.
    meshfolder: str
        Path to the computational mesh.
    meshname: str
        Name of the biventricular computational mesh, eventually with fast endocardial conductive layer(s).
    finalmeshfolder: str
        Path to the directory where the final computational mesh will be saved.
    finalmeshname: str
        Name of the torso computational mesh with new tags for each region.
    uvcfolder : str
        Path to folder with universal ventricular coordinates.
    _RVmyo: str
        Pre-existing tag for RV myocardium. Default: "34".
    _LVmyo: str
        Pre-existing tag for LV myocardium. Default: "35".
    _SE_RV: str
        Pre-existing tag for RV subendocardial layer. Default: "44".
    _SE_LV: str
        Pre-existing tag for LV subendocardial layer. Default: "45".
    _SE_septRV: str
        Pre-existing tag for septal subendocardial layer on the right side. Default: "46".
    factor_tag_myo : int
        Multiplicative factor for segmented regions with respect corresponding myocardium tag. Default: 100.
    factor_tag_SE : int
        Multiplicative factor for segmented regions with respect corresponding subendocardial layer(s) tag. Default: 100.
    Nz : int
        Number of regions along the apico-basal direction. Maximum value: 10.
    Nrho_LV : int
        Number of regions along the transmural direction in the LV. Maximum value: 5.
    Nrho_RV : int
        Number of regions along the transmural direction in the RV. Maximum value: 5.
    cleandir: bool
        If True, remove simulation output directory. Default: False.
    """

    # find executables
    meshtoolloc = CARP_INSTALL_DIR + "/meshtool/meshtool"

    # checks
    if Nz > 10:
        raise Exception(
            "Number of regions along apico-basal direction must be at most 10. If a higher number is desired, suppress this exception. NOTE: CARP does not allow more than 256 different regions!"
        )
    if Nrho_LV > 5 or Nrho_RV > 5:
        raise Exception(
            "Number of regions along transmural direction must be at most 5. If a higher number is desired, suppress this exception. NOTE: CARP does not allow more than 256 different regions!"
        )

    # IO directories
    simresult = finalmeshfolder + "/SEGMENTATION"
    if not os.path.exists(simresult):
        os.mkdir(simresult)

    # volumetric meshes
    meshloc = simresult + "/volumes"
    if not os.path.exists(meshloc):
        os.mkdir(meshloc)

    # volumetric meshes extracted from biv
    mesh_iloc = meshloc + "/biv"
    if not os.path.exists(mesh_iloc):
        os.mkdir(mesh_iloc)

    # extract biv (myo+SElayer(s)) from torso
    finalmesh_i = finalmeshname + "_i"
    if not meshtool.extract_mesh(
        meshtoolloc,
        meshfolder + "/" + meshname,
        _RVmyo + "," + _LVmyo + "," + _SE_RV + "," + _SE_LV + "," + _SE_septRV,
        meshloc + "/" + finalmesh_i,
    ):
        raise Exception("meshtool error. Abort!")

    # segment each ventricle independently
    segment_submesh(
        CARP_INSTALL_DIR,
        meshloc,
        finalmesh_i,
        mesh_iloc,
        "RVmyo",
        uvcfolder,
        _RVmyo,
        _SE_RV,
        _SE_septRV,
        factor_tag_myo,
        factor_tag_SE,
        Nz,
        Nrho_RV,
    )

    segment_submesh(
        CARP_INSTALL_DIR,
        meshloc,
        finalmesh_i,
        mesh_iloc,
        "LVmyo",
        uvcfolder,
        _LVmyo,
        _SE_LV,
        "",
        factor_tag_myo,
        factor_tag_SE,
        Nz,
        Nrho_LV,
    )

    # insert biv into torso
    if not meshtool.insert_submesh(
        meshtoolloc,
        meshfolder + "/" + meshname,
        meshloc + "/" + finalmesh_i,
        finalmeshfolder + "/" + finalmeshname,
    ):
        raise Exception("meshtool error. Abort!")

    # delete unnecessary files
    if cleandir:
        shutil.rmtree(simresult)

    return
