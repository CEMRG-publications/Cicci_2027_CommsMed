#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file is part of 'Cicci_2027_CommsMed'.

__author__ = Ludovica Cicci
__email__ = l.cicci@imperial.ac.uk
"""

# extra imports
import os
import json
import pandas as pd
import numpy as np
from natsort import natsorted


# --- TXT FILES --- #
def read_text(filename):
    """Read .txt file.

    Parameters
    ----------
    filename : str
        Path to txt file (without extension).
    """

    data = []
    with open(filename + ".txt", "r") as file:
        for line in file:
            data.append(line.replace("\n", ""))
        print("Successfully read {}".format(filename + ".txt"))

    return data


def read_labels(filename):
    """Read .txt file with the labels of input parameters, in alphabetical order.

    Labels are sorted by character codes (capital letters first), i.e., in the same order as the columns of the sample matrices.
    Use read_text for label files that are paired with data (e.g., xlabels of a data folder, ThetaArray_labels) to keep their order.

    Parameters
    ----------
    filename : str
        Path to txt file (without extension).
    """

    labels = read_text(filename)

    if "" in labels or any(label != label.strip() for label in labels):
        raise Exception(
            "{}.txt contains empty lines or labels with spaces. Remove them.".format(
                filename
            )
        )

    duplicates = sorted({label for label in labels if labels.count(label) > 1})
    if duplicates:
        raise Exception(
            "{}.txt contains duplicate labels: {}".format(filename, duplicates)
        )

    return sorted(labels)


def read_float_text(filename):
    """Read .txt file.

    Parameters
    ----------
    filename : str
        Path to txt file (without extension).
    """

    data = np.genfromtxt(filename + ".txt")
    print("Successfully read {}".format(filename + ".txt"))

    return data


def write_string_text(filename, data):
    """Write .txt file.

    Parameters
    ----------
    filename : str
        Path to txt file (without extension).
    data: str
        File content.
    """
    with open(filename + ".txt", "w") as file:
        for i in data:
            file.write("%s" % i)
        file.close()
        print("Data written to file {}.".format(filename + ".txt"))

    return


def write_float_text(filename, data):
    """Write .txt file.

    Parameters
    ----------
    filename : str
        Path to txt file (without extension).
    data: numpy array
        File content.
    """
    with open(filename + ".txt", "w") as file:
        for i in data:
            file.write("%f\n" % i)
        file.close()
        print("Data written to file {}.".format(filename + ".txt"))

    return


def write_array_text(filename, data, fmt="%.6e"):
    """Write .txt file.

    Parameters
    ----------
    filename : str
        Path to txt file (without extension).
    data: numpy array
        File content.
    """
    np.savetxt(filename + ".txt", data, fmt=fmt)
    print("Data written to file {}.".format(filename + ".txt"))

    return




# --- DAT FILES --- #
def read_data(filename):
    """Read .dat file.

    Parameters
    ----------
    filename:
        Path to dat file (without extension).
    """

    data = np.genfromtxt(filename + ".dat")
    print("Successfully read {}".format(filename + ".dat"))

    return data


def write_2D_data(filename, data):
    """Write .dat file.

    Parameters
    ----------
    filename:
        Path to dat file (without extension).
    data: numpy array
        File content.
    """

    data = np.savetxt(filename + ".dat", data)
    print("Data written to file {}".format(filename + ".dat"))

    return data


def write_data(filename, data):
    """Write .dat file.

    Parameters
    ----------
    filename:
        Path to dat file (without extension).
    data: numpy array
        File content.
    """

    with open(filename + ".dat", "w") as file:
        for i in data:
            file.write("%f\n" % i)
        file.close()
        print("Data written to file {}.".format(filename + ".dat"))

    return


def read_ecg_dat(filename):
    """Read 12-lead ECG .dat file.

    Parameters
    ----------
    filename : str
        Path to dat file (with extension). Columns: time, LI, LII, LIII, aVR, aVL, aVF, V1, V2, V3, V4, V5, V6.

    Returns
    -------
    ecg : pandas DataFrame
        One column per lead, indexed by time.
    """

    ecgdata = np.loadtxt(filename, dtype=float)
    leads = ["LI", "LII", "LIII", "aVR", "aVL", "aVF"] + [
        "V" + str(i) for i in range(1, 7)
    ]
    ecg = pd.DataFrame(ecgdata[:, 1:13], columns=leads, index=ecgdata[:, 0])
    ecg.index.name = "t"

    return ecg


# --- JSON FILES ---#
def read_json(filename):
    """Read .json file.

    Parameters
    ----------
    filename : str
        Path to json file (without extension).
    """

    # Retrieve JSON data from the file
    with open(filename + ".json", "r") as file:
        data = json.load(file)
        print("Successfully read {}".format(filename + ".json"))

    return data


def write_json(filename, data):
    """Write .json file.

    Parameters
    ----------
    filename:
        Path to json file (without extension).
    data: dict
        File content.
    """

    with open(filename + ".json", "w") as file:
        json.dump(data, file, indent="\t")
        file.close()
        print("Data written to file {}.".format(filename + ".json"))

    return


def append_json(filename, data):
    """Append data to a pre-existing .json file.

    Parameters
    ----------
    filename:
        Path to json file (without extension).
    data: dict
        File content.
    """

    with open(filename + ".json", "r") as file:
        predata = json.load(file)
    predata.update(data)
    write_json(filename, predata)

    return


# --- PTS FILES ---#
def read_pts(filename):
    """Read mesh vertices (.pts).

    Parameters
    ----------
    filename : str
        Path to pts file (without extension).
    """

    # read mesh files
    data = pd.read_csv(filename + ".pts", sep=" ", skiprows=1, header=None)
    print("Successfully read {}".format(filename + ".pts"))
    print("Mesh has {} nodes".format(len(data)))

    return data


# This function writes a pts file with a pre-defined header (from meshIO.write_pts)


# --- SURF FILES --- #
def read_surf(filename):
    """Read mesh surfaces (.surf).

    Parameters
    ----------
    filename : str
        Path to surf file (without extension).
    """

    surf = pd.read_csv(filename + ".surf", sep=" ", skiprows=1, header=None)
    print("Successfully read {}".format(filename + ".surf"))
    print("Mesh has {} surf".format(len(surf)))

    return surf


def read_surf_to_nodeList(filename):
    """Read surface and output a list of nodes (.surf).

    Parameters
    ----------
    filename : str
        Path to surf file (without extension).
    """

    surf = read_surf(filename)

    surf_nodes = []
    surf_nodes = np.append(surf[1], surf_nodes)
    surf_nodes = np.append(surf[2], surf_nodes)
    surf_nodes = np.append(surf[3], surf_nodes)
    surf_nodes = np.unique(surf_nodes)

    surf_nodes = surf_nodes.astype(int)

    return surf_nodes


def write_surf(filename, surf, shapes="Tr"):
    """Write surface mesh vertices (.surf).

    Parameters
    ----------
    filename : str
        Path to surf file (without extension).
    surf : DataFrame
        Mesh surface.
    shapes : str
        Shape of surface elements. Default: Tr (triangles).
    """

    with open(filename + ".surf", "w") as file:
        file.write("{}\n".format(len(surf)))
    if surf.iloc[0, 0] != "Tr":
        surf.insert(loc=0, value=shapes, column=0)
    surf.to_csv(filename + ".surf", sep=" ", header=False, index=False, mode="a")
    print("surf data written to file {}.".format(filename + ".surf"))
    del surf[0]  # remove added column to prevent cross-talk problems later
    file.close()

    return


# --- ELEM FILES --- #
def read_elem(filename):
    """Read mesh elements (.elem).

    Parameters
    ----------
    filename : str
        Path to elem file (without extension).
    """

    elem = pd.read_csv(
        filename + ".elem", sep=" ", skiprows=1, usecols=(1, 2, 3, 4, 5), header=None
    )
    print("Successfully read {}".format(filename + ".elem"))
    print("Mesh has {} elements".format(len(elem)))

    return elem


def write_elem(filename, elem, shapes="Tt"):
    """Write mesh elements (.elem).

    Parameters
    ----------
    filename : str
        Path to elem file (without extension).
    elem : DataFrame
        Mesh volume.
    shapes : str
        Shape of surface elements. Default: Tr (triangles).
    """

    with open(filename + ".elem", "w") as file:
        file.write("{}\n".format(len(elem)))
    elem.insert(loc=0, value=shapes, column=0, allow_duplicates=True)
    # TODO: replace 'filename + ".elem"' with 'file' and check it works okay
    elem.to_csv(filename + ".elem", sep=" ", header=False, index=False, mode="a")
    print("elem data written to file {}.".format(filename + ".elem"))
    del elem[0]  # Remove added column to prevent cross-talk problems later
    file.close()

    return


# --- VTX FILES --- #


def write_vtx(filename, vtx):
    """Write mesh vertix indices (.vtx).

    Parameters
    ----------
    filename : str
        Path to vtx file (without extension).
    vtx : DataFrame
        Vertex indices.
    """
    with open(filename + ".vtx", "+w") as file:
        file.write("%i\n" % len(vtx))
        file.write("extra\n")
        for i in range(len(vtx)):
            file.write("%i\n" % vtx.loc[i])
        file.close()
        print("Data written to file {}.".format(filename + ".vtx"))

    return


# --- LON FILES --- #
def read_lon(filename):
    """Read fibres (.lon).

    Parameters
    ----------
    filename : str
        Path to lon file (without extension).
    """

    lon = pd.read_csv(filename + ".lon", sep=" ", skiprows=1, header=None)
    print("Successfully read {}".format(filename + ".lon"))

    return lon


def write_lon(filename, lon):
    """Write fibres (.lon).

    Parameters
    ----------
    filename : str
        Path to lon file (without extension).
    vtx : DataFrame
        Fibre direction.
    """

    with open(filename + ".lon", "w") as file:
        file.write("1\n")
        lon.to_csv(
            file, sep=" ", header=False, index=False, mode="a"  # filename + ".lon",
        )
        file.close()
        print("Data written to file {}.".format(filename + ".lon"))

    return


# --- CSV FILES --- #
def write_MS_csv(in_folder, X_filename, xlabels_filename, out_folder, outputfile):
    X = np.genfromtxt(in_folder + "/" + X_filename, dtype=float)
    xlabels = np.genfromtxt(in_folder + "/" + xlabels_filename, dtype=str).tolist()

    df = pd.DataFrame(X)
    if len(df.columns) == 1:
        df = df.transpose()

    df.to_csv(
        out_folder + "/" + outputfile + ".csv",
        header=xlabels,
        index=False,
        float_format="%.6e",
    )

    return


def write_MS_csv_list(
    in_folder, X_filename_list, xlabels_filename_list, out_folder, outputfile
):
    for i, X_filename in enumerate(X_filename_list):
        if i == 0:
            X = np.genfromtxt(in_folder + "/" + X_filename, dtype=float).reshape(-1, 1)
            xlabels = np.genfromtxt(
                in_folder + "/" + xlabels_filename_list[i], dtype=str
            )
        else:
            X = np.hstack(
                (
                    X,
                    np.genfromtxt(in_folder + "/" + X_filename, dtype=float).reshape(
                        -1, 1
                    ),
                )
            )
            xlabels = np.hstack(
                (
                    xlabels,
                    np.genfromtxt(
                        in_folder + "/" + xlabels_filename_list[i], dtype=str
                    ),
                )
            )

    df = pd.DataFrame(X)
    df.to_csv(
        out_folder + "/" + outputfile + ".csv",
        header=xlabels,
        index=False,
        float_format="%.6e",
    )

    return


def write_csv(in_folder, out_folder, ID=""):
    print("----------------------------------------------------")
    print("Extracting data from " + in_folder + "...")

    # write inputs csv-file
    write_MS_csv(
        in_folder + "/X", "X.txt", "xlabels.txt", out_folder, "parameters" + str(ID)
    )

    # write output(s) csv-file
    y_list = []
    ylabel_list = []
    for l in natsorted(os.listdir(in_folder + "/y")):
        if "y_" in l:
            y_list.append(l)
        if "ylabels_" in l:
            ylabel_list.append(l)
    if len(y_list) == 1:
        write_MS_csv(
            in_folder + "/y", y_list[0], ylabel_list[0], out_folder, "outputs" + str(ID)
        )
    else:
        write_MS_csv_list(
            in_folder + "/y", y_list, ylabel_list, out_folder, "outputs" + str(ID)
        )

    print("Saving to " + out_folder)

    print("Done.")

    print("----------------------------------------------------")

    return


def write_txt(in_folder, out_folder):
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)

    par = pd.read_csv(in_folder + "/parameters.csv")
    outputs = pd.read_csv(in_folder + "/outputs.csv")

    print("----------------------------------------------------")
    print("Extracting data from " + in_folder)

    x_labels = list(par.columns)
    X = par.to_numpy(copy=True)

    y_labels = list(outputs.columns)
    Y = outputs.to_numpy(copy=True)

    print("Saving to " + out_folder + "...")

    with open(out_folder + "/xlabels.txt", "w") as f:
        for xl in x_labels:
            f.write(xl + "\n")

    with open(out_folder + "/ylabels.txt", "w") as f:
        for yl in y_labels:
            f.write(yl + "\n")

    np.savetxt(out_folder + "/X.txt", X, fmt="%.6e")
    np.savetxt(out_folder + "/Y.txt", Y, fmt="%.6e")

    print("Done.")

    print("----------------------------------------------------")

    return




