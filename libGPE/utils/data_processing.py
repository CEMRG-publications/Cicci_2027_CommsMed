#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file is part of 'Cicci_2027_CommsMed'.

__author__ = Ludovica Cicci
__email__ = l.cicci@imperial.ac.uk
"""

# dependencies

# other imports
import os
import numpy as np
import pandas as pd
import re
from natsort import natsorted


# --- READING --- #
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


# --- WRITING --- #
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

    # write inputs csv-file, if any (a clinical target has no inputs)
    if os.path.exists(in_folder + "/X"):
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

    outputs = pd.read_csv(in_folder + "/outputs.csv")

    print("----------------------------------------------------")
    print("Extracting data from " + in_folder)

    y_labels = list(outputs.columns)
    Y = outputs.to_numpy(copy=True)

    print("Saving to " + out_folder + "...")

    # inputs, if any (a clinical target has no inputs)
    if os.path.exists(in_folder + "/parameters.csv"):
        par = pd.read_csv(in_folder + "/parameters.csv")
        x_labels = list(par.columns)
        X = par.to_numpy(copy=True)

        with open(out_folder + "/xlabels.txt", "w") as f:
            for xl in x_labels:
                f.write(xl + "\n")

        np.savetxt(out_folder + "/X.txt", X, fmt="%.6e")

    with open(out_folder + "/ylabels.txt", "w") as f:
        for yl in y_labels:
            f.write(yl + "\n")

    np.savetxt(out_folder + "/Y.txt", Y, fmt="%.6e")

    features_idx = np.arange(Y.shape[1])
    np.savetxt(out_folder + "/features_idx_list.txt", features_idx, fmt="%d")

    print("Done.")

    print("----------------------------------------------------")

    return


def write_summary(out_filepath, out_label, in_folder, in_filename="training_summary"):
    # initialisation
    data = []
    headers = []

    # process each directory
    for i in natsorted(os.listdir(in_folder)):
        if not ".csv" in i:
            file_path = in_folder + "/" + str(i) + "/" + in_filename + ".txt"

            metrics, values = extract_lines(file_path)
            if values:
                # append the extracted line
                data.append([int(i)] + [out_label[int(i)]] + values)

                # use the header from the first file
                if not headers:
                    headers = ["folder"] + ["output"] + metrics

    # create a DataFrame
    df = pd.DataFrame(data, columns=headers)

    # save to .csv
    df.to_csv(out_filepath + ".csv", index=False)

    return




# --- EXTRA --- #
def extract_lines(file_path, text_pattern="*** Final GPE ***"):
    with open(file_path, "r") as file:
        lines = file.readlines()

    # initialisation
    capture = False
    metrics_list = []
    values_list = []
    for i, line in enumerate(lines):
        if text_pattern in line:
            capture = True
            continue
        if capture:
            target_line = line.strip()

            # Check if it contains a string followed by a number
            if re.search(r"\w+\s+\d+", target_line):
                split_target = target_line.split(" ")
                split_target = list(filter(None, split_target))
                metrics_list.append(split_target[0])
                values_list.append(split_target[1])

    return metrics_list, values_list


def combine_csv_files(in_filefolder, in_filename, N_files):

    # read the first file
    combined_df = pd.read_csv(in_filefolder + "/" + in_filename + "_0.csv")

    # read and append data from the other files
    for i in range(1, N_files):
        df = pd.read_csv(in_filefolder + "/" + in_filename + "_" + str(i) + ".csv")
        combined_df = pd.concat([combined_df, df], ignore_index=True)

    # save the combined data to the first file
    combined_df.to_csv(in_filefolder + "/" + in_filename + ".csv", index=False)

    # delete all other files
    for i in range(N_files):
        os.remove(in_filefolder + "/" + in_filename + "_" + str(i) + ".csv")

    return
