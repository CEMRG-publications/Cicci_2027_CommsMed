#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file is part of 'Cicci_2027_CommsMed'.

__author__ = Ludovica Cicci
__email__ = l.cicci@imperial.ac.uk
"""

# dependencies
from utils import IO

# extra imports
import os
import numpy as np


def lhd(I, n_samples):
    """Latin hypercube design (sampling)

    Parameters
    ----------
    I : numpy array
        Parameter space (hypercube).
    n_samples : int
        Number of samples.

    Results
    ----------
    H : numpy array
        Sample matrix with each row corresponding to a specific sample.
    """

    n_params = I.shape[0]
    H = np.zeros((n_samples, n_params))
    for j in range(n_params):
        u = np.random.random_sample(n_samples)
        b = np.arange(n_samples)
        np.random.shuffle(b)
        d = (b + u) / n_samples
        H[:, j] = I[j, 0] + (I[j, 1] - I[j, 0]) * d

    return H




def save_fixed_parameters(
    paramfolder,
    paramfilename,
    labelfolder,
    labelfilename_list,
    outputfolder,
    N_samples=500,
):
    """Save fixed input parameters.

    Parameters
    ----------
    paramfolder: str
        Path to the electrophysiology parameter datafile.
    paramfilename: str
        Name of .json file with names and reference values of all simulation parameters (without extension).
    labelfolder: str
        Path to the input parameter datafile.
    labelfilename_list: list
        Name(s) of all .txt file with lables of parameters of interest, i.e., both scalar and spatially varying, if any (without extension).
    outputfolder: str
        Path to the directory where the sample datafile(s) will be saved.
    N_samples: int
        Number of samples. Default: 500.
    """

    # TODO: it might make more sense to change it in "in not type(labelfilename_list) == list" [same everywhere else in the library]
    if type(labelfilename_list) == str:
        labelfilename_list = [labelfilename_list]

    # set path to parameters .json file and to input labels .txt file
    paramfile = paramfolder + "/" + paramfilename
    labelfile_list = [
        labelfolder + "/" + labelfilename for labelfilename in labelfilename_list
    ]

    # check folder structure
    print("\nChecking folder structure...")
    to_check = [paramfile + ".json"] + [
        labelfile + ".txt" for labelfile in labelfile_list
    ]
    for f in to_check:
        if not os.path.exists(f):
            raise Exception("Cannot find {}. Check path to file.".format(f))
        else:
            print(f + " found.")

    # load dictionary with all simulation parameters (including those that are NOT considered as inputs)
    par_all_dct = IO.read_json(paramfile)

    # find parameters of interest as given in labelfile
    par_dct = {}
    for labelfile in labelfile_list:
        for key in IO.read_labels(labelfile):
            if key in list(par_dct.keys()):
                raise Exception(
                    "{} can't be both scalar and spatially varying. Remove it from one of the label files.".format(
                        key
                    )
                )
            par_dct[key] = par_all_dct[key]

    # TODO: check if it works correctly in all cases!!!
    # CALIBRATION: find parameters with at least one among theta, grad_rho and grad_z in labelfile
    all_spatial_inputs = []
    for key in list(par_dct.keys()):
        if "_theta" in key or "_grad_rho" in key or "_grad_z" in key:
            all_spatial_inputs.append(
                key.replace("_theta", "")
                .replace("_grad_rho", "")
                .replace("_grad_z", "")
            )

    # TODO: check if it works correctly in all cases!!!
    # CALIBRATION: remove spatially varying parameters that are not inputs of interest, unless the gradient is different to 1
    if len(all_spatial_inputs):
        tmp = {}
        for _, key in enumerate(par_all_dct):
            if (
                key.replace("_theta", "")
                .replace("_grad_rho", "")
                .replace("_grad_z", "")
                in all_spatial_inputs
            ):
                tmp[key] = par_all_dct[key]
            elif "_theta" in key:
                key_tmp = key.replace("_theta", "")
                if (
                    par_all_dct[key_tmp + "_grad_rho"] == 1
                    or par_all_dct[key_tmp + "_grad_z"] == 1
                ):
                    tmp[key_tmp] = par_all_dct[key]
                else:
                    tmp[key] = par_all_dct[key]
                    tmp[key_tmp + "_grad_rho"] = par_all_dct[key_tmp + "_grad_rho"]
                    tmp[key_tmp + "_grad_z"] = par_all_dct[key_tmp + "_grad_z"]
            elif not any(["_theta" in key, "_grad_rho" in key, "_grad_z" in key]):
                tmp[key] = par_all_dct[key]
        par_all_dct = tmp

    # find fixed parameters
    par_fixed_dct = {}
    for key in sorted(list(set(par_all_dct.keys()) - set(par_dct.keys()))):
        par_fixed_dct[key] = par_all_dct[key]

    # for each sample, generete a .json file with all simulation parameters (including fixed ones)
    print("\nGenerating json file for each sample...")
    for i in range(N_samples):
        # build current sample output directory
        SIMID = outputfolder + "/" + str(i)
        if not os.path.exists(SIMID):
            os.mkdir(SIMID)

        # write parameters in json file
        pfile = SIMID + "/param"
        IO.write_json(pfile, par_fixed_dct)

    return


def assemble_simulation_param(
    paramfolder,
    paramfilename,
    labelfolder,
    labelfilename,
    simulfolder,
    simulfilename,
    outputfolder,
):
    """Save both fixed parameters and parameters of interest.

    Parameters
    ----------
    paramfolder: str
        Path to the electrophysiology parameter datafile.
    paramfilename: str
        Name of .json file with names and reference values of all simulation parameters (without extension).
    labelfolder: str
        Path to the input parameter datafile.
    labelfilename: str
        Name of .txt file with lables of parameters of interest (without extension).
    simulfolder: str
        Path to the input values datafile.
    simulfilename: str
        Name of .txt file with values of parameters of interest (without extension).
    outputfolder: str
        Path to the directory where the sample datafile(s) will be saved.
    """

    # load dictionary with all simulation parameters (including those that are NOT considered as inputs)
    paramfile = paramfolder + "/" + paramfilename
    par_all_dct = IO.read_json(paramfile)

    # find parameters of interest as given in labelfile
    labelfile = labelfolder + "/" + labelfilename
    par_dct = {}

    ## TODO: check if it works correctley
    for key in IO.read_labels(labelfile):
        if key in list(par_dct.keys()):
            raise Exception(
                "{} can't be both scalar and spatially varying. Remove it from one of the label files.".format(
                    key
                )
            )
        par_dct[key] = par_all_dct[key]
    par_keys = list(par_dct.keys())

    # find fixed parameters
    par_fixed_dct = {}
    for key in sorted(list(set(par_all_dct.keys()) - set(par_keys))):
        par_fixed_dct[key] = par_all_dct[key]

    # read simulation parameter values from text file
    simulfile = simulfolder + "/" + simulfilename
    X_simul = IO.read_float_text(simulfile)

    # for each sample, generete a .json file with all simulation parameters (including fixed ones)
    print("\nGenerating json file for each sample...")
    for i in range(np.size(X_simul, 0)):
        # build current sample output directory
        SIMID = outputfolder + "/" + str(i)
        if not os.path.exists(SIMID):
            os.mkdir(SIMID)

        # write fixed parameters in json file
        pfile = SIMID + "/param"
        IO.write_json(pfile, par_fixed_dct)

        # get current input values
        p_temp = {}
        for j, key in enumerate(par_keys):
            p_temp[key] = X_simul[i, j]

        # append current input parameters to json file
        IO.append_json(pfile, p_temp)

    return


def sample_scalar_param(
    paramfolder,
    paramfilename,
    labelfolder,
    labelfilename,
    outputfolder,
    N_samples=500,
    nfold_range=[0.5, 2],
    usr_range_json="",
):
    """Sample scalar input parameters for global sensitivity analysis.

    Parameters
    ----------
    paramfolder: str
        Path to the electrophysiology parameter datafile.
    paramfilename: str
        Name of .json file with names and reference values of all simulation parameters (without extension).
    labelfolder: str
        Path to the input parameter datafile.
    labelfilename: str
        Name of .txt file with lables of parameters of interest, i.e., GSA inputs, (without extension).
    outputfolder: str
        Path to the directory where the sample datafile(s) will be saved.
    N_samples: int
        Number of samples. Default: 500.
    nfold_range: list
        Range of n-fold variation of each input with respect to its reference value. Default: [0.5, 2].
    usr_range_json: str
        Path to json file (without extension) with user-defined range of variation, if any. Otherwise use empty string. Default: "".
        Note: The json file must be structured as {"name_par":[min,max], ...}.
    """

    # set path to parameters .json file and to input labels .txt file
    paramfile = paramfolder + "/" + paramfilename
    labelfile = labelfolder + "/" + labelfilename

    # check folder structure
    print("\nChecking folder structure...")
    to_check = [paramfile + ".json", labelfile + ".txt"]
    for f in to_check:
        if not os.path.exists(f):
            raise Exception("Cannot find {}. Check path to file.".format(f))
        else:
            print(f + " found.")

    # check usr_range
    if usr_range_json:
        usr_range_dct = IO.read_json(usr_range_json)
        usr_range = [list(usr_range_dct.keys()), list(usr_range_dct.values())]
        if not len(usr_range) == 2 or not len(usr_range[0]) == len(usr_range[1]):
            raise Exception("User-defined ranges not valid.")
    else:
        usr_range = []

    # set parameter ranges and print information
    print(
        "\nUsing "
        + str(nfold_range[0])
        + "-fold to "
        + str(nfold_range[1])
        + "-fold range of variation..."
    )
    if len(usr_range) > 0:
        print("...except for {}".format(usr_range[0]))

    # load dictionary with all simulation parameters (including those that are NOT considered as inputs)
    par_all_dct = IO.read_json(paramfile)

    # find parameters of interest as given in labelfile
    par_dct = {}
    for key in IO.read_labels(labelfile):
        par_dct[key] = par_all_dct[key]
    par_keys = list(par_dct.keys())
    par_values = list(par_dct.values())

    # matrix of parameter ranges
    par_range = np.zeros((len(par_keys), 2))
    for i, val in enumerate(par_values):
        if len(usr_range) > 0:
            if par_keys[i] in usr_range[0]:
                for j, pr in enumerate(usr_range[0]):
                    if pr == par_keys[i]:
                        par_range[i, :] = np.array(
                            [usr_range[1][j][0], usr_range[1][j][1]]
                        )
            else:
                par_range[i, :] = np.sort(
                    np.array([val * nfold_range[0], val * nfold_range[1]])
                )
        else:
            par_range[i, :] = np.sort(
                np.array([val * nfold_range[0], val * nfold_range[1]])
            )

        # impose constrains on specific parameters
        bounds = []
        if (
            "Zmax_SE" in par_keys[i]
            or "_Z" in par_keys[i]
            or "fractionSS" in par_keys[i]
        ):
            bounds = [0, 1]
        elif "_PHI" in par_keys[i]:
            bounds = [-np.pi, np.pi]
        if bounds:
            if par_range[i, 0] < bounds[0] or par_range[i, 1] > bounds[1]:
                print(par_keys[i] + " limited within {}".format(bounds))
                par_range[i, 0] = max(par_range[i, 0], bounds[0])
                par_range[i, 1] = min(par_range[i, 1], bounds[1])

    # print parameter ranges
    print("\nFinal parameter ranges:")
    for i in range(len(par_values)):
        print("{} \t {}".format(par_keys[i], par_range[i, :]))

    # matrix of samples within the parameter ranges
    print("\nGenerating {} samples...".format(N_samples))
    H = lhd(par_range, N_samples)

    # for each sample, generete a .json file with all simulation parameters (including fixed ones)
    print("\nGenerating json file for each sample...")
    for i in range(N_samples):
        # build current sample output directory
        SIMID = outputfolder + "/" + str(i)
        if not os.path.exists(SIMID):
            raise Exception("Run save_fixed_parameters first.")

        p_temp = {}
        for j, key in enumerate(par_keys):
            p_temp[key] = H[i, j]

        # write parameters in json file
        pfile = SIMID + "/param"
        if not os.path.exists(pfile + ".json"):
            IO.write_json(pfile, p_temp)
        else:
            IO.append_json(pfile, p_temp)

    return




def sample_constraint_gradient(
    var_name,
    range_val,
    range_grad_rho,
    range_grad_z,
    N_samples,
    batch_size=100,
    max_iter=150,
):
    """
    Parameters
    ----------
    var_name: str
        Name of spatially-varying parameter.
    range_val: list
        Range of variation for spatially-varying parameter baseline value. To be structured as [min,max].
    range_grad_rho: list
        Range of variation for spatially-varying parameter transmural gradient. To be structured as [min,max].
    range_grad_z: list
        Range of variation for spatially-varying parameter longitudinal gradient. To be structured as [min,max].
    N_samples: int
        Number of samples.
    batch_size: int
        Size of current batch. Default: 100.
    max_iter:
        Maximum number of iterations allowed to find a sufficient number of accepted samples. Default: 150.
    """

    # initialisation
    H_tmp = np.zeros((N_samples + batch_size, 3))
    n_accept = 0
    n_iter = 0

    # loop until the desired number of samples has been reached
    while n_accept < N_samples:
        # update iterations
        n_iter = n_iter + 1

        # sample from discrete uniform distributions
        val = range_val[0] + (range_val[1] - range_val[0]) * np.random.sample(
            batch_size
        )
        grad_rho = range_grad_rho[0] + (
            range_grad_rho[1] - range_grad_rho[0]
        ) * np.random.sample(batch_size)
        grad_z = range_grad_z[0] + (
            range_grad_z[1] - range_grad_z[0]
        ) * np.random.sample(batch_size)

        # check where constraints are fulfilled
        evaluate_constraints = (
            (val * grad_rho >= range_val[0])
            & (val * grad_rho <= range_val[1])
            & (val * grad_z >= range_val[0])
            & (val * grad_z <= range_val[1])
            & (val * grad_rho * grad_z >= range_val[0])
            & (val * grad_rho * grad_z <= range_val[1])
        )

        # number of accepted combinations in current batch
        n_accept_update = n_accept + sum(evaluate_constraints)

        # update accepted combinations
        H_tmp[n_accept:n_accept_update] = np.stack(
            (
                val[evaluate_constraints],
                grad_rho[evaluate_constraints],
                grad_z[evaluate_constraints],
            )
        ).T
        n_accept = n_accept_update

        if n_iter > max_iter:
            raise Exception(
                "Reached maximum number of iterations for {}. Cannot find enough samples satisfying the constraints.".format(
                    var_name
                )
            )

    return H_tmp[:N_samples]


def sample_gradient_param(
    paramfolder,
    paramfilename,
    labelfolder,
    labelfilename,
    imfilename,
    outputfolder,
    N_samples=500,
    nfold_range=[0.5, 2],
    usr_range_json="",
    usr_grad_json="",
):
    """Sample spatially varying (apico-basally and transmurally) input parameters for global sensitivity analysis.

    Parameters
    ----------
    paramfolder: str
        Path to the electrophysiology parameter datafile.
    paramfilename: str
        Name of .json file with names and reference values of all simulation parameters (without extension).
    labelfolder: str
        Path to the input parameter datafile.
    labelfilename: str
        Name of .txt file with lables of parameters of interest, i.e., GSA inputs, (without extension).
    imfilename: str
        Name of .txt file with labels of all ionic model parameters (without extension).
    outputfolder: str
        Path to the directory where the sample datafile(s) will be saved.
    N_samples: int
        Number of samples. Default: 500.
    nfold_range: list
        Range of n-fold variation of each input with respect to its reference value. Default: [0.5, 2].
    usr_range_json: str
        Path to json file (without extension) with user-defined range of variation, if any. Otherwise use empty string. Default: "".
        Note: The json file must be structured as {"name_par":[min,max], ...}.
    usr_grad_json: str
        Path to json file (without extension) with user-defined range of variation for the transmural (endo-to-epi) and the longitudinal (base-to-apex) gradient., if any. Otherwise use empty string. Default: "".
        Note: The json file must be structured as {name_par:[[min_grad_rho,max_grad_rho],[min_grad_z,max_grad_z], ...}.
    """

    # Define
    #   1- apex_endo/apex_epi = base_endo/base_epi =: grad_rho
    #   2- base_epi/apex_epi = base_endo/apex_endo =: grad_z
    # such that
    #   1- apex_epi =: theta
    #   2- apex_endo = theta * grad_rho
    #   3- base_epi = theta * grad_z
    #   4- base_endo = theta * grad_rho * grad_z

    # NOTE: since theta > 0 (always),
    #   - grad_rho > 1 ==> ENDO > EPI
    #   - grad_z > 1   ==> BASE > APEX

    # set path to parameters .json file and to input labels .txt file
    paramfile = paramfolder + "/" + paramfilename
    labelfile = labelfolder + "/" + labelfilename
    imfile = labelfolder + "/" + imfilename

    # check folder structure
    print("\nChecking folder structure...")
    to_check = [paramfile + ".json", labelfile + ".txt", imfile + ".txt"]
    for f in to_check:
        if not os.path.exists(f):
            raise Exception("Cannot find {}. Check path to file.".format(f))
        else:
            print(f + " found.")

    if not IO.read_text(labelfile):
        print(
            "NOTE: {} is empty. There are no spatially varying parameters.".format(
                labelfile
            )
        )
    else:
        # check usr_range
        if usr_range_json:
            usr_range_dct = IO.read_json(usr_range_json)
            usr_range = [list(usr_range_dct.keys()), list(usr_range_dct.values())]
            if not len(usr_range) == 2 or not len(usr_range[0]) == len(usr_range[1]):
                raise Exception("User-defined ranges not valid.")
        else:
            usr_range = []

        # check usr_grad
        if usr_grad_json:
            usr_grad_dct = IO.read_json(usr_grad_json)
            grad_var = list(usr_grad_dct.keys())  # variable
            grad_dir = list(usr_grad_dct.values())  # direction
            list_usr_grad_dct_values = []
            for i, var in enumerate(grad_var):
                list_usr_grad_dct_values.append(
                    [grad_dir[i]["grad_rho"], grad_dir[i]["grad_z"]]
                )
                usr_grad_dct[var] = [grad_dir[i]["grad_rho"], grad_dir[i]["grad_z"]]
            usr_grad = [grad_var, list_usr_grad_dct_values]
            if not len(usr_grad) == 2 or not len(usr_grad[0]) == len(usr_grad[1]):
                raise Exception("User-defined ranges not valid.")
        else:
            usr_grad = []

        # set parameter ranges and print information
        print(
            "\nUsing "
            + str(nfold_range[0])
            + "-fold to "
            + str(nfold_range[1])
            + "-fold range of variation..."
        )
        if len(usr_range) > 0:
            print("...except for {}".format(usr_range[0]))

        if len(usr_grad) > 0:
            print(
                "Restrict transmural and longitudinal gradients for {}\n".format(
                    usr_grad[0]
                )
            )

        # load dictionary with all simulation parameters (including those that are NOT considered as inputs)
        par_all_dct = IO.read_json(paramfile)

        # find parameters of interest as given in labelfile
        par_dct = {}
        for key in IO.read_labels(labelfile):
            if not key == "vel_f" and not key in list(IO.read_text(imfile)):
                raise Exception(
                    "Only vel_f or ionic model parameters can be varied spatially.\nPlease remove {} from {}.".format(
                        key, labelfile
                    )
                )
            par_dct[key] = par_all_dct[key]
        par_keys = list(par_dct.keys())
        par_values = list(par_dct.values())

        # matrix of parameter ranges
        par_range = np.zeros((len(par_keys), 2))
        for i, val in enumerate(par_values):
            if len(usr_range) > 0:
                if par_keys[i] in usr_range[0]:
                    for j, pr in enumerate(usr_range[0]):
                        if pr == par_keys[i]:
                            par_range[i, :] = np.array(
                                [usr_range[1][j][0], usr_range[1][j][1]]
                            )
                else:
                    par_range[i, :] = np.sort(
                        np.array([val * nfold_range[0], val * nfold_range[1]])
                    )
                    usr_range_dct[par_keys[i]] = [
                        val * nfold_range[0],
                        val * nfold_range[1],
                    ]
            else:
                par_range[i, :] = np.sort(
                    np.array([val * nfold_range[0], val * nfold_range[1]])
                )
                usr_range_dct[par_keys[i]] = [
                    val * nfold_range[0],
                    val * nfold_range[1],
                ]

            # impose constrains on specific parameters
            bounds = []
            if "fractionSS" in par_keys[i]:
                bounds = [0, 1]
            if bounds:
                if par_range[i, 0] < bounds[0] or par_range[i, 1] > bounds[1]:
                    print(par_keys[i] + " limited within {}".format(bounds))
                    par_range[i, 0] = max(par_range[i, 0], bounds[0])
                    par_range[i, 1] = min(par_range[i, 1], bounds[1])

        # print parameter ranges
        print("\nFinal parameter ranges:")
        for i in range(len(par_values)):
            print("{} \t {}".format(par_keys[i], par_range[i, :]))

        # matrix of samples within the parameter ranges
        print("\nGenerating {} samples...".format(N_samples))

        n_biv_ch = 1
        biv_ch = ["_"]

        H = np.zeros((N_samples, len(par_keys), n_biv_ch, 3))
        for k in range(n_biv_ch):
            H[:, :, k, 0] = lhd(par_range, N_samples)  # theta
            H[:, :, k, 1] = lhd(par_range, N_samples) / H[:, :, k, 0]  # grad_rho
            H[:, :, k, 2] = lhd(par_range, N_samples) / H[:, :, k, 0]  # grad_z

        # TODO: check that it works correctly for NON usr-defined parameter ranges
        # constrain ranges for gradients
        if len(usr_grad) > 0:
            for j, key in enumerate(par_keys):
                for k, ch in enumerate(biv_ch):
                    if key in usr_grad[0]:
                        H[:, j, k, :] = sample_constraint_gradient(
                            key,
                            usr_range_dct[key],
                            usr_grad_dct[key][0],
                            usr_grad_dct[key][1],
                            N_samples,
                        )
                    else:
                        H[:, j, k, :] = sample_constraint_gradient(
                            key,
                            usr_range_dct[key],
                            nfold_range,
                            nfold_range,
                            N_samples,
                        )
        else:
            for j, key in enumerate(par_keys):
                for k, ch in enumerate(biv_ch):
                    H[:, j, k, :] = sample_constraint_gradient(
                        key,
                        usr_range_dct[key],
                        nfold_range,
                        nfold_range,
                        N_samples,
                    )

        # for each sample, generete a .json file with all simulation parameters (including fixed ones)
        print("\nGenerating json file for each sample...")
        for i in range(N_samples):
            # build current sample output directory
            SIMID = outputfolder + "/" + str(i)
            if not os.path.exists(SIMID):
                raise Exception("Run save_fixed_parameters first.")

            p_temp = {}
            for j, key in enumerate(par_keys):
                for k, ch in enumerate(biv_ch):
                    p_temp[key + ch + "theta"] = H[i, j, k, 0]
                    p_temp[key + ch + "grad_rho"] = H[i, j, k, 1]
                    p_temp[key + ch + "grad_z"] = H[i, j, k, 2]

            # write parameters in json file
            pfile = SIMID + "/param"
            if not os.path.exists(pfile + ".json"):
                IO.write_json(pfile, p_temp)
            else:
                IO.append_json(pfile, p_temp)

    return


def sample_fitting(
    paramfolder,
    paramfilename,
    labelfolder,
    labelfilename,
    outputfolder,
    nfold_range=[0.5, 2],
    usr_range_dct={},
    usr_range_json="",
    normpar_list=[
        "alpha_endo",
        "alpha_epi",
        "Zmax_SE",
        "LV_af_PHI",
        "LV_af_Z",
        "LV_af_t",
        "LV_pf_PHI",
        "LV_pf_Z",
        "LV_pf_t",
        "LV_sf_PHI",
        "LV_sf_Z",
        "LV_sf_t",
        "RV_mod_PHI",
        "RV_mod_Z",
        "RV_mod_t",
        "RV_sf_PHI",
        "RV_sf_Z",
        "RV_sf_t",
        "aniso_ratio",
        "CV_SE_ratio",
        "vel_f_theta",
        "vel_f_grad_rho",
        "vel_f_grad_z",
    ],
):
    """Define input parameters for polynomial fitting.
    Sample matrix H is defined as follows,
        for p = 1:Nparam
            - H[2p,p] = min(par[p])
            - H[2p+1,p] = max(par[p])
            for k != p:
                if par[p] is normally distributed:
                    - H[2p,k] = H[2p+1,k] = mean([min(par[k]), max(par[k])])
                else (i.e., if par[p] is log-normally distributed):
                    - H[2p,k] = H[2p+1,k] = exp(mean([ln(min(par[k])), ln(max(par[k]))]))
    Parameters
    ----------
    paramfolder: str
        Path to the electrophysiology parameter datafile.
    paramfilename: str
        Name of .json file with names and reference values of all simulation parameters (without extension).
    labelfolder: str
        Path to the input parameter datafile.
    labelfilename: str
        Name of .txt file with lables of all parameters of interest (without extension).
        For each spatially varying parameter, define "*_theta", "*_grad_rho" and "*_grad_z" (for each chamber, if needed).
    outputfolder: str
        Path to the directory where the sample datafile(s) will be saved.
    usr_range_json: str
        Path to json file (without extension) with user-defined range of variation, if any. Otherwise use empty string. Default: "".
        Note: The json file must be structured as {"name_par":[min,max], ...}.
    normpar_list: list
        List of normally distributed parameters. All other parameters are assumed to be log-normally distributed. Default: fibres angles, HPS, anisotropy ratios.
    """

    # set path to parameters .json file and to input labels .txt file
    paramfile = paramfolder + "/" + paramfilename
    labelfile = labelfolder + "/" + labelfilename

    # check folder structure
    print("\nChecking folder structure...")
    to_check = [paramfile + ".json", labelfile + ".txt"]
    for f in to_check:
        if not os.path.exists(f):
            raise Exception("Cannot find {}. Check path to file.".format(f))
        else:
            print(f + " found.")

    # check usr_range
    if usr_range_json:
        usr_range_dct = IO.read_json(usr_range_json)
        usr_range = [list(usr_range_dct.keys()), list(usr_range_dct.values())]
        if not len(usr_range) == 2 or not len(usr_range[0]) == len(usr_range[1]):
            raise Exception("User-defined ranges not valid.")
    else:
        usr_range = []

    # set parameter ranges and print information
    print(
        "\nUsing "
        + str(nfold_range[0])
        + "-fold to "
        + str(nfold_range[1])
        + "-fold range of variation..."
    )
    if len(usr_range) > 0:
        print("...except for {}".format(usr_range[0]))

    # load dictionary with all simulation parameters (including those that are NOT considered as inputs)
    par_all_dct = IO.read_json(paramfile)

    # find parameters of interest as given in labelfile
    par_dct = {}
    for key in IO.read_labels(labelfile):
        par_dct[key] = par_all_dct[key]
    par_keys = list(par_dct.keys())
    par_values = list(par_dct.values())

    # matrix of parameter ranges
    par_range = np.zeros((len(par_keys), 2))
    for i, val in enumerate(par_values):
        if len(usr_range) > 0:
            if par_keys[i] in usr_range[0]:
                for j, pr in enumerate(usr_range[0]):
                    if pr == par_keys[i]:
                        par_range[i, :] = np.array(
                            [usr_range[1][j][0], usr_range[1][j][1]]
                        )
            else:
                par_range[i, :] = np.sort(
                    np.array([val * nfold_range[0], val * nfold_range[1]])
                )
        else:
            par_range[i, :] = np.sort(
                np.array([val * nfold_range[0], val * nfold_range[1]])
            )

        # impose constrains on specific parameters
        bounds = []
        if (
            "Zmax_SE" in par_keys[i]
            or "_Z" in par_keys[i]
            or "fractionSS" in par_keys[i]
        ):
            bounds = [0, 1]
        elif "_PHI" in par_keys[i]:
            bounds = [-np.pi, np.pi]
        if bounds:
            if par_range[i, 0] < bounds[0] or par_range[i, 1] > bounds[1]:
                print(par_keys[i] + " limited within {}".format(bounds))
                par_range[i, 0] = max(par_range[i, 0], bounds[0])
                par_range[i, 1] = min(par_range[i, 1], bounds[1])

    # print parameter ranges
    print("\nFinal parameter ranges:")
    for i in range(len(par_values)):
        print("{} \t {}".format(par_keys[i], par_range[i, :]))

    # find extended-CoM
    extCoM = np.zeros((len(par_keys), 1))
    for k in range(len(par_keys)):
        if par_keys[k] in normpar_list:
            extCoM[k] = (par_range[k, 1] + par_range[k, 0]) / 2
        else:
            extCoM[k] = np.exp((np.log(par_range[k, 1]) + np.log(par_range[k, 0])) / 2)

    # matrix of samples within the parameter ranges
    print("\nGenerating {} samples...".format(2 * len(par_keys)))
    H_labels = []
    H = np.zeros((2 * len(par_keys), len(par_keys)))
    for k in range(len(par_keys)):
        H_labels.append(par_keys[k] + "\n")
        for p in range(len(par_keys)):
            if p == k:
                H[2 * k, p] = par_range[p, 0]
                H[2 * k + 1, p] = par_range[p, 1]
            else:
                H[2 * k, p] = extCoM[p]
                H[2 * k + 1, p] = extCoM[p]

    # save sample matrix
    IO.write_array_text(outputfolder + "/ThetaArray", H)

    # save xlabels used to build ThetaArray
    IO.write_string_text(outputfolder + "/ThetaArray_labels", H_labels)

    # for each sample, generete a .json file with all simulation parameters (including fixed ones)
    print("\nGenerating json file for each sample...")
    for i in range(2 * len(par_keys)):
        # build current sample output directory
        SIMID = outputfolder + "/" + str(i)
        if not os.path.exists(SIMID):
            raise Exception("Run save_fixed_parameters first.")

        p_temp = {}
        for j, key in enumerate(par_keys):
            p_temp[key] = H[i, j]

        # write parameters in json file
        pfile = SIMID + "/param"
        if not os.path.exists(pfile + ".json"):
            IO.write_json(pfile, p_temp)
        else:
            IO.append_json(pfile, p_temp)

    return


def sample_constraint_calibration(
    var_name,
    target_range_val,
    range_val,
    range_grad_rho,
    range_grad_z,
    N_samples,
    batch_size=100,
    max_iter=150,
):
    """
    Parameters
    ----------
    var_name: str
        Name of spatially-varying parameter.
    target_range_val: list
        Range of variation for spatially-varying parameter baseline value. To be structured as [min,max].
    range_val: list
        Range of variation for spatially-varying parameter baseline value. To be structured as [min,max].
    range_grad_rho: list
        Range of variation for spatially-varying parameter transmural gradient. To be structured as [min,max].
    range_grad_z: list
        Range of variation for spatially-varying parameter longitudinal gradient. To be structured as [min,max].
    N_samples: int
        Number of samples.
    batch_size: int
        Size of current batch. Default: 100.
    max_iter:
        Maximum number of iterations allowed to find a sufficient number of accepted samples. Default: 150.
    """

    # initialisation
    H_tmp = np.zeros((N_samples + batch_size, 3))
    n_accept = 0
    n_iter = 0

    # loop until the desired number of samples has been reached
    while n_accept < N_samples:
        # update iterations
        n_iter = n_iter + 1

        # sample from discrete uniform distributions
        val = range_val[0] + (range_val[1] - range_val[0]) * np.random.sample(
            batch_size
        )
        grad_rho = range_grad_rho[0] + (
            range_grad_rho[1] - range_grad_rho[0]
        ) * np.random.sample(batch_size)
        grad_z = range_grad_z[0] + (
            range_grad_z[1] - range_grad_z[0]
        ) * np.random.sample(batch_size)

        # check where constraints are fulfilled
        evaluate_constraints = (
            (val * grad_rho >= target_range_val[0])
            & (val * grad_rho <= target_range_val[1])
            & (val * grad_z >= target_range_val[0])
            & (val * grad_z <= target_range_val[1])
            & (val * grad_rho * grad_z >= target_range_val[0])
            & (val * grad_rho * grad_z <= target_range_val[1])
        )

        # number of accepted combinations in current batch
        n_accept_update = n_accept + sum(evaluate_constraints)

        # update accepted combinations
        H_tmp[n_accept:n_accept_update] = np.stack(
            (
                val[evaluate_constraints],
                grad_rho[evaluate_constraints],
                grad_z[evaluate_constraints],
            )
        ).T
        n_accept = n_accept_update

        if n_iter > max_iter:
            raise Exception(
                "Reached maximum number of iterations for {}. Cannot find enough samples satisfying the constraints.".format(
                    var_name
                )
            )

    return H_tmp[:N_samples]


def sample_calibration_param(
    paramfolder,
    paramfilename,
    labelfolder,
    labelfilename,
    outputfolder,
    N_samples=500,
    nfold_range=[0.5, 2],
    usr_range_json="",
):
    """Sample (scalar) input parameters for calibration. Gradients are sampled independently.

    Parameters
    ----------
    paramfolder: str
        Path to the electrophysiology parameter datafile.
    paramfilename: str
        Name of .json file with names and reference values of all simulation parameters (without extension).
    labelfolder: str
        Path to the input parameter datafile.
    labelfilename: str
        Name of .txt file with lables of parameters of interest, i.e., GSA inputs, (without extension).
    outputfolder: str
        Path to the directory where the sample datafile(s) will be saved.
    N_samples: int
        Number of samples. Default: 500.
    nfold_range: list
        Range of n-fold variation of each input with respect to its reference value. Default: [0.5, 2].
    usr_range_json: str
        Path to json file (without extension) with user-defined range of variation, if any. Otherwise use empty string. Default: "".
        Note: The json file must be structured as {"name_par":[min,max], ...}.
    """

    # set path to parameters .json file and to input labels .txt file
    paramfile = paramfolder + "/" + paramfilename
    labelfile = labelfolder + "/" + labelfilename

    # check folder structure
    print("\nChecking folder structure...")
    to_check = [paramfile + ".json", labelfile + ".txt"]
    for f in to_check:
        if not os.path.exists(f):
            raise Exception("Cannot find {}. Check path to file.".format(f))
        else:
            print(f + " found.")

    # check usr_range
    if usr_range_json:
        usr_range_dct = IO.read_json(usr_range_json)
        usr_range = [list(usr_range_dct.keys()), list(usr_range_dct.values())]
        if not len(usr_range) == 2 or not len(usr_range[0]) == len(usr_range[1]):
            raise Exception("User-defined ranges not valid.")
    else:
        usr_range = []

    # set parameter ranges and print information
    print(
        "\nUsing "
        + str(nfold_range[0])
        + "-fold to "
        + str(nfold_range[1])
        + "-fold range of variation..."
    )
    if len(usr_range) > 0:
        print("...except for {}".format(usr_range[0]))

    # load dictionary with all simulation parameters (including those that are NOT considered as inputs)
    par_all_dct = IO.read_json(paramfile)

    # find parameters of interest as given in labelfile
    par_dct = {}
    for key in IO.read_labels(labelfile):
        par_dct[key] = par_all_dct[key]
    par_keys = list(par_dct.keys())
    par_values = list(par_dct.values())

    # matrix of parameter ranges
    par_range = np.zeros((len(par_keys), 2))
    par_range_ind_dct = {}
    for i, val in enumerate(par_values):
        if len(usr_range) > 0:
            if par_keys[i] in usr_range[0]:
                for j, pr in enumerate(usr_range[0]):
                    if pr == par_keys[i]:
                        par_range[i, :] = np.array(
                            [usr_range[1][j][0], usr_range[1][j][1]]
                        )
            else:
                par_range[i, :] = np.sort(
                    np.array([val * nfold_range[0], val * nfold_range[1]])
                )
        else:
            par_range[i, :] = np.sort(
                np.array([val * nfold_range[0], val * nfold_range[1]])
            )

        # impose constrains on specific parameters
        bounds = []
        if (
            "Zmax_SE" in par_keys[i]
            or "_Z" in par_keys[i]
            or "fractionSS" in par_keys[i]
        ):
            bounds = [0, 1]
        elif "_PHI" in par_keys[i]:
            bounds = [-np.pi, np.pi]
        if bounds:
            if par_range[i, 0] < bounds[0] or par_range[i, 1] > bounds[1]:
                print(par_keys[i] + " limited within {}".format(bounds))
                par_range[i, 0] = max(par_range[i, 0], bounds[0])
                par_range[i, 1] = min(par_range[i, 1], bounds[1])

        # save also dictionary of parameter ranges
        par_range_ind_dct[par_keys[i]] = i  # par_range[i, :]

    # print parameter ranges
    print("\nFinal parameter ranges:")
    for i in range(len(par_values)):
        print("{} \t {}".format(par_keys[i], par_range[i, :]))

    # matrix of samples within the parameter ranges
    print("\nGenerating {} samples...".format(N_samples))
    H = lhd(par_range, N_samples)

    # find parameters such that at least two among theta, grad_rho and grad_z are inputs of interest
    spatial_inputs = []
    tmp = []
    for key in par_keys:
        if "_theta" in key or "_grad_rho" in key or "_grad_z" in key:
            tmp.append(
                key.replace("_theta", "")
                .replace("_grad_rho", "")
                .replace("_grad_z", "")
            )
    for key in list(set(tmp)):
        if [
            key + "_theta" in par_keys,
            key + "_grad_rho" in par_keys,
            key + "_grad_z" in par_keys,
        ].count(True) > 1:
            spatial_inputs.append(key)

    # if only one among theta, grad_rho and grad_z (for each input and each chamber) is an input of interest, there is no need to re-sample
    if len(spatial_inputs):
        # matrix of theta parameter ranges
        par_theta_range_dct = {}
        for key in list(par_all_dct.keys()):
            if "theta" in key:
                if len(usr_range) > 0:
                    if key in usr_range[0]:
                        for j, pr in enumerate(usr_range[0]):
                            if pr == key:
                                par_theta_range_dct[key] = np.array(
                                    [usr_range[1][j][0], usr_range[1][j][1]]
                                )
                    else:
                        par_theta_range_dct[key] = np.sort(
                            np.array(
                                [
                                    par_all_dct[key] * nfold_range[0],
                                    par_all_dct[key] * nfold_range[1],
                                ]
                            )
                        )
                else:
                    par_theta_range_dct[key] = np.sort(
                        np.array(
                            [
                                par_all_dct[key] * nfold_range[0],
                                par_all_dct[key] * nfold_range[1],
                            ]
                        )
                    )

        # dictionary of spatially varying parameters
        par_spatial_dct = {}
        par_spatial_lo = []  # left-out component of spatially varying parameter, if any
        for key in spatial_inputs:
            tmp = [
                key + "_theta" in par_keys,
                key + "_grad_rho" in par_keys,
                key + "_grad_z" in par_keys,
            ]
            if tmp.count(True) == 2:
                if tmp[0]:
                    par_spatial_dct[key + "_theta"] = par_all_dct[key + "_theta"]
                else:
                    par_spatial_lo.append(key + "_theta")
                if tmp[1]:
                    par_spatial_dct[key + "_grad_rho"] = par_all_dct[key + "_grad_rho"]
                else:
                    par_spatial_lo.append(key + "_grad_rho")
                if tmp[2]:
                    par_spatial_dct[key + "_grad_z"] = par_all_dct[key + "_grad_z"]
                else:
                    par_spatial_lo.append(key + "_grad_z")
            else:
                par_spatial_dct[key + "_theta"] = par_all_dct[key + "_theta"]
                par_spatial_dct[key + "_grad_rho"] = par_all_dct[key + "_grad_rho"]
                par_spatial_dct[key + "_grad_z"] = par_all_dct[key + "_grad_z"]
        par_spatial_values = list(par_spatial_dct.values())
        par_spatial_keys = list(par_spatial_dct.keys())

        # matrix of spatially varying parameter ranges
        par_spatial_range_dct = {}
        for i, val in enumerate(par_spatial_values):
            if len(usr_range) > 0:
                if par_spatial_keys[i] in usr_range[0]:
                    for j, pr in enumerate(usr_range[0]):
                        if pr == par_spatial_keys[i]:
                            par_spatial_range_dct[par_spatial_keys[i]] = np.array(
                                [usr_range[1][j][0], usr_range[1][j][1]]
                            )
                else:
                    par_spatial_range_dct[par_spatial_keys[i]] = np.sort(
                        np.array([val * nfold_range[0], val * nfold_range[1]])
                    )
            else:
                par_spatial_range_dct[par_spatial_keys[i]] = np.sort(
                    np.array([val * nfold_range[0], val * nfold_range[1]])
                )
        # for the left-out componet, define a "fake" range of variantion
        for key in par_spatial_lo:
            par_spatial_range_dct[key] = np.array([par_all_dct[key], par_all_dct[key]])

        # print parameter ranges
        print("\nSpatial parameter ranges:")
        for _, key in enumerate(par_spatial_range_dct):
            print("{} \t {}".format(key, par_spatial_range_dct[key]))

        # resample if necessary
        for sp_key in spatial_inputs:
            idx = []
            idx.append(
                par_range_ind_dct[sp_key + "_theta"]
                if sp_key + "_theta" in list(par_range_ind_dct.keys())
                else None
            )
            idx.append(
                par_range_ind_dct[sp_key + "_grad_rho"]
                if sp_key + "_grad_rho" in list(par_range_ind_dct.keys())
                else None
            )
            idx.append(
                par_range_ind_dct[sp_key + "_grad_z"]
                if sp_key + "_grad_z" in list(par_range_ind_dct.keys())
                else None
            )

            H_tmp = sample_constraint_calibration(
                sp_key,
                par_theta_range_dct[sp_key + "_theta"],
                par_spatial_range_dct[sp_key + "_theta"],
                par_spatial_range_dct[sp_key + "_grad_rho"],
                par_spatial_range_dct[sp_key + "_grad_z"],
                N_samples,
                batch_size=1000,
                max_iter=1500,
            )

            for i in range(3):
                if not idx[i] == None:
                    H[:, idx[i]] = H_tmp[:, i]

    # for each sample, generete a .json file with all simulation parameters (including fixed ones)
    print("\nGenerating json file for each sample...")
    for i in range(N_samples):
        # build current sample output directory
        SIMID = outputfolder + "/" + str(i)
        if not os.path.exists(SIMID):
            raise Exception("Run save_fixed_parameters first.")

        p_temp = {}
        for j, key in enumerate(par_keys):
            p_temp[key] = H[i, j]

        # write parameters in json file
        pfile = SIMID + "/param"
        if not os.path.exists(pfile + ".json"):
            IO.write_json(pfile, p_temp)
        else:
            IO.append_json(pfile, p_temp)

    return
