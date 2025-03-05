#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# __main__.py
#
# 04 de gen. 2025  <adria@molevol-OptiPlex-9020>

help_msg = """
Describe here the goal, input, output, etc. of the script
as a multi-line block of text.
"""

import sys
import pprint
import moments, numpy as np

from . import reading, predefined_models as predef, Optimize_Functions


def func_create_sfs(args):
    """
    Create an SFS from polymorphism data.
    """

    # Avoid having the file open in memory through 'argparse'.
    if args.intersect_sites:
        args.intersect_sites.close()
        # Assign the name of the input file.
        intersect_sites = args.intersect_sites.name
    else:
        intersect_sites = None

    # Record args into local variables. These are zipped lists.
    pop_labels = [str(a.split("=")[0]) for a in args.data_input]
    nchroms = [int(a.split("=")[1]) for a in args.data_input]
    fins = [str(a.split("=")[2]) for a in args.data_input]

    # Create the moments data dictionary.
    moments_dd = reading.build_moments_data_dict(
        fins=fins,
        pop_labels=pop_labels,
        nchroms=nchroms,
        intersect_sites=intersect_sites)

    sfs = moments.Spectrum.from_data_dict(
        moments_dd, polarized=False,
        pop_ids=pop_labels, projections=nchroms,
        # Remember to mask the corners before an expensive analysis?
        mask_corners=False)

    return {"sfs": sfs, "moments_dd": moments_dd}

def func_tosfs(args):
    """
    Exports a previously created SFS.
    """

    data = func_create_sfs(args)
    # Export the intermediate datadict for debugging reasons.
    with open("output.debug.datadict", "x") as fout:
        pprint.pp(data["moments_dd"], stream=fout)
    # Export the newly created SFS into a file. First makes sure the file does
    # not exist already. Can be reused within another `moments' instance.
    with open("output.sfs", "x") as fout:
        pass
    data["sfs"].to_file("output.sfs")

def func_optim(args):
    """
    Runs an optimization algorithm to fit SFS data.
    """

    # Check command line parameters are given correctly.
    if len(args.replicates) != int(args.rounds):
        raise Exception(
            "The list 'REPLICATES'"
            + " ({}) ".format(args.replicates)
            + " should be as long as the number of rounds "
            + "({}).".format(int(args.rounds)))
    if args.max_iters:
        if len(args.max_iters) != int(args.rounds):
            raise Exception(
                "The length of the list 'MAX_ITERS'"
                + " ({}) ".format(len(args.max_iters))
                + " should be the number of rounds"
                + " ({}).".format(int(args.rounds)))
    else:
        # Setup default max_iters if no value is given:
        args.max_iters = [50]
        while len(args.max_iters) < int(args.rounds):
            args.max_iters.append(np.ceil(
                # Geometric series starting from 50.
                args.max_iters[0] * ((1/2) ** len(args.max_iters))
            ))
        args.max_iters = args.max_iters[::-1]
    if args.fold_algor:
        if len(args.fold_algor) != int(args.rounds):
            raise Exception(
                "The length of the list 'FOLD' " +
                "({}) ".format(len(args.fold_algor)) +
                "should be the number of rounds " +
                "({}).".format(int(args.rounds)))
    else:
        # Setup default fold:
        args.fold_algor = list(range(int(args.rounds), 0, -1))
    # The arguments 'independent_runs' and 'manual_exec_i' are
    # mutually exclusive.
    if args.independent_runs:
        running = range(0, int(args.independent_runs))
    elif args.manual_exec_i:
        # In the case of 'manual_exec_i', a single execution will be computed.
        i = int(args.manual_exec_i)
        running = range(i, i + 1)

    # Read SFS file.
    sfs = func_create_sfs(args)["sfs"]
    # Make sure that the given SFS ends up being bidimensional;
    # first, at least a pair of pops are required.
    if len(sfs.pop_ids) < 2:
        raise Exception(
            "Population data (positional parameter) includes only a"
            + " single population. Make sure at least a pair of them"
            + " are given for the optimization of 2D-SFS.")
    # If more than 2 pop data, the SFS should be marginalized over 2.
    elif len(sfs.pop_ids) > 2:
        if not args.pop_lab:
            raise Exception(
                "Please provide a pair of pop. labels through"
                + " '--pop_lab' in order to isolate these and create a"
                + " bidimensional SFS from a multidimensional SFS.")
        if len(args.pop_lab) != 2:
            raise Exception(
                "Population data (positional parameter) includes more"
                + " than 2 populations. However, a pair of pop. labels"
                + " should be provided to marginalize these and create"
                + " a bidimensional SFS (see '--pop_lab' within '--help').")
        # Marginalize a pair of pops with the user-given pop labels.
        # `moments` cannot marginalize based on `pop_ids`; must translate
        # labels (str) to "list" indice (int) for slicing.
        margin_pop_indice = [sfs.pop_ids.index(label)
                             for label in sfs.pop_ids
                             if label not in args.pop_lab]
        # Marginalize (remove) population data we are not interested in.
        # Keep population data labelled as in "args.pop_lab".
        sfs = sfs.marginalize(margin_pop_indice)
        print("INFO: More than a pair of population data was included."
              + f" Keeping {args.pop_lab} while removing the rest.")

    # If an outfile prefix was not given through CLI, obtain it from a
    # combination of population IDs.
    if not args.out_prefix:
        args.out_prefix = "_".join(sorted(sfs.pop_ids))

    # Print summary (useful info. of data).
    print_summary(args, sfs)

    # Mask the corners; does it speed up computation?
    sfs.mask_corners()

    # Start the algor.
    for mod in args.fit_models:
        for run in running:

            Optimize_Functions.Optimize_Routine(
                fs=sfs,
                outfile=str(args.out_prefix) + f"_Exec{run}",
                model_name=str(mod),
                func=predef.models[mod]["func"],
                rounds=int(args.rounds),
                param_number=int(predef.models[mod]["param"]),
                fs_folded=sfs.folded,
                reps=args.replicates,
                maxiters=args.max_iters,
                folds=args.fold_algor,
                param_labels=predef.models[mod]["labels"],
                in_upper=predef.models[mod]["upper"],
                in_lower=predef.models[mod]["lower"],
            )

    return None

def print_summary(args, sfs):
    """
    Print useful info about the AFS or join-site FS.
    """

    print("\n# Overview of the 2D-SFS")
    print("========================")
    # Check if the SFS is folded or unfolded.
    if sfs.folded: kind_fold = "folded"
    else: kind_fold = "unfolded"
    print("· The 2D-SFS data has been interpreted"
          + " to be {}.".format(kind_fold))
    print("· Populations labels or IDs.: {}".format(sfs.pop_ids))
    print("· In the same order of pop. ids., their sample sizes"
          + " (nchr.): {}".format( sfs.sample_sizes))
    # Sum of 'sites' and 'segregating sites' rounded to two decimal places.
    print("· Sum of sites (including monomorph.) in the SFS: {}".format(
        np.around(sfs.sum(), 2)))
    print("· Sum of segregating sites: {}".format(
        np.around(sfs.S(), 2)))

    print("· The Fst between pops is '{}'.".format(sfs.Fst()))
    print("\n# Overview of algor. params.")
    print("============================")
    print("· User-specified steps of the optimization algorithm, " +
          "per round ('--max_iters'): " +
          str(args.max_iters))
    print("· User-specified perturbation of the starting parameters, " +
          "per round ('--fold_algor'): " +
          str(args.fold_algor))
    if args.independent_runs:
        print("· Performing {} independent runs.".format(
            args.independent_runs))
    elif args.manual_exec_i:
        print("· Performing a single run with ID/number " +
              "'{}'.".format(args.manual_exec_i))
    print("· Using the output file prefix '{}'".format(args.out_prefix))

    return None


if __name__ == '__main__':
    import argparse
    # Create the top-level parser.
    parser = argparse.ArgumentParser(
        "ephemeral",
        description=help_msg,
        # make sure the 'help_msg' is not automatically
        # wrapped at 80 characters (manually assign newlines).
        formatter_class=argparse.RawTextHelpFormatter)
    subparsers = parser.add_subparsers(
        required=True, help="Available subcommands.")

    # "TOSFS" SUBCOMMAND
    # ------------------
    parser_tosfs = subparsers.add_parser(
        "toSFS", # subcommand name.
        description="Create a multidimensional SFS from"
        + " polymorphism data files.",
        help="Create a multidimensional SFS from polymorphism"
        + " data files.", )
    parser_tosfs.set_defaults(func=func_tosfs)
    parser_tosfs.add_argument(
        "--intersect", type=argparse.FileType("r"), required=False,
        metavar="BED-LIKE", dest="intersect_sites", default=None,
        help="A BED-like file with sites which will intersect the"
        + " SNP files before creating an SFS.")
    parser_tosfs.add_argument(
        "data_input", type=str, nargs="+", metavar="LABEL=NCHR=FILE",
        help="Label of the population, number of chromosomes,"
        + " and file-name, separated by equal characters. Add as many"
        + " population data as required (will output multidim. SFS).")

#>    parser_tosfs.add_argument(
#>        "--polarized", action="store_const", const=True,
#>        default=False,
#>        help="Flag to create a polarized SFS input. If unset the"
#>        + " SFS will be unpolarized (folded)."
#>        + " DO NOT USE, NOT IMPLEMENTED WELL.")

    # "OPTIM" SUBCOMMAND
    # ------------------
    parser_optim = subparsers.add_parser(
        "optim", # subcommand name.
        description="Optimize population variables to fit an"
        + " observed bidimensional SFS.",
        help="Optimize population variables to fit an observed"
        + " bidimensional SFS.", )
    parser_optim.set_defaults(func=func_optim)
    # Files with polymorphism calling data from which an SFS will be created.
    parser_optim.add_argument(
        "data_input", type=str, nargs="+", metavar="LABEL=NCHR=FILE",
        help="Label of the population, number of chromosomes,"
        + " and file-name, separated by equal characters. Add as many"
        + " population data as required (will output multidim. SFS)."
        + " Make sure to provide at least 2 data files."
        + " Moreover, if more than 2 files are provided, it will"
        + " be required to marginalize and leave only a pair"
        + " (through '--pop_lab').")
    parser_optim.add_argument(
        "--fit_models", required=True, nargs="+",
        choices=list(predef.models.keys()), metavar="MODEL",
        help="The name of the model(s) desired to be fit."
        + " Run with '--fit_models h' to display available models."
        + " Lower and upper bounds for optimization can be tweaked "
        + " within the submodule 'predefined_models.py'.")
    parser_optim.add_argument(
        "--rounds", required=True, type=int,
        help="Amount of rounds of optimization.")
    parser_optim.add_argument(
        "--replicates", required=True, nargs="+", type=int,
        help="Replicates for each round; must be a list of ints "
        + "of length 'ROUNDS'.")

    execucions = parser_optim.add_mutually_exclusive_group(required=True)
    execucions.add_argument(
        "--independent_runs", type=int, metavar="INT",
        help="The number of independent runs of optimization performed.")
    execucions.add_argument(
        "--manual_exec_i", type=int, metavar="INT",
        help="If the independent runs are performed manually or on behalf"
        + " of other software, use this option for the"
        + " output files to be labeled accordingly.")

    # Optionally, marginalize data to obtain only a pair of pops.
    parser_optim.add_argument(
        "--pop_lab", required=False, nargs="+",
        metavar="ID", type=str, default=None,
        help="Population labels (matching those provided in the positional"
        + " argument). Only used to isolate a pair of pops. when there are"
        + " more than 2 polymorphism data files (ie. positional args);"
        + " otherwise not required. Make sure exactly 2 labels are provided.")
    # Prefix of the output file (followed by model and manual exec.).
    parser_optim.add_argument(
        "--out_prefix", required=False, metavar="FOUT", type=str, default=None,
        help="The prefix name of the output file(s). Naming convention:"
        + r" 'out_prefix' + 'Exec$N' + 'model' + '.log' (default: using"
        + " population labels joined by an undercase).")
    # Optionally, intersect sites provided within a BED file.
    parser_optim.add_argument(
        "--intersect", type=argparse.FileType("r"), required=False,
        metavar="BED-LIKE", dest="intersect_sites", default=None,
        help="A BED-like file with sites which will intersect the"
        + " SNP files before creating an SFS.")
    parser_optim.add_argument(
        "--max_iters", required=False, type=int, nargs="+", default=None,
        help="Steps of the optimization algorithm (default: the reverse"
        + " of a geometric series starting with '50').")
    parser_optim.add_argument(
        "--fold_algor", required=False, type=int, nargs="+", default=None,
        help="The perturbation of the starting parameters; must be a list"
        + " of ints of length 'ROUNDS' (default: from 'ROUNDS' to '1').")


    # To call a value, use `args.value` or `args.filename`.
    args = parser.parse_args()
    # Pass the args to the function assigned to each subcommand.
    args.func(args)

