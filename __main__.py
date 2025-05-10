#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# __main__.py
#
# 04 de gen. 2025  <adria@molevol-OptiPlex-9020>

help_msg = """
Create an SFS from pooled polymorphism data files and iteratively find some
population variables which fit these SFS.
"""

import os.path
import logging

import moments
import numpy as np
import matplotlib.pyplot as plt

from . import reading, predefined_models as predef, Optimize_Functions


def func_tosfs(args):
    """
    Create one or multiple SFS from polymorphism data. Exports it to an output
    directory.

    Args
    ----

    args.intersect_sites : str
    Input BED file with genomic coordinates.

    args.data_input : list
    A list with pairs of pop. label, num. of chromosomes and polymorphism data
    file.

    args.output_dir : str
    Inexistent directory name; create a new directory and write output there.
    """

    # Avoid having the file open in memory through 'argparse'.
    if args.intersect_sites:
        args.intersect_sites.close()
        # Assign the name of the input file.
        intersect_sites = args.intersect_sites.name
    else:
        intersect_sites = None
    # Check whether "args.output_dir" exists to avoid overwritting.
    if os.path.exists(args.output_dir):
        raise FileExistsError(f"The given output dir. '{args.output_dir}'"
                              + " already exists.")

    # Record args into local variables. These are zipped lists.
    pop_labels = [str(a.split("=")[0]) for a in args.data_input]
    nchroms = [int(a.split("=")[1]) for a in args.data_input]
    fins = [str(a.split("=")[2]) for a in args.data_input]

    # Create the requested single or multiple SFS.
    sfs_dict, logs = reading.main(
        finputs=fins, pop_labels=pop_labels, nchroms=nchroms,
        intersect_sites=intersect_sites)

    # Create a directory for writing output to.
    os.mkdir(args.output_dir)
    os.chdir(args.output_dir)
    # Store the logs.
    with open("input.log", "x") as fhandle:
        # Write a header.
        line = "\t".join(["filename"] + [str(k) for k in logs.keys()])
        line = "#" + str(line) + "\n"
        fhandle.write(line)

        # Write rows for each input data file (for each pop. lab.).
        for i in range(len(logs["poplabel"])):
            pl = logs["poplabel"][i]
            if pl in pop_labels:
                filename = fins[pop_labels.index(pl)]
            else:
                filename = None
            values = [str(val[i]) for val in logs.values()]
            line = "\t".join([str(filename)] + list(values))
            fhandle.write(line + "\n")

    # Store each SFS and information/logs (num. of sites, Fst, etc.).
    for pls, sfs in sfs_dict.items():
        # Write SFS to a file.
        outname = write_sfs(pls, sfs)
        # Plot SFS to a file.
        plot_sfs(pls, sfs)
        # Change the extension of the outname to ".log" instead of ".sfs" for
        # writing log information to.
        logname = outname.strip(".sfs") + ".log"
        # Create a summary.
        summary = sfs_summary(sfs)
        with open(logname, "x") as fhandle:
            fhandle.write(summary)

    # Plot all of the SFS in a single file (for comparing them).
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, layout="constrained")
    [ax.set_xlabel("Minor (folded) allele freq.")
     for ax in (ax1, ax2)]
    ax1.set_ylabel("Density")
    ax2.set_ylabel("Count of sites")
    for pls, sfs in sfs_dict.items():
        if len(pls) > 1: continue
        # Mask corners, i.e. remove sites with minor allele freq. = 0
        sfs.mask_corners()
        # Plot a density and total count graphs.
        ax1.plot(sfs/sfs.S(), label=str(pls[0]))
        ax2.plot(sfs, label=str(pls[0]))
    [ax.legend() for ax in (ax1, ax2)]
    [ax.set_xlim(left=1) for ax in (ax1, ax2)]
    [ax.set_xticks(list(ax.get_xticks())[1:] + [1]) for ax in (ax1, ax2)]
    plt.savefig("comparison.sfs.png")

    return None

def write_sfs(pop_labels, moments_sfs):
    if len(pop_labels) == 1:
        outname = "onedim." + ".".join(pop_labels) + ".sfs"
    elif len(pop_labels) == 2:
        outname = "bidim." + ".".join(pop_labels) + ".sfs"

    # Check file does not exist.
    with open(outname, "x"):
        pass
    # Write SFS to disk.
    moments_sfs.to_file(outname)

    return outname

def plot_sfs(pop_labels, moments_sfs):
    if len(pop_labels) == 1:
        moments_sfs.mask_corners()
        outname = "onedim." + ".".join(pop_labels) + ".sfs.png"
        fig, ax = plt.subplots(1, 1, figsize=(8, 4))
        ax.plot(moments_sfs, label=pop_labels)
        ax.legend()
        ax.set_xlim(left=1)
        ax.set_xticks(list(ax.get_xticks())[1:] + [1])
    elif len(pop_labels) == 2:
        outname = "bidim." + ".".join(pop_labels) + ".sfs.png"
        fig, ax = plt.subplots(1, 1, figsize=(8, 4))
        moments.Plotting.plot_single_2d_sfs(moments_sfs, ax=ax)

    plt.savefig(outname)
    plt.close()

    return outname

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

    # Import sfs file as a 'moments.Spectrum' object.
    sfs = moments.Spectrum.fromfile(args.sfs_input, mask_corners=False)

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
                + " '--pop-lab' in order to isolate these and create a"
                + " bidimensional SFS from a multidimensional SFS.")
        if len(args.pop_lab) != 2:
            raise Exception(
                "Population data (positional parameter) includes more"
                + " than 2 populations. However, a pair of pop. labels"
                + " should be provided to marginalize these and create"
                + " a bidimensional SFS (see '--pop-lab' within '--help').")

        print("INFO: More than a pair of population data was included.")
        sfs = marginalize_sfs(sfs, args.pop_lab)

    # If an outfile prefix was not given through CLI, obtain it from a
    # combination of population IDs.
    if not args.out_prefix:
        args.out_prefix = "_".join(sorted(sfs.pop_ids))

    # Print summary (useful info. of data).
    print(sfs_summary(sfs))

    # Print the options/config. selected for optim algorithm.
    print("\n# Overview of algor. params.")
    print("============================")
    print("· User-specified steps of the optimization algorithm, " +
          "per round ('--max-iters'): " +
          str(args.max_iters))
    print("· User-specified perturbation of the starting parameters, " +
          "per round ('--fold-algor'): " +
          str(args.fold_algor))
    if args.independent_runs:
        print("· Performing {} independent runs.".format(
            args.independent_runs))
    elif args.manual_exec_i:
        print("· Performing a single run with ID/number " +
              "'{}'.".format(args.manual_exec_i))
    print("· Using the output file prefix '{}'".format(args.out_prefix))

    # Mask the corners; does it speed up computation? Is it necessary?
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

def sfs_summary(sfs: moments.Spectrum):
    """
    For printing useful info about the AFS or join-site FS.
    """

    # Check if the SFS is folded or unfolded.
    if sfs.folded: kind_fold = "folded"
    else: kind_fold = "unfolded"

    summary = list([
        f"# Overview of the {len(sfs.pop_ids)}D-SFS",
        "========================",
        f"· The data has been interpreted to be {kind_fold}.",
        f"· Populations labels or IDs.: {sfs.pop_ids}",
        "· Sample sizes for these populations (~ nchr.), respectively: "
        + f"{sfs.sample_sizes}",
        "· Sum of sites (including monomorph.) in the SFS: "
        + f"{np.around(sfs.sum(), 2)}.",
        f"· Sum of segregating sites: {np.around(sfs.S())}.",
    ])
    # Summary stats only defined in 1D.
    if len(sfs.pop_ids) == 1:
        summary.append(f"· The Fst between pops: nan.")
        summary.append(f"· Diversity: {sfs.pi()}")
        summary.append(f"· Watterson's theta: {sfs.Watterson_theta()}")
    elif len(sfs.pop_ids) > 1:
        summary.append(f"· The Fst between pops: {sfs.Fst()}.")
        summary.append(f"· Diversity: nan.")
        summary.append(f"· Watterson's theta: nan.")
    # Add a citation.
    summary.append("")
    summary.append(
        "See"
        + " <https://momentsld.github.io/moments/sfs/sfs.html#computing-summary-statistics>"
    )

    # Join lines in the list with newline characters.
    summary = "\n".join(summary)

    return summary

def marginalize_sfs(sfs: moments.Spectrum, pop_lab: list):
    """
    Marginalize some pops. with the user-given pop labels.
    """

    # `moments` cannot marginalize based on `pop_ids`; must translate
    # labels (str) to "list" indice (int) for slicing.
    pop_indice_to_remove = [
        sfs.pop_ids.index(label) for label in sfs.pop_ids
        if label not in pop_lab]
    # Marginalize (remove) population data we are not interested in.
    # Keep population data with the labels found in "pop_lab".
    print(f"INFO: Keeping {pop_lab} while removing the rest.")
    sfs = sfs.marginalize(pop_indice_to_remove)

    return sfs

def func_margin(args):
    """
    Wrapper which uses 'marginalize_sfs' to marginalize and then writes to
    file.
    """

    # Import sfs file as a 'moments.Spectrum' object.
    sfs = moments.Spectrum.fromfile(args.sfs_input, mask_corners=False)

    sfs = marginalize_sfs(sfs, args.pop_lab)
    print("INFO: Writing to 'output.marginalized.sfs'.")
    with open("output.marginalized.sfs", "x") as fout:
        pass
    sfs.to_file("output.marginalized.sfs")

    return None


if __name__ == '__main__':
    import argparse
    # Create the top-level parser.
    parser = argparse.ArgumentParser(
        "ephemeral",
        description=help_msg,
        # Make sure the 'help_msg' is not automatically
        # wrapped at 80 characters. Newlines must be added manually.
        formatter_class=argparse.RawTextHelpFormatter)

    # The top-level parser does not read args., a call to a subcommand/subparser
    # is always required.
    subparsers = parser.add_subparsers(
        required=True, help="Available subcommands.")

    # "TOSFS" SUBPARSER
    # -----------------
    parser_tosfs = subparsers.add_parser(
        # Subparser name.
        "toSFS",
        description="Create 1D or/and 2D SFS from pooled polymorphism"
        + " data files. Accepts SYNC and NGSPOOL (from ngsJulia) formats.",
        help="Create 1D or/and 2D SFS from pooled polymorphism"
        + " data files. Accepts SYNC and NGSPOOL (from ngsJulia) formats.")
    parser_tosfs.set_defaults(func=func_tosfs)
    parser_tosfs.add_argument(
        "--intersect-sites", type=argparse.FileType("r"), required=False,
        metavar="BED-LIKE", dest="intersect_sites", default=None,
        help="A BED-like file with sites which will intersect the"
        + " SNP files before creating an SFS.")
    parser_tosfs.add_argument(
        "data_input", type=str, nargs="+", metavar="LABEL=NCHR=FILE",
        help="Label of the population, number of chromosomes,"
        + " and file-name, separated by equal characters. Add as many"
        + " population data as required (will output multidim. SFS).")
    parser_tosfs.add_argument(
        "--output-dir", required=True, metavar="DIR", type=str,
        help="New directory name in which the output will be"
        + " written to.")
    parser_tosfs.add_argument(
        "-v", "--verbose", action="store_true")

#>    parser_tosfs.add_argument(
#>        "--polarized", action="store_const", const=True,
#>        default=False,
#>        help="Flag to create a polarized SFS input. If unset the"
#>        + " SFS will be unpolarized (folded)."
#>        + " DO NOT USE, NOT IMPLEMENTED WELL.")

    # "MARGIN" SUBPARSER
    # ------------------
    parser_margin = subparsers.add_parser(
        "margin", # subcommand name.
        description="Marginalize (filter out) populations from a"
        + " multidimensional SFS (reducing dimensionality).",
        help="Marginalize (filter out) populations from a"
        + " multidimensional SFS (reducing dimensionality).", )
    parser_margin.set_defaults(func=func_margin)
    # SFS file which will be imported with 'moments' to be filtered.
    parser_margin.add_argument(
        "sfs_input", type=str, metavar="SFS",
        help="'Site Frequency Spectrum' obtained from the other"
        + " subcommand 'toSFS'.")
    parser_margin.add_argument(
        "--pop-lab", required=True, nargs="+", metavar="ID", type=str,
        help="Population labels (matching those provided within the SFS)"
        + " which will be isolated to a new SFS of reduced dimensionality.")

    # "OPTIM" SUBPARSER
    # -----------------
    parser_optim = subparsers.add_parser(
        "optim", # subcommand name.
        description="Optimize population variables to fit an"
        + " observed bidimensional SFS.",
        help="Optimize population variables to fit an observed"
        + " bidimensional SFS.", )
    parser_optim.set_defaults(func=func_optim)
    # SFS file which will be imported with 'moments' to be filtered.
    parser_optim.add_argument(
        "sfs_input", type=str, metavar="SFS",
        help="'Site Frequency Spectrum' obtained from the other"
        + " subcommand 'toSFS'. If it is above 2D (multidimensional),"
        + " then it is required to marginalize, leaving only 2 pops"
        + " (through specifying '--pop-lab').")
    parser_optim.add_argument(
        "--fit-models", required=True, nargs="+",
        choices=list(predef.models.keys()), metavar="MODEL",
        help="The name of the model(s) desired to be fit."
        + " Run with '--fit-models h' to display available models."
        + " Lower and upper bounds for optimization can be tweaked"
        + " within the submodule 'predefined_models.py'.")
    parser_optim.add_argument(
        "--rounds", required=True, type=int,
        help="Amount of rounds of optimization.")
    parser_optim.add_argument(
        "--replicates", required=True, nargs="+", type=int,
        help="Replicates for each round; must be a list of ints"
        + " of length 'ROUNDS'.")

    execucions = parser_optim.add_mutually_exclusive_group(required=True)
    execucions.add_argument(
        "--independent-runs", type=int, metavar="INT",
        help="The number of independent runs of optimization performed.")
    execucions.add_argument(
        "--manual-exec-i", type=int, metavar="INT",
        help="If the independent runs are performed manually or on behalf"
        + " of other software, use this option for the"
        + " output files to be labeled accordingly.")

    # Optionally, marginalize data to obtain only a pair of pops.
    parser_optim.add_argument(
        "--pop-lab", required=False, nargs="+",
        metavar="ID", type=str, default=None,
        help="Population labels (matching those provided in the positional"
        + " argument). Only used to isolate a pair of pops. when there are"
        + " more than 2 polymorphism data files (ie. positional args);"
        + " otherwise not required. Make sure exactly 2 labels are provided.")
    # Prefix of the output file (followed by model and manual exec.).
    parser_optim.add_argument(
        "--out-prefix", required=False, metavar="FOUT", type=str, default=None,
        help="The prefix name of the output file(s). Naming convention:"
        + r" 'out_prefix' + 'Exec$N' + 'model' + '.log' (default: using"
        + " population labels joined by an undercase).")
    parser_optim.add_argument(
        "--max-iters", required=False, type=int, nargs="+", default=None,
        help="Steps of the optimization algorithm (default: the reverse"
        + " of a geometric series starting with '50').")
    parser_optim.add_argument(
        "--fold-algor", required=False, type=int, nargs="+", default=None,
        help="The perturbation of the starting parameters; must be a list"
        + " of ints of length 'ROUNDS' (default: from 'ROUNDS' to '1').")



    # To call a value, use `args.value` or `args.filename`.
    args = parser.parse_args()

    # Configure logging. args.verbose might not exist.
    try:
        level = logging.DEBUG if args.verbose else logging.INFO
        logging.basicConfig(
            level=level,
            format='%(levelname)s.%(name)s.%(asctime)s.  %(message)s',
            datefmt='%H:%M:%S',
            force=True)
    except:
        logging.basicConfig(
            level=logging.INFO,
            format='%(levelname)s.%(name)s.%(asctime)s.  %(message)s',
            datefmt='%H:%M:%S',
            force=True)

    # Pass the args to the function assigned to each subcommand.
    args.func(args)

