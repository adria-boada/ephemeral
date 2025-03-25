#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# reading.v2.py
#
# 19 de març 2025  <adria@molevol-OptiPlex-9020>

help_msg = """
Module which reads an input file with called polymorphisms from sequencing data
and uses it to create input compatible with `moments`.

Available formats for conversion: *sync or *ngsPool.out (output from ngsJulia).
"""

import sys
# Reading gzipped files.
import gzip
# Assess file sizes.
import os


# Decide the correct function to open an input file depending on extension.
opener = lambda finput: gzip.open(finput, "rt") if finput.endswith(".gz") \
    else open(finput, "rt")
# Strip and split lines read from a file.
ss = lambda line, sep=None: line.strip("\n").split(sep)


def generator_loci(finput: str, nchr: int, lrt_snp_thresh: float=6.64):
    """
    Read a file with polymorphism data. Can be formatted as an "*ngsPool.out" or
    a "*sync" (possibly gzipped). Yields dictionaries with data per row.

    Usage: 'for data in generator_loci(file, nchr): analyze data'

    Input
    -----

    finput : str
    The path/filename to an "*.ngsPool.out" or a "*.sync" file, either
    uncompressed or gzipped.

    nchr : int
    The number of sets of sequenced chromosomes (depth). It is twice the amount
    of sequenced individuals for diploid organisms.

    lrt_snp_thresh : float
    Threshold to call SNPs if "finput" is formatted as an "*.ngsPool.out". See
    `ngsJulia` documentation. Otherwise not used for "*.sync".
    """

    # Function which creates a list of [0/nchr, 1/nchr, etc. nchr/nchr].
    # It will help in translating allele frequency to allele counts.
    to_sfs_indices = lambda nchr: [i / nchr for i in range(0, nchr + 1)]
    # Initialize a list of frequencies within the SFS with the previous
    # function.
    sfs_indice = to_sfs_indices(nchr)
    # Compute the difference between site allele frequency `saf` and each indice
    # within the list `sfs_indice`.
    to_abs_diff = lambda saf, sfs_indice: [abs(saf - i) for i in sfs_indice]
    # Finally, return the allele count (minimum of `abs_diff`).
    to_allele_count = lambda abs_diff_list: abs_diff_list.index(
        min(abs_diff_list))

    if "sync" in finput:
        print("INFO: Reading the file '" + str(finput)
              + "' as a SYNC formatted file.")

        # Initialize a function to parse the "sync string" into a list of pairs
        # of alleles (nucleotide) and allele counts.
        from_sync_to_allele_depths = lambda sync_string: sorted((
            (nuc, int(count)) for count, nuc in zip(
                sync_string.split(":")[0:4], ("A", "T", "C", "G"))
            if int(count) > 0), key=lambda t: t[1])

        with opener(finput) as fhandle:
            for line in fhandle:
                # Split the line by blank space into a list of fields.
                line = ss(line)
                # First, check that the reference is unambigous. Otherwise the
                # "sync string" will not be parseable, i.e. ".:.:.:.".
                if str(line[2]) not in ("A", "T", "C", "G"):
                    print("WARNING: Missing data or ambigous reference at "
                          + str((fields["chr"], fields["pos"]))
                          + " in file '" + str(finput) + "'; skipping it.")
                    continue
                # Parse the information of interest.
                fields = {
                    "chr": str(line[0]), "pos": int(line[1]),
                    "allele_depths": from_sync_to_allele_depths(line[3]), }
                # If the site is triallelic (more than 2 allele counts), or the site
                # has missing data (no allele counts), skip the site.
                if len(fields["allele_depths"]) > 2:
                    print("WARNING: Triallelic site at "
                          + str((fields["chr"], fields["pos"]))
                          + " in file '" + str(finput) + "'; skipping it.")
                    continue
                elif len(fields["allele_depths"]) == 0:
                    print("WARNING: Missing data at "
                          + str((fields["chr"], fields["pos"]))
                          + " in file '" + str(finput) + "'; skipping it.")
                    continue
                # Yield only one allele (monomorphism).
                elif len(fields["allele_depths"]) == 1:
                    nuc, depth = fields["allele_depths"][0]
                    yield {
                        "chr": fields["chr"],
                        "pos": fields["pos"],
                        "maj_allele": nuc,
                        "maj_depth": depth,
                        "maj_freq": 1,
                        "maj_count": nchr,
                        "min_allele": "N",
                        "min_depth": 0,
                        "min_freq": 0,
                        "min_count": 0, }
                # Yield two alleles (polymorphic site).
                elif len(fields["allele_depths"]) == 2:
                    (min_nuc, min_depth), (maj_nuc, maj_depth) = \
                        fields["allele_depths"]
                    min_freq = min_depth / sum([maj_depth, min_depth])
                    min_count = to_allele_count(
                        to_abs_diff(min_freq, sfs_indice))
                    # For each row in the file "finput", yield a dictionary with
                    # the data in this row.
                    yield {
                        "chr": fields["chr"],
                        "pos": fields["pos"],
                        "maj_allele": maj_nuc,
                        "maj_depth": maj_depth,
                        "maj_freq": float(1 - min_freq),
                        "maj_count": int(nchr - min_count),
                        "min_allele": min_nuc,
                        "min_depth": min_depth,
                        "min_freq": float(min_freq),
                        "min_count": int(min_count), }

    elif "out" in finput or "ngsPool" in finput:
        print("INFO: Reading the file '" + str(finput)
              + "' as an ngsPool file from ngsJulia.")
        # Create a function which categorises the fields from a given line.
        to_fields = lambda ssline: {
            "chr": str(ssline[0]), "pos": int(ssline[1]),
            "maj_allele": str(ssline[4]), "min_allele": str(ssline[5]),
            "min_saf_mle": float(ssline[10]), "lrt_snp": float(ssline[6]), }

        with opener(finput) as fhandle:
            # Read the first line, which may be a header.
            line = fhandle.readline()
            # If it is not a header, parse it.
            if not "##chr" in line:
                line = ss(line)
                fields = to_fields(line)
                # If Likelihood-Ratio-Test (LRT) for SNPs is higher than the
                # threshold, accept it as a polymorphism.
                if fields["lrt_snp"] > lrt_snp_thresh:
                    min_count = to_allele_count(
                        to_abs_diff(fields["min_saf_mle"], sfs_indice))
                    min_freq = fields["min_saf_mle"]
                    min_nuc = fields["min_allele"]
                else:
                    min_count = 0
                    min_freq = 0
                    min_nuc = "N"
                yield {
                    "chr": fields["chr"],
                    "pos": fields["pos"],
                    "maj_allele": fields["maj_allele"],
                    "maj_depth": None,
                    "maj_freq": float(1 - min_freq),
                    "maj_count": int(nchr - min_count),
                    "min_allele": min_nuc,
                    "min_depth": None,
                    "min_freq": min_freq,
                    "min_count": min_count, }
            # Read and parse the rest of the lines in the file.
            for line in fhandle:
                line = ss(line)
                fields = to_fields(line)
                # If Likelihood-Ratio-Test (LRT) for SNPs is higher than the
                # threshold, accept it as a polymorphism.
                if fields["lrt_snp"] > lrt_snp_thresh:
                    min_count = to_allele_count(
                        to_abs_diff(fields["min_saf_mle"], sfs_indice))
                    min_freq = fields["min_saf_mle"]
                    min_nuc = fields["min_allele"]
                else:
                    min_count = 0
                    min_freq = 0
                    min_nuc = "N"
                yield {
                    "chr": fields["chr"],
                    "pos": fields["pos"],
                    "maj_allele": fields["maj_allele"],
                    "maj_depth": None,
                    "maj_freq": float(1 - min_freq),
                    "maj_count": int(nchr - min_count),
                    "min_allele": min_nuc,
                    "min_depth": None,
                    "min_freq": min_freq,
                    "min_count": min_count, }

def generator_sites(finput: str, chr_column: int=0, pos_column: int=1):
    """
    Read a BED-like file. Yield pairs of "chromosome" and "site" per row.

    Input
    -----

    finput : str
    File input name.

    chr_column : int
    The index (zero-based) of the column within "finput" with "chromosome" info.

    pos_column : int
    The index (zero-based) of the column within "finput" with "site" info.
    """

    print(f"INFO: Reading sites ('chr' and 'pos' pairs) from '{finput}'.")
    with opener(finput) as fhandle:
        for line in fhandle:
            line = ss(line)
            # Make sure this line is not a header or a comment.
            if "#" in line[0]: continue
            # Slice the line by provided column indices.
            try: chromosome, pos = str(line[chr_column]), int(line[pos_column])
            except ValueError as e:
                raise e(f"'{pos}' was retrieved as a chromosome position"
                        + f" from '{finput}'. It could not be tipified as an"
                        + " 'int'; check if the correct column indice is"
                        + f" '{pos_column}' and also if there are comments"
                        + " and headers which do not start with '#' as"
                        + " their first character.")
            yield chromosome, pos










def intersect_sites(sites: dict, finput: str,
                    chr_column: int=0, pos_column: int=1):
    """
    Find the intersection between the sites in 'sites' and the sites in 'finp'.

    Input
    -----

    sites : dict
    A dictionary returned from 'read_bedlike_sites'.

    finp : str
    Name of an input file.

    chr_column : int
    The index (zero-based) of the column of 'finp' with 'chromosome'.

    pos_column : int
    The index (zero-based) of the column of 'finp' with 'position'.
    """

    new_sites = dict()
    for chromosome, pos in generator_sites(finput, chr_column, pos_column):
        # If the chromosome or position is not in 'sites', skip it
        # (this site does not belong to the intersection).
        if chromosome not in sites.keys(): continue
        elif pos not in sites[chromosome]: continue
        # Otherwise, the site belongs to the intersection. Add it.
        if chromosome not in new_sites.keys():
            new_sites[chromosome] = set([pos])
        else:
            new_sites[chromosome].add(pos)

    return new_sites

def generator_intersecting(first_file: str, other_files: list,
                           chunksize: int=int(1e6)):
    """
    Create an iterable generator composed of groups of sites shared between all
    given files.
    """

    sites = dict()
    memory_warning = True

    for numline, (chromosome, pos) in enumerate(generator_sites(first_file)):
        if chromosome not in sites.keys():
            sites[chromosome] = set([pos])
        else:
            sites[chromosome].add(pos)
        # Avoid loading too many sites to RAM. Analyze in chunks.
        if numline >= chunksize:
            if memory_warning:
                print(f"INFO: Potentially more than '{chunksize}' shared sites"
                      + " between files; will analyze them chunk by chunk"
                      + " to avoid taking up too much memory.")
                memory_warning = False
            for finput in other_files:
                print("INFO: Finding shared sites between in-memory and"
                      + f" '{finput}'.")
                sites = intersect_sites(sites, finput)

            yield sites

    # Yield remaining sites.
    for finput in other_files:
        print("INFO: Finding shared sites between in-memory and"
              + f" '{finput}'.")
        sites = intersect_sites(sites, finput)

    yield sites








def find_intersection_sites(finps: list, intersect_bedlike: str=None):
    """
    With the sites belonging to 'finps' and 'intersect_bedlike', find their
    intersection; sites appearing in all files, shared between all files.

    Input
    -----

    finps : list
    List with at least one input file name.

    intersect_bedlike : str
    A 'BED-like' file which will be added to the requirements to intersect.
    """

    # If 'intersect_bedlike' parameter was not provided, some starting sites
    # will have to be initialized.
    if not intersect_bedlike:
        # Find the name of the file with the smallest size. It will be the
        # starting 'sites'.
        file_sizes = tuple(map(os.path.getsize, finps))
        index_smallest_file = file_sizes.index(min(file_sizes))
        smallest_file_name = finps[index_smallest_file]
        print(f"INFO: Reading the sites in {smallest_file_name}.")
        # Both ngsPool and sync files have coordinates in first and second
        # columns.
        sites = read_bedlike_sites(smallest_file_name)
        # Remove these sites from the input files (already initialized).
        finps = [filename for filename in finps
                 if filename != smallest_file_name]
    else:
        print(f"INFO: Reading the sites in {intersect_bedlike}.")
        sites = read_bedlike_sites(intersect_bedlike)

    # Compare the sites in the first file with the second. Then with the third,
    # etc. Remove sites that are not shared.
    for filename in finps:
        print("INFO: Finding intersection of sites with previous sites"
              + f" and the sites in the file {filename}.")
        # Replace 'sites' with the intersection between itself and 'filename'.
        sites = intersect_sites(sites, filename)

    return sites

def initialize_moments_dd(sites: dict):
    """
    Initialize a moments 'data_dictionary' with blank fields/keys: 'context',
    'outgroup_context' and 'outgroup_allele' are empty (i.e. equal to '---').
    """


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Si es crida com a script:
if __name__ == '__main__':

    ### PARSE ARGUMENTS ###
    import argparse
    parser = argparse.ArgumentParser(description=help_msg,
        # make sure the 'help_msg' is not automatically
        # wrapped at 80 characters (manually assign newlines).
        formatter_class=argparse.RawTextHelpFormatter)
    # file-name: positional arg.
#    parser.add_argument('filename', type=str, help='Path to ... file-name')
    # integer argument
#    parser.add_argument('-a', '--numero_a', type=int, help='Paràmetre "a"')
    # choices argument
#    parser.add_argument('-o', '--operacio', 
#            type=str, choices=['suma', 'resta', 'multiplicacio'],
#            default='suma', required=False,
#            help='')

    args = parser.parse_args()
    # call a value: args.operacio or args.filename.


