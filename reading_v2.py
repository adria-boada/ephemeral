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
opener = lambda fin: gzip.open(fin, "rt") if fin.endswith(".gz") \
    else open(fin, "rt")


def generator_sites(finput: str, chr_column: int=0, pos_column: int=1):
    """
    Create an iterable generator composed of pairs of "chromosome" and "site"
    from a BED-like file.

    Input
    -----

    finput : str
    File input name.

    chr_column : int
    The index (zero-based) of the column within "finput" with "chromosome" info.

    pos_column : int
    The index (zero-based) of the column within "finput" with "site" info.
    """

    # Create a local function that splits the rows of the input file.
    ss = lambda line: line.strip("\n").split()
    
    with opener(finput) as fhandle:
        for line in fhandle:
            line = ss(line)
            # Slice the line by provided column indices.
            chrom, pos = str(line[chr_column]), int(line[pos_column])

            yield chrom, pos

def intersect_sites(sites: dict, finp: str,
                    chrom_column: int=0, pos_column: int=1):
    """
    Find the intersection between the sites in 'sites' and the sites in 'finp'.

    Input
    -----

    sites : dict
    A dictionary returned from 'read_bedlike_sites'.

    finp : str
    Name of an input file.

    chrom_column : int
    The index (zero-based) of the column of 'finp' with 'chromosome'.

    pos_column : int
    The index (zero-based) of the column of 'finp' with 'position'.
    """

    # Create a local func that splits the rows of the input file.
    ss = lambda line: line.strip("\n").split()
    # Initialize an empty dict.
    new_sites = dict()

    with opener(finp) as fhandle:
        # If 'finp' is ngsPool, skip its header.
        if finp.endswith("ngsPool.out") or finp.endswith("ngsPool.out.gz"):
            fhandle.readline()
        for line in fhandle:
            # Strip newline and split by blank spaces.
            line = ss(line)
            # Slice line by column indices.
            chrom, pos = str(line[chrom_column]), int(line[pos_column])
            # If the chromosome or position is not in 'sites', skip it
            # (this site does not belong to the intersection).
            if chrom not in sites.keys():
                continue
            elif pos not in sites[chrom]:
                continue
            # Otherwise, the site belongs to the intersection. Add it.
            if chrom not in new_sites.keys():
                new_sites[chrom] = set()
            new_sites[chrom].add(pos)

    return new_sites

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


