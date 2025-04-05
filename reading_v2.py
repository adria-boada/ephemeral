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
# Create sqlite databases as temporary files.
import sqlite3
from tempfile import NamedTemporaryFile
from pathlib import Path
# Sanitize sqlite input instead of using placeholders.
import re
import os

import moments


# Decide the correct function to open an input file depending on extension.
opener = lambda finput: gzip.open(finput, "rt") if finput.endswith(".gz") \
    else open(finput, "rt")
# Strip and split lines read from a file.
ss = lambda line, sep=None: line.strip("\n").split(sep)


def generator_loci(finput: str, nchr: int, lrt_snp_thresh: float=6.64):
    """
    Read a file with POOLED polymorphism data. Can be formatted as an
    "*ngsPool.out" or a "*sync" (possibly gzipped). Yields dictionaries with
    data per row.

    Usage: 'for row_data in generator_loci(file, nchr): process row_data'

    Input
    -----

    finput : str
    The path/filename to an "*.ngsPool.out" or a "*.sync" file, either
    uncompressed or gzipped.

    nchr : int
    The number of sets of chromosomes in the pool (depth). It is twice the
    amount of sequenced individuals for diploid organisms.

    lrt_snp_thresh : float
    Threshold to call SNPs if "finput" is formatted as an "*.ngsPool.out". See
    `ngsJulia` documentation. Otherwise not used for "*.sync".
    """

    # Function which creates a list of [0/nchr, 1/nchr, etc. nchr/nchr], which
    # is a discrete distribution with "nchr + 1" bins of possible SFS
    # frequencies. It will help in translating allele freq. to allele counts.
    to_sfs_indices = lambda nchr: [i / nchr for i in range(0, nchr + 1)]
    # Initialize it.
    sfs_indice = to_sfs_indices(nchr)
    # Compute the difference between site allele frequency "saf" and each indice
    # within the list "sfs_indice".
    to_abs_diff = lambda saf, sfs_indice: [abs(saf - i) for i in sfs_indice]
    # Finally, return the allele count: the indice of the minimum of "abs_diff".
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
            for numline, line in enumerate(fhandle):
                # Split the line by blank space into a list of fields.
                line = ss(line)
                # Parse the site "(chr, pos)".
                fields = {"chr": str(line[0]), "pos": int(line[1]), }
                # First, check that the reference is unambigous. Otherwise the
                # "sync string" will not be parseable, i.e. ".:.:.:.".
                if str(line[2]) not in ("A", "T", "C", "G"):
                    print("WARNING: Missing data or ambigous reference at "
                          + str((fields["chr"], fields["pos"]))
                          + " in file '" + str(finput) + "'; skipping site.",
                          file=sys.stderr)
                    continue
                # Parse the polymorphism data.
                fields["allele_depths"] = from_sync_to_allele_depths(line[3])
                # If the site is triallelic (more than 2 allele counts), or the
                # site has missing data (no allele counts), skip the site.
                if len(fields["allele_depths"]) > 2:
                    print("WARNING: Triallelic site at "
                          + str((fields["chr"], fields["pos"]))
                          + " in file '" + str(finput) + "'; skipping site.",
                          file=sys.stderr)
                    continue
                elif len(fields["allele_depths"]) == 0:
                    print("WARNING: Missing data at "
                          + str((fields["chr"], fields["pos"]))
                          + " in file '" + str(finput) + "'; skipping site.",
                          file=sys.stderr)
                    continue
                # If there is one allele, yield it as a dictionary.
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
                # If there are two alleles (polymorphic sites), yield them as a
                # dictionary.
                elif len(fields["allele_depths"]) == 2:
                    (min_nuc, min_depth), (maj_nuc, maj_depth) = \
                        fields["allele_depths"]
                    min_freq = min_depth / (maj_depth + min_depth)
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

        # Finally, print how many lines have been read.
        print(f"INFO: Finished! Read {numline + 1} sites/lines"
              + f" from '{finput}' (inluding triallelic or"
              + " missing data sites).")

    elif "out" in finput or "ngsPool" in finput:
        print("INFO: Reading the file '" + str(finput)
              + "' as an ngsPool file from ngsJulia.")
        # Create a function which categorises the fields from a given line.
        to_fields = lambda ssline: {
            "chr": str(ssline[0]), "pos": int(ssline[1]),
            "maj_allele": str(ssline[4]), "min_allele": str(ssline[5]),
            "min_saf_mle": float(ssline[10]), "lrt_snp": float(ssline[6]), }

        def ngsjulia_to_datadict(line):
            line = ss(line)
            fields = to_fields(line)
            # If Likelihood-Ratio-Test (LRT) for SNPs is higher than the
            # threshold, accept it as a polymorphism (see ngsJulia docs).
            if fields["lrt_snp"] > lrt_snp_thresh:
                min_count = to_allele_count(
                    to_abs_diff(fields["min_saf_mle"], sfs_indice))
                min_freq = fields["min_saf_mle"]
                min_nuc = fields["min_allele"]
            # If LRT is lower than the threshold the site is monomorphic.
            else:
                min_count = 0
                min_freq = 0
                min_nuc = "N"
            return {
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

        with opener(finput) as fhandle:
            # Read the first line, which may be a header.
            line = fhandle.readline()
            # If it is not a header, parse it.
            add_to_numlines = 0
            if not "##chr" in line:
                add_to_numlines = 1
                yield ngsjulia_to_datadict(line)
            # Start reading and parsing the rest of the lines in the file.
            for numline, line in enumerate(fhandle):
                yield ngsjulia_to_datadict(line)
            numline += add_to_numlines

        # Finally, print how many lines have been read.
        print(f"INFO: Finished! Read {numline + 1} sites/lines"
              + f" from '{finput}'.")

def generator_sites(finput: str, chr_column: int=0, pos_column: int=1):
    """
    Read a BED file. Yield pairs of "chromosome" and "site" per row.

    Input
    -----

    finput : str
    File input name. Comments and headers within the BED file must start
    with "#".

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
            yield {
                "chr": chromosome,
                "pos": pos, }

def sanitized_sqlite(parameter: str):
    if not re.fullmatch( r"[a-zA-Z0-9_]+", parameter):
        raise Exception(f"The parameter '{parameter}' is not"
                        + " sanitized for it to be used with"
                        + " sqlite3; it must be a string"
                        + " matching '[a-zA-Z0-9_]+'.")

def populate_new_sqlite_table(generator: iter, tblname: str, db_sqlite: str):
    """
    Input
    -----

    finput : str
    File input name. Accepts BED, *ngsPool.out and *sync. Can be uncompressed or
    gzipped.

    tblname : str
    Name of the newly created table within the sqlite database file.

    db_sqlite : str
    File name of the sqlite database which this function will write to.
    """

    # Open sqlite database.
    con = sqlite3.connect(db_sqlite)
    cur = con.cursor()
    # Look at the fields/columns returned by "generator".
    first_row = next(generator)
    # Make sure they include required "chromosome" and "pos" fields.
    if "chr" not in first_row.keys() or "pos" not in first_row.keys():
        raise Exception("'Generator' from input file does not have"
                        + " 'chr' or 'pos' fields.")
    # Try to avoid SQL injections...
    [sanitized_sqlite(param) for param in (tblname, *first_row.keys())]
    # Create a table with the name "tblname" and the fields/columns in
    # "first_row.keys()".
    cur.execute(f"CREATE TABLE {tblname} ("
                + ", ".join(first_row.keys()) + ", "
                + " PRIMARY KEY(chr, pos) )")
    # Add the data in "first_row" to the newly created table.
    command_insert_sqlite = str(
        f"INSERT INTO {tblname} VALUES ("
        + ", ".join([":" + key for key in first_row.keys()]) + ")")
    cur.execute(command_insert_sqlite, first_row)
    # Add the data in the rest of the rows of the file input.
    for row in generator:
        cur.execute(command_insert_sqlite, row)
    # Commit these INSERT changes.
    con.commit()

def _generator_moments_dd(finputs: list, pop_labels: list, nchroms: list,
                          intersect_sites: str=None):
    """
    Yield moments data dictionary, allowing little footprint on memory usage
    when creating a moments.Spectrum.
    """

    # Check that the three lists in parameters are of the same length.
    if len(finputs) != len(pop_labels) or len(finputs) != len(nchroms):
        raise Exception("The lists of file inputs, pop. labels and num."
                        + " of chromosomes must be of the same length.")
    # Create a temporary file where an sqlite3 database will be stored.
    # Avoid using the directory/partition "/tmp", because the database file
    # could become too big to fit within there. Instead, place this file in the
    # home directory.
    home = Path.home()
    temp_file = NamedTemporaryFile(dir=home, prefix="ephemeral.tmp")
    print("INFO: Storing an sqlite3 database with the polymorphism data"
          + f" in a temporary file at '{temp_file.name}'.")
    # Create a table within the database for each input file.
    for fi, nc, pl in zip(finputs, nchroms, pop_labels):
        generator = generator_loci(fi, nc)
        populate_new_sqlite_table(generator, pl, temp_file.name)
    # Add a table with sites from "intersect sites".
    if intersect_sites:
        generator = generator_sites(intersect_sites)
        # Remove directory and extensions from filename so it can be used as
        # table name.
        intersect_tblname = os.path.split(intersect_sites)[1].split(".")[0]
        populate_new_sqlite_table(generator, intersect_tblname,
                                  temp_file.name)

    # Select an inner join of all the tables, equivalent to sites shared across
    # all files/tables.
    con = sqlite3.connect(temp_file.name)
    cur = con.cursor()
    command_intersection_sqlite = str(
        "SELECT ROW_NUMBER() OVER(ORDER BY chr, pos), *"
        + f" FROM {pop_labels[0]} "
        + " ".join(["INNER JOIN " + pl + " USING (chr, pos)"
                    for pl in pop_labels[1:]]))
    if intersect_sites:
        command_intersection_sqlite += str(
            f" INNER JOIN {intersect_tblname} USING (chr, pos)")
    print("INFO: Executing the inner join (finding sites shared across"
          + " all of the given files).")
    #> print(f"DEBUG: {command_intersection_sqlite}")

    # Execute the command to the open sqlite3 database.
    # The selected rows with polymorphism data are stored in "res".
    res = cur.execute(command_intersection_sqlite)

    # Initialize a dictionary which will be later read by moments and converted
    # to an SFS.
    moments_dd = dict()
    print("INFO: Building a 'moments data dictionary'.")
    # Iterate across each row (sites) of the inner joined tables.
    for row in res:
        # (chr, pos) in second and third fields.
        site = tuple(row[1:3])
        # Major allele in the fourth field, and repeats every eight.
        maj_alleles = tuple([
            str(row[i]) for i in range(3, 8 * len(pop_labels) + 3, 8)])
        # Minor allele in the eighth field, and repeats every eight.
        min_alleles = tuple([
            str(row[i]) for i in range(7, 8 * len(pop_labels) + 7, 8)])
        # Obtain all of the observed alleles at this site.
        alleles = list()
        for a in maj_alleles + min_alleles:
            if a not in alleles and a != "N": alleles.append(a)
        alleles = tuple(alleles)
        # If the site is triallelic, skip it.
        if len(alleles) > 2: continue
        # If the site is monomorphic, add an ambigous nucleotide for formatting
        # reasons.
        if len(alleles) == 1:
            alleles = (alleles[0], "N")
        # Major allele count in the seventh field, and repeats every eight.
        maj_counts = tuple([
            int(row[i]) for i in range(6, 8 * len(pop_labels) + 6, 8)])
        # Minor allele count in the eleventh field, and repeats every eight.
        min_counts = tuple([
            int(row[i]) for i in range(10, 8 * len(pop_labels) + 10, 8)])
        # Compute calls for each population.
        calls = {
            pop_lab:
            (maj_count, min_count) if maj_allele == alleles[0]
            else (min_count, maj_count)
            for pop_lab, maj_allele, maj_count, min_count in zip(
                pop_labels, maj_alleles, maj_counts, min_counts)
        }
        # Store the information for this row in "moments_dd".
        moments_dd[site] = {
            "context": "---",
            "outgroup_context": "---",
            "outgroup_allele": "-",
            "segregating": alleles,
            "calls": calls, }
        # Yield moments_dd every 10 000 rows to avoid overflowing memory with a
        # massive dictionary.
        if row[0] % 10000 == 0:
            print(f"INFO: Reached {row[0]} sites.")
            yield moments_dd
            moments_dd = dict()

    # Yield the last moments_dd.
    print(f"INFO: Reached the end at {row[0]} sites.")
    if moments_dd: yield moments_dd

def build_moments_sfs(finputs: list, pop_labels: list, nchroms: list,
                      intersect_sites: str=None):
    """
    Build a 'moments.Spectrum' from calling data.
    """

    # Associate pop. labels with their num. of chr.
    poplab_to_nchrom = dict()
    for pl, nc in zip(pop_labels, nchroms):
        poplab_to_nchrom[str(pl)] = int(nc)

    # Compute combinations of populations to create 1D and 2D (joint) SFS for
    # each one of them.
    def combinations(arr):
        if len(arr) == 0: return [[]]
        combs = []
        for c in combinations(arr[1:]):
            combs += [c, c + [arr[0]]]
        return combs
    # Only combinations of 1 or 2 pops. allowed.
    combinations_pop = [c for c in combinations(pop_labels)
                            if len(c)==2 or len(c)==1]
    print("INFO: Creating {} 1D-SFS and {} 2D-SFS.".format(
        len([c for c in combinations_pop if len(c) == 1]),
        len([c for c in combinations_pop if len(c) == 2])))

    moments_sfs = dict()
    for poplabs in combinations_pop:
        moments_sfs[tuple(poplabs)] = list()

    # Retrieve polymorphism data.
    gen_moments_dd = _generator_moments_dd(
        finputs, pop_labels, nchroms, intersect_sites)

    # Create all possible combinations of 1D or 2D SFS with this polymorphism
    # data.
    for moments_dd in gen_moments_dd:
        for poplabs in moments_sfs.keys():
            # Filter "moments_dd" by removing calls not matching "poplabs".
            filtered_moments_dd = dict()
            for site, data in moments_dd.items():
                filtered_calls = {pl: calls
                                  for pl, calls in data["calls"].items()
                                  if pl in poplabs}
                filtered_moments_dd[site] = data.copy()
                filtered_moments_dd[site]["calls"] = filtered_calls

            # From the filtered "moments_dd" create an SFS and append it to the
            # "moments_sfs" dictionary with the key "poplabs".
            moments_sfs[poplabs].append(
                moments.Spectrum.from_data_dict(
                    filtered_moments_dd, poplabs,
                    [poplab_to_nchrom[pl] for pl in poplabs],
                    mask_corners=False, polarized=False))
            # Summing the sites of a pair of SFS into a single SFS.
            # It reduces the amount of data in memory.
            moments_sfs[poplabs] = [sum(moments_sfs[poplabs])]

    # PRINT STATS OF THE NEWLY CREATED SFS?

    # Return all of the SFS. Remember to remove the temporary list the SFS are
    # in (i.e. slice and return the first item).
    return {key: val[0] for key, val in moments_sfs.items()}



    # RETALLS!!!
    ...

    # Create an output SFS filename from a list/tuple of pop. labels.
    sfs_filename = lambda poplabs: (
        str(poplabs[0]) + ".sfs" if len(poplabs) == 2
        else
        "joint." + str(poplabs[0]) + "." + str(poplabs[1]) + ".sfs")


    # Create a list of moments.Spectrum for pops in "moments_sfs.keys()".

    moments_dd = next(gen_moments_dd)
    # Initialize SFS for each population and for each pair of pops.
    for comb in combinations_sfs:
        ...
        # CREATE FILENAME:
    moments_sfs = moments.Spectrum.from_data_dict(
        moments_dd, pop_labels, nchroms,
        mask_corners=False, polarized=False)

    for moments_dd in gen_moments_dd:
        moments_sfs += moments.Spectrum.from_data_dict(
            moments_dd, pop_labels, nchroms,
            mask_corners=False, polarized=False)

    # Create a folder to store definitive SFS.
    try:
        os.makedirs(outfolder)
    except OSError:
        raise OSError(f"The path to create '{outfolder}' already exists.")

    return moments_sfs











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


