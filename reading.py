#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Reading.py
#
# 20 de des. 2024  <adria@molevol-OptiPlex-9020>

help_msg = """
Module which reads an input file with called polymorphisms from sequencing data
and uses it to create a Python3 dictionary, which is compatible with `moments`.

Available formats for conversion: *sync or *ngsPool.out (output from ngsJulia).
"""

import sys
# Reading gzipped files.
import gzip
# Subsample data and create pseudoreplicates through bootstrapping.
import random


# Function which creates a list of [0/nchr, 1/nchr, etc. nchr/nchr]
to_sfs_indices = lambda nchr: [i / nchr for i in range(0, nchr + 1)]
# Compute the difference between site allele frequency `saf` and each indice
# within the list `sfs_idx`.
to_abs_diff = lambda saf, sfs_idx: [abs(saf - i) for i in sfs_idx]
# Compute the allele count (minimum of `abs_diff`).
to_allele_count = lambda abs_diff: abs_diff.index(min(abs_diff))
# Decide the correct function to open an input file depending on extension.
opener = lambda fin: gzip.open(fin, "rt") if fin.endswith(".gz") else open(fin)


def read_bedlike_sites(fin: str):
    """
    Read a BED-like file and return a `set` of sites with ("Chr", "Pos") pairs.
    "Chr" must be in the first column, while "Pos" must be in the second column.
    """

    # Create a local function that splits the rows of the input file.
    ss = lambda line: line.strip("\n").split()
    # Initialise an empty set to fill in with sites.
    sites = set()

    with opener(fin) as f:
        for line in f:
            crm, pos = ss(line)[0:2]
            site = (str(crm), int(pos))
            sites.add(site)

    return sites

def sort_loci_dict(loci_dict, intersect_sites: set=None):
    """
    Sort a loci dictionary by "Chr" and "Pos" keys.
    """

    # Sort fields by "Chr", then "Pos", then "Maj_allele", etc.
    sorted_fields = sorted(zip(
        loci_dict["Chr"], loci_dict["Pos"],
        loci_dict["Maj_allele"], loci_dict["Min_allele"],
        loci_dict["Maj_count"], loci_dict["Min_count"]))
    # Remove sites which are not found in `intersect_sites`
    if intersect_sites:
        sorted_fields = list(filter(
            lambda i: (i[0], i[1]) in intersect_sites,
            sorted_fields))
    # Recombine tuples into a dictionary of lists.
    sorted_loci = {"Chr": [f[0] for f in sorted_fields],
            "Pos": [f[1] for f in sorted_fields],
            "Maj_allele": [f[2] for f in sorted_fields],
            "Min_allele": [f[3] for f in sorted_fields],
            "Maj_count": [f[4] for f in sorted_fields],
            "Min_count": [f[5] for f in sorted_fields], }

    return sorted_loci

def from_sync_to_loci_dict(fin: str, nchr: int):
    """
    Read a file "*sync" (possibly gzipped) and convert it to a python
    dictionary.

    Input
    -----

    fin : str
    The path/filename to a "*sync" file, either uncompressed or gzipped.

    nchr : int
    The number of sets of sequenced chromosomes (depth). It is twice the amount
    of sequenced individuals for diploid organisms.

    Output
    ------

    A dictionary which stores "chr", "pos", "alleles" and "allele counts".
    """

    # Initialise a function to parse the "sync string" into allele counts.
    from_sync_to_allele_count = lambda sync_string: sorted((
        (nuc, int(count)) for count, nuc in zip(
            sync_string.split(":")[0:4], ("A", "T", "C", "G"))
        if int(count) > 0), key=lambda t: t[1])
    # Initialise the output dictionary. Monomorphisms will be stored as
    # "Min_allele": "N", "Min_count": 0 (as if the minor allele were ambiguous
    # with an allele count of zero).
    loci = {"Chr": [], "Pos": [], "Maj_allele": [], "Min_allele": [],
            # Counts of allele1 (NA1, major) and allele2 (NA2, minor).
            "Maj_count": [], "Min_count": [], }
    # Create a local function that splits the rows of the input file.
    ss = lambda line: line.strip("\n").split()
    # Create a function which categorises the fields of a given line.
    to_fields = lambda line: {
        "Chr": str(ss(line)[0]), "Pos": int(ss(line)[1]),
        "Allele_counts": from_sync_to_allele_count(ss(line)[3]), }
    # Initialise the `sfs_idx`.
    sfs_idx = to_sfs_indices(nchr)

    with opener(fin) as f:
        for line in f:
            # First of all, check that the reference is unambigous. Otherwise
            # the `sync_string` will not be parseable, i.e. ".:.:.:.".
            if str(ss(line)[2]) not in ("A", "T", "C", "G"):
                continue
            # Parse its fields.
            fields = to_fields(line)
            # If the site is triallelic (more than 2 allele counts), or the site
            # has missing data (no allele counts), skip the site.
            if len(fields["Allele_counts"]) not in (1, 2):
                continue
            elif len(fields["Allele_counts"]) == 1:
                nuc, count = fields["Allele_counts"][0]
                loci["Maj_allele"].append(nuc)
                loci["Maj_count"].append(nchr)
                loci["Min_allele"].append("N")
                loci["Min_count"].append(0)
            elif len(fields["Allele_counts"]) == 2:
                (maj_nuc, maj_count), (min_nuc, min_count) = fields["Allele_counts"]
                min_freq = min_count / (sum([maj_count, min_count]))
                minor_all_count = to_allele_count(
                    to_abs_diff(min_freq, sfs_idx))
                loci["Maj_allele"].append(maj_nuc)
                loci["Maj_count"].append(int(minor_all_count))
                loci["Min_allele"].append(min_nuc)
                loci["Min_count"].append(int(nchr - minor_all_count))
            # Finally, add the "Chr" and "Pos".
            loci["Chr"].append(fields["Chr"])
            loci["Pos"].append(fields["Pos"])

    return sort_loci_dict(loci)

def from_ngsPool_out_to_loci_dict(fin: str, nchr: int,
                                  lrt_snp_thresh: float=6.64):
    """
    Read a file "*ngsPool.out" (possibly gzipped) and convert it to a python
    dictionary.

    Input
    -----

    fin : str
    The path/filename to an "*.ngsPool.out" file, either uncompressed or
    gzipped.

    nchr : int
    The number of sets of chromosomes. It is twice the amount of individuals for
    diploid organisms.

    lrt_snp_thresh : float
    Threshold to call SNPs. See `ngsJulia` documentation.

    Output
    ------

    A dictionary which stores "chr", "pos", "alleles" and "allele counts".
    """

    # Initialise the output dictionary. Monomorphisms will be stored as
    # "Min_allele": "N", "Min_count": 0.
    loci = {"Chr": [], "Pos": [], "Maj_allele": [], "Min_allele": [],
            # Counts of allele1 (NA1, major) and allele2 (NA2, minor).
            "Maj_count": [], "Min_count": [], }
    # Create a local function that splits the rows of the input file.
    ss = lambda line: line.strip("\n").split()
    # Create a function which categorises the fields of a given line.
    to_fields = lambda line: {
        "Chr": str(ss(line)[0]), "Pos": int(ss(line)[1]),
        "Maj_allele": str(ss(line)[4]), "Min_allele": str(ss(line)[5]),
        "Minor_saf_mle": float(ss(line)[10]), "lrt_snp": float(ss(line)[6]), }
    # Initialise the `sfs_idx`.
    sfs_idx = to_sfs_indices(nchr)

    with opener(fin) as f:
        # Read the first line, which may be a header.
        line = f.readline()
        # If it is a header, skip it. If it is not a header, append it.
        if not "##chromosome" in line:
            # Parse its fields.
            fields = to_fields(line)
            # Append fields to the lists of a dict.
            loci["Chr"].append(fields["Chr"])
            loci["Pos"].append(fields["Pos"])
            loci["Min_allele"].append(fields["Min_allele"])
            loci["Maj_allele"].append(fields["Maj_allele"])
            # Compute the minor allele count from minor allele freq.
            if fields["lrt_snp"] > lrt_snp_thresh:
                minor_all_count = to_allele_count(
                    to_abs_diff(fields["Minor_saf_mle"], sfs_idx))
            else:
                minor_all_count = 0
            # Set minor and then major allele counts.
            loci["Min_count"].append(int(minor_all_count))
            loci["Maj_count"].append(int(nchr - minor_all_count))

        # Repeat the same procedure for the rest of the file lines after the
        # first one.
        for line in f:
            # Parse its fields.
            fields = to_fields(line)
            # Append fields to the lists of a dict.
            loci["Chr"].append(fields["Chr"])
            loci["Pos"].append(fields["Pos"])
            loci["Min_allele"].append(fields["Min_allele"])
            loci["Maj_allele"].append(fields["Maj_allele"])
            # Compute the minor allele count from minor allele freq.
            if fields["lrt_snp"] > lrt_snp_thresh:
                minor_all_count = to_allele_count(
                    to_abs_diff(fields["Minor_saf_mle"], sfs_idx))
            else:
                minor_all_count = 0
            # Set minor and then major allele counts.
            loci["Min_count"].append(int(minor_all_count))
            loci["Maj_count"].append(int(nchr - minor_all_count))

    return sort_loci_dict(loci)

def shared_sites_loci_dict(*loci_dicts: dict, intersect_sites: set=None):
    """
    Find the intersection between multiple `loci` dictionaries.
    """

    # If there is a single dict, it is trivial to compute shared sites.
    if len(loci_dicts) == 1:
        raise UserWarning("A single `loci_dict` was provided; to find shared"
                          + " sites at least a pair (or more) must be"
                          + " provided.")
    # Local function to get a set of sites from a `loci_dict`.
    get_sites = lambda loci: set({(crm, pos)
                                  for crm, pos in zip(
                                      loci["Chr"], loci["Pos"])})
    # Get the sites from the first `loci_dict`.
    sites = get_sites(loci_dicts[0])
    # Get the intersection between the first `loci_dict` and all the following
    # `loci_dict`. At the end we will obtain an intersection of all files.
    for loci in loci_dicts[1:]:
        # Get the intersection between already intersected sites and next sites.
        sites = sites.intersection(get_sites(loci))

    # If they also inputted some `intersect_sites`, find this intersection.
    if intersect_sites:
        sites = sites.intersection(intersect_sites)

    # Check that there is at least a single match.
    if len(sites) == 0:
        raise UserWarning("There are no matching sites between ALL of the"
                          + " queried files.")
    # Now, filter all of the `loci_dict` using these shared sites.
    answer = list()
    for loci in loci_dicts:
        sorted_and_filtered_loci = sort_loci_dict(loci, sites)
        answer.append(sorted_and_filtered_loci)

    return answer

def equal_length_dict_values(*dicts: dict):
    """
    Given a dictionary where all of its values are lists, check that all lists
    are of the same length.
    """

    agglutinative_set = set()
    for d in dicts:
        agglutinative_set = agglutinative_set.union(
            set([len(l) for l in d.values()]))

    if len(agglutinative_set) == 1: return agglutinative_set.pop()
    else: raise UserWarning("The given files may"
                            + " be of different lengths"
                            + " (number of rows); make sure"
                            + " the sites have one-to-one"
                            + " correspondence.")

def from_loci_dict_to_moments_dict(loci_dicts: list, pop_labels: list):
    """
    Builds a data dict with the `moments` or `dadi` format, which is the
    following:

    data_dict = {
     ('chr', 'pos'): {
        'context': '---',
        'outgroup_context': '---',
        'outgroup_allele': '-',
        'segregating': ('nuc1', 'nuc2'),
        'calls': {'label-pop1': (nuc1, nuc2),
                  'label-pop2': (nuc1, nuc2),
                  ...} }
     ('chr', 'pos'): { ... }
     ...
    }

    "Context", "outgroup context" and "outgroup allele" are left in blank (unset
    or unknown, as in a folded SFS).
    """

    # Initialise the output dictionary.
    moments_dd = dict()
    # Make sure all of the input `dicts` are sorted by ("Chr", "Pos").
    loci_dicts = [sort_loci_dict(d) for d in loci_dicts]

    # For all the loci (sites) in the input dicts...
    for idx_locus in range(equal_length_dict_values(*loci_dicts)):
        # Check that the locus at `idx_locus` is found at the same ("Chr",
        # "Pos") across all dicts.
        sites = set([(d["Chr"][idx_locus], d["Pos"][idx_locus])
                     for d in loci_dicts])
        if len(sites) > 1:
            raise UserWarning("The input files do not have one-to-one"
                              + " correspondence.")
        # Pop out the only site in sites of len == 1.
        else: site = sites.pop()
        # Obtain alleles zipped with their counts.
        alleles = [(d["Min_allele"][idx_locus],
                    d["Min_count"][idx_locus]) for d in loci_dicts]
        alleles.extend([(d["Maj_allele"][idx_locus],
                         d["Maj_count"][idx_locus]) for d in loci_dicts])
        # Filter out alleles with zero counts (these may be errors).
        alleles = list(filter(lambda a: a[1] > 0, alleles))
        alleles_nucleotides = list(set([a[0] for a in alleles]))
        # If there are more than two alleles, this is a triallelic site. Skip
        # it and warn.
        if len(alleles_nucleotides) > 2:
            print("WARNING: Skipping triallelic site at", site,
                  "with inferred alleles", alleles_nucleotides)
            continue

        # Create an entry for the output dictionary.
        locus = {
            # Outgroup information is not available.
            "context": "---",
            "outgroup_context": "---",
            "outgroup_allele": "-",
            "segregating": tuple(alleles_nucleotides)
                if len(alleles_nucleotides) == 2
                else (alleles_nucleotides[0], "N"),
            "calls": dict(), }

        # Append "calls" for each input/population dict.
        for d, label in zip(loci_dicts, pop_labels):
            if d["Maj_allele"][idx_locus] == alleles_nucleotides[0] or \
                    d["Min_allele"][idx_locus] == alleles_nucleotides[1]:
                locus["calls"][label] = (d["Maj_count"][idx_locus],
                                         d["Min_count"][idx_locus])

            elif d["Maj_allele"][idx_locus] == alleles_nucleotides[1] or \
                    d["Min_allele"][idx_locus] == alleles_nucleotides[1]:
                locus["calls"][label] = (d["Maj_count"][idx_locus],
                                         d["Min_count"][idx_locus])

            else:
                print(alleles_nucleotides)
                print(d["Allele2"][i], d["Allele1"][i])
                raise UserWarning("One allele did not match?")
        # Append this locus to `moments_dd`.
        moments_dd[site] = locus

    return moments_dd

def build_moments_data_dict(fins: list, pop_labels: list, nchroms: list,
                            intersect_sites: str=None, ):
    """
    Build a `moments` data dict, a dict which represents a multidimensional SFS.

    Input
    -----

    fins : list
    A list of input filenames.

    pop_labels : list
    A list of population labels which will be zipped with the list of input
    filenames.

    nchroms : list
    A list with the num. of chrs. which will be zipped with the list of input
    filenames. It will define the highest index of the SFS.

    intersect_sites : str
    A filename that directs to a BED-like file with coordinates arranged as
    "chr" in first and "pos" in second columns.
    """

    # Read the provided files in `fins`.
    loci_dicts = list()
    for fin, nchr in zip(fins, nchroms):
        # Read if it is a SYNC file.
        if "sync" in fin:
            print("INFO: Reading the file `"
                  + str(fin) + "` as a SYNC file.")
            loci = from_sync_to_loci_dict(fin, nchr)
            loci_dicts.append(loci)
        # Read if it is an ngsPool file from ngsJulia.
        elif "out" in fin or "ngsPool" in fin:
            print("INFO: Reading the file `"
                  + str(fin) + "` as an ngsPool file from ngsJulia.")
            loci = from_ngsPool_out_to_loci_dict(fin, nchr)
            loci_dicts.append(loci)
        else:
            raise UserWarning("The file `" + str(fin) + "` has neither"
                              + " 'sync', 'out' nor 'ngsPool' extensions.")

    # Read the intersecting files, if provided.
    if intersect_sites:
        print("INFO: Reading the given bed-like coordinates.")
        intersect_sites = read_bedlike_sites(intersect_sites)
        print("INFO: Finding the intersection between files"
              + " and the provided `intersect_sites`.")
    else:
        print("INFO: Finding the intersection between files.")

    # Make sure all of the input files contain the same loci.
    loci_dicts = shared_sites_loci_dict(*loci_dicts,
                                        intersect_sites=intersect_sites)

    # Build the moments data dict.
    print("INFO: Building the `moments` data dict.")
    moments_dd = from_loci_dict_to_moments_dict(loci_dicts, pop_labels)

    return moments_dd

def distance_thin(moments_dd: dict, pop_label: str, unlinked_distance: int):
    """
    Iterate across loci in a moments' "data_dict" to make sure that all
    polymorphisms are further away than "unlinked_distance". If two
    polymorphisms are too close-by, one of them will be removed from the output
    list of loci.

    Input
    -----

    moments_dd : dict
    Output from "from_loci_dict_to_moments_dict".

    pop_label : str
    Population label within the "calls" key of the aforementioned "moments_dd".

    unlinked_distance : int
    The distance presumed to be enough for two loci to be unlinked (with
    independent recombination).

    Output
    ------

    A "moments_dd" with loci at enough distance for them to be assumed unlinked.
    """

    unlinked_moments_dd = dict()
    last_sequid = 0
    last_pos = 0

    for locus, data in moments_dd.items():
        # Find polymorphic sites (both alleles are at freq. > 0).
        if 0 not in data["calls"][pop_label] and (
                # Moreover, either the last and present sites are in different
                # sequids/chrs. and thus they are unlinked...
                last_sequid != locus[0] or
                # ...Otherwise, the distance within the sequid is big enough.
                last_pos + unlinked_distance <= locus[1]):
            # Add the sites which pass both tests.
            unlinked_moments_dd[locus] = data
            # Update last sequid and last pos.
            last_sequid, last_pos = locus

    return unlinked_moments_dd

def remove_monomorphic_sites(moments_dd: dict):
    """
    Removes sites where all populations have a monomorphism.
    """

    polym_moments_dd = dict()

    for locus, data in moments_dd.items():
        # Monomorphisms are composed of a segregating 'ATGC' and an 'N'.
        # So, if it is not a monomorphism, keep it in the filtered dict.
        if "N" not in data["segregating"]:
            polym_moments_dd[locus] = data

    return polym_moments_dd

def bootstrap(moments_dd: dict, n_loci: int, n_pseudoreps: int):
    """
    Create a list of pseudoreplicate `moments_dd` by randomly sampling subsets.

    Input
    -----

    moments_dd : dict
    Output from "from_loci_dict_to_moments_dict".

    n_loci : int
    asd

    n_pseudoreps : int
    asd

    Output
    ------

    A list of pseudoreplicates obtained by subsampling.
    """

    # Initialise a list of data-dicts which will be returned.
    pseudoreplicates = list()
    loci_available = [locus for locus in moments_dd.keys()]
    if len(loci_available) <= n_loci:
        print("WARNING: Amount of loci requested per pseudoreplicate is"
              + " higher than or equal to the amount of loci in the"
              + " original dataset.")
        print("WARNING: Num. of loci in the original dataset:"
              f" {len(loci_available)}")
        print("WARNING: Num. of loci requested per pseudoreplicate:"
              f" {n_loci}")

    for _ in range(n_pseudoreps):
        # Create each instance of a pseudoreplicate.
        pseudorep_moments_dd = dict()
        loci_sampled_with_replacement = [
            # Randomly draw as many as "n_loci" loci.
            random.choice(loci_available) for _ in range(n_loci)]
        # Due to the fact that the keys of a dict. cannot be duplicated, and the
        # same site may be sampled more than once, assign the loci an arbitrary "i"
        # key.
        for i in range(len(loci_sampled_with_replacement)):
            seq, pos = loci_sampled_with_replacement[i]
            # Keys are made up of "seq" and "pos", but also the arbitrary "i".
            pseudorep_moments_dd[(seq, pos, i)] = moments_dd[seq, pos]
        # Append this pseudoreplicate to a list.
        pseudoreplicates.append(pseudorep_moments_dd)

    return pseudoreplicates


if __name__ == '__main__':
    # This module might be easier to test and debug through the script
    # '__main__.py'.
    import pprint

    # Input strings must be formatted as "label=nchr=file-input".
    pop_labels = [arg.split("=")[0] for arg in sys.argv[1:]]
    nchroms = [int(arg.split("=")[1]) for arg in sys.argv[1:]]
    fins = [arg.split("=")[2] for arg in sys.argv[1:]]

    moments_dd = build_moments_data_dict(fins, pop_labels, nchroms)
    pprint.pprint(moments_dd)

    bootstrapped_moments_dd = bootstrap(moments_dd,
                                        n_loci=int(len(moments_dd)/2),
                                        n_pseudoreps=3)
    pprint.pprint(bootstrapped_moments_dd)

