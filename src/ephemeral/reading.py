#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# reading.py (v3)
#
# 20 de des. 2024  <adria@molevol-OptiPlex-9020>

help_msg = """
Module responsible for reading input files with called polymorphisms from
pooled sequencing data with the goal of creating a `moments.Spectrum` --a Site
Frequency Spectrum implementation in Python.

Available formats for conversion to `moments.Spectrum`:
- *sync
- *ngsPool.out (output from ngsJulia).
"""

# Prepare for reading compressed (gzipped) input files.
import gzip
# Create sqlite3 databases as temporary data files in the home directory.
import sqlite3
import tempfile
import pathlib
# Sanitize sqlite3 input with regex instead of using placeholders.
import re
# To assign a name to the sqlite3 table created with the BED intersection sites,
# we will use its filename.
import os.path
# Logging.
import logging

# Conversion from data dictionaries to `moments.Spectrum`.
import moments.Spectrum_mod


# Prepare logging.
logger = logging.getLogger(__name__)

# Decide if a file is compressed or uncompressed and open it with the correct
# function depending on its extension.
opener = lambda finput: (
    gzip.open(finput, "rt") if finput.endswith(".gz")
    else open(finput, "rt"))
# Strip and split lines read from a file.
ss = lambda line, sep=None: line.strip("\n").split(sep)
# Sanitize sqlite3 input, because table names cannot be created with
# placeholders, which could allow sqlite injection attacks.
# If the regex doesn't match, raise Exception.
def sanitized_sqlite(parameter: str):
    allowed_regex = r"[\+a-zA-Z0-9_-]+"
    if not re.fullmatch(allowed_regex, parameter):
        raise Exception(f"The parameter '{parameter}' is not"
                        + " sanitized for it to be used with"
                        + " sqlite3; it must be a string"
                        + f" matching '{allowed_regex}'.")
# Parse a "sync string" (its fourth column) into a list of pairs of alleles
# (nucleotides) and allele counts.
def from_syncstring_to_allele_depths(sync_string):
    return sorted(((nuc, int(count)) for count, nuc in zip(
            sync_string.split(":")[0:4], ("A", "T", "C", "G"))
        # Do not call/include if depth is zero.
        if int(count) > 0),
        # Sorting key is allele count.
        key=lambda tup: tup[1])


class ParserGenData:
    """
    Parse genomic data. Converts a genomic data file to a python iterable.
    """

    def __init__(self, filename, func_generator, nchr=None,
                 lrt_snp_thresh=None):

        self.filename = filename
        self._func_generator = func_generator
        if nchr:
            self.nchr = nchr
            # Create a list of [0/nchr, 1/nchr, etc. nchr/nchr], which is a
            # discrete distribution with "nchr + 1" bins of the possible allele
            # frequencies in an SFS with "nchr" sequences/chromosomes/haploid
            # samples. It will help in translating continuous allele frequencies
            # to discrete allele counts.
            self.sfs_frequencies = [i / nchr for i in range(0, nchr + 1)]
        if lrt_snp_thresh:
            self.lrt_snp_thresh = lrt_snp_thresh

        return None

    def __iter__(self):
        return iter(self._func_generator(self))

    def saf_to_allele_count(self, saf):
        """
        Compute the absolute differences between a site allele frequency "saf"
        and each index within the list "sfs_frequencies".
        Return a discrete allele count from the starting continuous allele
        frequency, which is the minimum in a list of absolute differences.
        """

        abs_diff_list = [abs(saf - i) for i in self.sfs_frequencies]
        return abs_diff_list.index(min(abs_diff_list))

    @classmethod
    def from_bed(cls, filename):
        """
        Read a BED file and yield pairs of "chromosome" and "site" per row.

        Input
        -----

        filename : str
        File input name. Comments and headers within the BED file must start
        with "#".
        """

        def generator(self):
            # Reset counters after each call to the iterable.
            self.skipped_comment = 0

            with opener(self.filename) as fhandle:
                for numline, line in enumerate(fhandle):
                    # Make sure this line is not a header or a comment.
                    if "#" in line[0]:
                        self.skipped_comment += 1
                        continue
                    line = ss(line)
                    # Slice the line by the column indices of a BED file.
                    chr_column, pos_column = 0, 1
                    try: chromosome, pos = (str(line[chr_column]),
                                            int(line[pos_column]))
                    except ValueError as err:
                        raise err(f"'{pos}' was retrieved as a chromosome"
                                  + f" position from '{finput}' at line"
                                  + f" num. '{numline + 1}'. It could"
                                  + " not be cast as an integer type;"
                                  + " check if the correct column indice"
                                  + f" for this data is '{pos_column}' and"
                                  + " also if there are comments or headers"
                                  + " which do not start with '#' as their"
                                  + " first character.")
                    yield {
                        "chr": chromosome,
                        "pos": pos, }

            # Store num. of lines (sites in data) after reaching EOF.
            self.numsites = numline + 1 - self.skipped_comment

        return cls(filename, generator)

    @classmethod
    def from_sync(cls, filename, nchr):
        """
        Read a SYNC file with POOLED polymorphism data. It can be gzipped or
        uncompressed. Then, yields a data dictionary for each row.

        Input
        -----

        filename : str
        File input name. Comments and headers within the SYNC file must start
        with "#".

        nchr : int
        The number of sets of chromosomes in the pool (sequencing depth). It is
        twice the amount of sequenced individuals for diploid organisms.
        """

        def generator(self):
            # Reset counters after each call to the iterable.
            self.skipped_comment = 0
            self.skipped_triallelic = 0
            self.skipped_missing_data = 0
            self.skipped_ambiguous_ref = 0

            # Start by creating a function which categorises fields.
            def sync_to_datadict(line):
                # Split and slice line.
                ssline = ss(line)
                # First, check that the reference is unambiguous. Otherwise
                # the "sync string" will not be parseable, i.e. ".:.:.:.".
                if str(ssline[2]) not in ("A", "T", "C", "G"):
                    self.skipped_ambiguous_ref += 1
                    return None
                # Locate data by column indice.
                fields = {
                    "chr": str(ssline[0]), "pos": int(ssline[1]),
                    "allele_depths": from_syncstring_to_allele_depths(
                        ssline[3]),
                    }
                # If the site is triallelic (more than 2 allele counts), or the
                # site has missing data (no allele counts), skip the site.
                if len(fields["allele_depths"]) > 2:
                    self.skipped_triallelic += 1
                    return None
                # Missing data, all alleles at zero depth.
                elif len(fields["allele_depths"]) == 0:
                    self.skipped_missing_data += 1
                    return None
                # If there is one allele, yield it as a dictionary.
                elif len(fields["allele_depths"]) == 1:
                    nuc, depth = fields["allele_depths"][0]
                    return {
                        "chr": fields["chr"],
                        "pos": fields["pos"],
                        "maj_allele": nuc,
                        "maj_depth": depth,
                        "maj_freq": 1,
                        "maj_count": self.nchr,
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
                    min_count = self.saf_to_allele_count(min_freq)
                    # For each row in the file "finput", yield a dictionary with
                    # the data in this row.
                    return {
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

            # Open the input file.
            with opener(self.filename) as fhandle:
                # Read the first lines, which may be comments or a header.
                line = fhandle.readline()
                while "#" in line[0]:
                    self.skipped_comment += 1
                    line = fhandle.readline()
                # After the first "#" in lines, start reading data.
                # Data dictionaries "dd" can be "None" because of missing data;
                # skip them, do not yield.
                dd = sync_to_datadict(line)
                if dd: yield dd
                for numline, line in enumerate(fhandle):
                    dd = sync_to_datadict(line)
                    if dd: yield dd

            # "enumerate" starts at 0, and we process the first line before the
            # "enumerate" for-loop; we need to add 2 to obtain num. of sites.
            self.numsites = numline + 2

        return cls(filename, generator, nchr=nchr)

    @classmethod
    def from_ngsPool(cls, filename, nchr, lrt_snp_thresh: float=6.64):
        """
        Read an ngsPool (from ngsJulia) file with POOLED polymorphism data.
        It can be gzipped or uncompressed. Then, yields a data dictionary for
        each row.

        Input
        -----

        filename : str
        File input name. Comments and headers within the SYNC file must start
        with "#".

        nchr : int
        The number of sets of chromosomes in the pool (sequencing depth). It is
        twice the amount of sequenced individuals for diploid organisms.

        lrt_snp_thresh : float
        Threshold to call SNPs. See `ngsJulia` documentation.
        """

        def generator(self):
            # Reset counters after each call to the iterable.
            self.skipped_comment = 0
            self.skipped_below_lrt_snp = 0

            # Start by creating a function which categorises fields.
            def ngsjulia_to_datadict(line):
                # Split and slice line.
                ssline = ss(line)
                # Locate data by column indice.
                fields = {
                    "chr": str(ssline[0]), "pos": int(ssline[1]),
                    "maj_allele": str(ssline[4]),
                    "min_allele": str(ssline[5]),
                    "min_saf_mle": float(ssline[10]),
                    "lrt_snp": float(ssline[6]), }
                # If the Likelihood-Ratio-Test (LRT) for the SNP is higher than
                # the threshold, accept the hypothesis that this loci is
                # polymorphic (see ngsJulia docs).
                if fields["lrt_snp"] > self.lrt_snp_thresh:
                    min_freq = fields["min_saf_mle"]
                    min_count = self.saf_to_allele_count(min_freq)
                    min_nuc = fields["min_allele"]
                # Otherwise, below the threshold the loci is monomorphic.
                else:
                    # Counter of loci with maximum likelihood estimate "mle" of
                    # site allele frequency "saf" above zero which have been
                    # rejected by lrt threshold.
                    if fields["min_saf_mle"] > 0.0:
                        self.skipped_below_lrt_snp += 1
                    min_freq = 0
                    min_count = 0
                    min_nuc = "N"
                # Return the datadict created from the ngsjulia line.
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

            # Open the input file.
            with opener(self.filename) as fhandle:
                # Read the first lines, which may be comments or a header.
                line = fhandle.readline()
                while "#" in line[0]:
                    self.skipped_comment += 1
                    line = fhandle.readline()
                # After the first "#" in lines, start reading data.
                yield ngsjulia_to_datadict(line)
                for numline, line in enumerate(fhandle):
                    yield ngsjulia_to_datadict(line)

            # "enumerate" starts at 0, and we process the first line before the
            # "enumerate" for-loop; we need to add 2 to obtain num. of sites.
            self.numsites = numline + 2

        return cls(filename, generator, nchr=nchr,
                   lrt_snp_thresh=lrt_snp_thresh)

    def populate_new_sqlite_table(self, tblname: str, db_sqlite: str):
        """
        Create a new sqlite3 table and populate it with the values parsed from
        this instance.

        Input
        -----

        tblname : str
        Name of the newly created table within the sqlite database file.
    
        db_sqlite : str
        File name of the sqlite database which this function will write to.
        """

        # Open a connection and a cursor to the sqlite database.
        con = sqlite3.connect(db_sqlite)
        cur = con.cursor()
        generator = iter(self)
        # Take a look at the fields/columns returned by the generator.
        first_row = next(generator)
        # Make sure they include required "chromosome" and "pos" fields.
        if "chr" not in first_row.keys() or "pos" not in first_row.keys():
            raise Exception("'Generator' from input file does not have"
                            + " 'chr' or 'pos' fields.")
        # Try to avoid SQL injections by sanitizing input.
        [sanitized_sqlite(param) for param in (tblname, *first_row.keys())]
        # Create a table with the name "tblname" and the fields/columns in
        # "first_row.keys()".
        cur.execute(f"CREATE TABLE {tblname} ("
                    + ", ".join(first_row.keys()) + ", "
                    + " PRIMARY KEY(chr, pos) )")
        # Create a command to insert data in the sqlite table.
        command_insert_sqlite = str(
            f"INSERT INTO {tblname} VALUES ("
            + ", ".join([":" + key for key in first_row.keys()]) + ")")
        # Add the data in "first_row" to the newly created table.
        cur.execute(command_insert_sqlite, first_row)
        # Add the data of the rest of the generator.
        for row in generator:
            logger.debug(row)
            cur.execute(command_insert_sqlite, row)
        # Commit these INSERT changes.
        con.commit()
        con.close()
    
        return None

def select_shared_loci(db_sqlite: str):
    """
    Open an sqlite3 database and select the inner join (shared or intersection
    of sites) of all of its tables.

    Input
    -----

    db_sqlite : str
    File name of the sqlite database which this function will read from.
    """

    # Open a connection and a cursor to the sqlite database.
    con = sqlite3.connect(db_sqlite)
    cur = con.cursor()
    # Get a list of all the tables within "db_sqlite".
    db_tblnames = cur.execute(
        "SELECT name FROM sqlite_master WHERE type='table'")
    # It returns an sqlite.Cursor (iterable of tuples with a single item "name")
    # which must be cast into a list.
    db_tblnames = [tup[0] for tup in db_tblnames]
    # Select an inner join of all the tables, equivalent to sites/loci shared
    # across all files/tables.
    command_intersection_sqlite = str(
        "SELECT ROW_NUMBER() OVER(ORDER BY chr, pos), *"
        + f" FROM {db_tblnames[0]} "
        + " ".join(["INNER JOIN " + t + " USING (chr, pos)"
                    for t in db_tblnames[1:]])
    )
    # Execute the command.
    # The selected rows with polymorphism data are stored in "res".
    db_res = cur.execute(command_intersection_sqlite)

    return db_res, db_tblnames, con

def select_median(colname: str, tblname: str, db_sqlite: str):
    """
    Returns the median of a numeric column in a sqlite3 database.
    """

    # Open a connection and a cursor to the sqlite database.
    con = sqlite3.connect(db_sqlite)
    cur = con.cursor()
    # Try to avoid SQL injections by sanitizing input.
    [sanitized_sqlite(param) for param in (tblname, colname)]

    command_median = str(
        "SELECT AVG(x) FROM"
        + f" (SELECT {colname} AS x FROM {tblname} ORDER BY x"
        + f" LIMIT 2 - (SELECT COUNT(*) FROM {tblname}) % 2"
        + f" OFFSET (SELECT (COUNT(*) - 1) / 2 FROM {tblname}))"
    )
    median = cur.execute(command_median)
    con.close()

    # "median" is an iterable of tuples (database query); return only the first
    # item of the first tuple (a number).
    return float(median.fetchone()[0])

def select_average(colname: str, tblname: str, db_sqlite: str):
    """
    Returns the average of a numeric column in a sqlite3 database.
    """

    # Open a connection and a cursor to the sqlite database.
    con = sqlite3.connect(db_sqlite)
    cur = con.cursor()
    # Try to avoid SQL injections by sanitizing input.
    [sanitized_sqlite(param) for param in (tblname, colname)]

    command_average = str(f"SELECT AVG({colname}) FROM {tblname}")
    average = cur.execute(command_average)
    con.close()

    # "average" is an iterable of tuples (database query); return only the first
    # item of the first tuple (a number).
    return float(average.fetchone()[0])

def select_average_expected_heterozygosity(tblname: str, db_sqlite: str,
                                           colname_minfreq: str="min_freq"):
    """
    Returns the average across all sites of the expected heterozygosity
    (as in Nei's genetic diversity estimator).
    """

    # Open a connection and a cursor to the sqlite database.
    con = sqlite3.connect(db_sqlite)
    cur = con.cursor()
    # Try to avoid SQL injections by sanitizing input.
    [sanitized_sqlite(param) for param in (tblname, colname_minfreq)]

    command_average_expected_heterozygosity = str(
        "SELECT SUM(heteroz) / COUNT(heteroz) FROM ("
            + "SELECT 1 - (x * x) - ((1 - x) * (1 - x)) AS heteroz FROM ("
                + f"SELECT {colname_minfreq} AS x FROM {tblname}))"
    )
    avg_exp_hetero = cur.execute(command_average_expected_heterozygosity)
    con.close()

    # "avg_exp_hetero" is an iterable of tuples (database query); return only
    # the first item of the first tuple (a number).
    return float(avg_exp_hetero.fetchone()[0])

def from_dbrow_to_moments_dd(db_row, pop_labels):
    """
    Convert a row from a database created with the method
    "ParserGenData.populate_new_sqlite_table" into a moments data dictionary.
    """

    # (chr, pos) in the second and third fields of "db_row".
    site = tuple(db_row[1:3])
    # Compute and store the length of "db_row".
    len_dbrow = len(db_row)
    # Major allele in the fourth field "3", which repeats every eight (for the
    # rest of the populations).
    maj_alleles = tuple([str(db_row[i]) for i in range(3, len_dbrow, 8)])
    # Minor allele in the eighth field "7", and repeats every eight.
    min_alleles = tuple([str(db_row[i]) for i in range(7, len_dbrow, 8)])
    # Obtain all of the observed alleles at this site.
    alleles = list(set([a for a in maj_alleles + min_alleles
                        if a != "N"]))
    # If the site is triallelic, skip it (return None).
    if len(alleles) > 2: return site, None
    # If the site is monomorphic, add an ambiguous nucleotide for formatting
    # reasons.
    if len(alleles) == 1:
        alleles = (alleles[0], "N")
    # Major allele counts in the seventh field "6", and repeats every eight.
    maj_counts = tuple([int(db_row[i]) for i in range(6, len_dbrow, 8)])
    # Minor allele counts in the eleventh field "10", and repeats every eight.
    min_counts = tuple([int(db_row[i]) for i in range(10, len_dbrow, 8)])
    # Compute calls for each population.
    calls = {
        pop_lab:
        (maj_count, min_count) if maj_allele == alleles[0]
        else (min_count, maj_count)
        for pop_lab, maj_allele, maj_count, min_count in zip(
                pop_labels, maj_alleles, maj_counts, min_counts)
        }
    data_dict = {
            "context": "---",
            "outgroup_context": "---",
            "outgroup_allele": "-",
            "segregating": alleles,
            "calls": calls, }

    return site, data_dict

def combinations_populations(pop_labels: list):
    """
    Compute possible combinations of a single or a pair of populations, in order
    to compute all possible 1D and 2D SFS.
    """

    def combine(arr):
        if len(arr) == 0: return [[]]
        combs = []
        for c in combine(arr[1:]):
            combs += [c, c + [arr[0]]]
        return combs
    # Only combinations of one or two items (pops) allowed.
    pop_combinations = [c for c in combine(pop_labels)
                        if len(c) == 1 or len(c) == 2]

    return {tuple(combo): list() for combo in pop_combinations}

def from_moments_dd_to_sfs_dict(moments_dd, sfs_dict, poplab_to_nchrom):
    """
    Modifies in place "sfs_dict" by adding information from sites in
    "moments_dd".
    """

    for poplabs in sfs_dict.keys():
        # Now, "moments_dd" has all of the provided population calling data.
        # Filter the pop. ids. we are interested in.
        filtered_moments_dd = dict()
        for site, data in moments_dd.items():
            filtered_calls = {pl: calls
                              for pl, calls in data["calls"].items()
                              if pl in poplabs}

            filtered_moments_dd[site] = data.copy()
            filtered_moments_dd[site]["calls"] = filtered_calls

        # From the filtered "moments_dd" create an SFS and append it to the
        # "sfs_dict" dictionary with the key "poplabs".
        sfs_dict[poplabs].append(
            moments.Spectrum.from_data_dict(
                filtered_moments_dd, pop_ids=poplabs,
                projections=[int(poplab_to_nchrom[pl]) for pl in poplabs],
                mask_corners=False, polarized=False))
        # Summing the sites of a pair of SFS into a single SFS.
        # It reduces the amount of data in memory.
        sfs_dict[poplabs] = [sum(sfs_dict[poplabs])]

    return None

def write_stats(parser, stats_dict, poplabel, db_sqlite):
    """
    parser : instance of ParserGenData
    """

    if not stats_dict:
        # Initialize a dictionary to store "stats_dict" regarding the analysis.
        stats_dict = {
            "poplabel": list(),
            "sites_total": list(),
            # Included sites are all lines with data minus skipped lines, which
            # are "None" for the ngsPool and BED parsers; only applicable to
            # sync files.
            "sites_included": list(),
            "skipped_comment": list(),
            "sites_skipped_ambiguous_ref": list(),
            "sites_skipped_missing_data": list(),
            "sites_skipped_triallelic": list(),
            # ngsPool does not represent "depth", only "frequency" directly.
            "median_depth": list(),
            "mean_depth": list(),
            "average_expected_heterozygosity": list(),
        }

    # Append stats_dict to these lists.
    stats_dict["poplabel"].append(poplabel)
    stats_dict["sites_total"].append(parser.numsites)
    # The sync parser has more stats_dict than the ngspool, which does not
    # consider triallelic or missing data sites.
    try:
        stats_dict["sites_included"].append(
                parser.numsites - (parser.skipped_ambiguous_ref
                    + parser.skipped_missing_data + parser.skipped_triallelic)
            )
    except: stats_dict["sites_included"].append(parser.numsites)
    stats_dict["skipped_comment"].append(parser.skipped_comment)
    try: stats_dict["sites_skipped_ambiguous_ref"].append(
            parser.skipped_ambiguous_ref)
    except: stats_dict["sites_skipped_ambiguous_ref"].append(None)
    try: stats_dict["sites_skipped_missing_data"].append(
        parser.skipped_missing_data)
    except: stats_dict["sites_skipped_missing_data"].append(None)
    try: stats_dict["sites_skipped_triallelic"].append(
        parser.skipped_triallelic)
    except: stats_dict["sites_skipped_triallelic"].append(None)
    try:
        stats_dict["median_depth"].append(
            select_median("maj_depth+min_depth", poplabel, db_sqlite))
    except: stats_dict["median_depth"].append(None)
    try:
        stats_dict["mean_depth"].append(
            select_average("maj_depth+min_depth", poplabel, db_sqlite))
    except: stats_dict["mean_depth"].append(None)
    try:
        stats_dict["average_expected_heterozygosity"].append(
            select_average_expected_heterozygosity(poplabel, db_sqlite))
    except: stats_dict["average_expected_heterozygosity"].append(None)

    return stats_dict

def main(finputs: list, pop_labels:list, nchroms: list,
         intersect_sites: str=None):
    """
    Create multiple "moments.Spectrum" from a list of input files with variant
    call data.
    """

    # Associate pop. labels with their num. of chr.
    poplab_to_nchrom = dict()
    for pl, nc in zip(pop_labels, nchroms):
        poplab_to_nchrom[str(pl)] = int(nc)

    # Initialize combinations of pop. labels. An SFS will be computed for each
    # single population and pair of populations.
    sfs_dict = combinations_populations(pop_labels)
    amount_sfs_onedim = len([k for k in sfs_dict.keys() if len(k) == 1])
    amount_sfs_bidim  = len([k for k in sfs_dict.keys() if len(k) == 2])
    logger.info("As many as"
                + " '{}' 1D-SFS and '{}' 2D-SFS will be created.".format(
                amount_sfs_onedim, amount_sfs_bidim))

    # Tell which files will be read as SYNC or ngsPool. Raise an Exception if
    # any extension does not match the expected formats.
    for fi in finputs:
        if "sync" in str(fi).lower():
            logger.info(f"Detected that '{fi}' is a SYNC-formatted file.")
        elif "out" in str(fi).lower() or "ngspool" in str(fi).lower():
            logger.info(f"Detected that '{fi}' is an ngsPool-formatted"
                  + " file.")
        else:
            raise Exception("Could not detect the correct file"
                            + f" extension for '{fi}'; accepts"
                            + " 'sync' for SYNC files and"
                            + " 'out' or 'ngsPool' for ngsPool files"
                            + " (either lower or upper case).")
    if intersect_sites:
        logger.info(f"Using '{intersect_sites}' as an intersecting BED.")

    # Check that the three lists in parameters are of the same length; each
    # input filename must also have a pop. label and a n. of chrom.
    if len(finputs) != len(pop_labels) or len(finputs) != len(nchroms):
        raise Exception("The lists of file inputs, pop. labels and num."
                        + " of chromosomes must be of the same length.")

    # Create a temporary file where an sqlite3 database will be stored.
    # Avoid using the directory/partition "/tmp", because the database file
    # could become too big to fit within there. Instead, place this file in the
    # home directory.
    home = pathlib.Path.home()
    temp_file = tempfile.NamedTemporaryFile(
        dir=home, prefix="ephemeral.tmp", suffix=".sqlite3")
    logger.info("Storing the polymorphism data in a temporary"
          + f" sqlite3 database file at '{temp_file.name}'.")

    # Initialize a "stats" dictionary for storing num. of sites and other
    # information of the input files.
    stats = dict()

    # Create a table within the database for each input file.
    for fi, nc, pl in zip(finputs, nchroms, pop_labels):
        logger.info(f"Reading '{fi}' and writing data to database.")
        # Try to match a "sync", "ngspool" or "out" extension.

        if "sync" in str(fi).lower():
            parser = ParserGenData.from_sync(fi, nc)
            parser.populate_new_sqlite_table(pl, temp_file.name)
            # Compute stats of this data file.
            stats = write_stats(parser, stats, pl, temp_file.name)
            # Print these stats.
            for key, vals in stats.items():
                print("        * {}: {}".format(key, vals[-1]))

        elif "out" in str(fi).lower() or "ngspool" in str(fi).lower():
            parser = ParserGenData.from_ngsPool(fi, nc)
            parser.populate_new_sqlite_table(pl, temp_file.name)
            # Compute stats of this data file.
            stats = write_stats(parser, stats, pl, temp_file.name)
            # Print these stats.
            for key, vals in stats.items():
                print("        * {}: {}".format(key, vals[-1]))

    if intersect_sites:
        logger.info(f"Reading the BED '{intersect_sites}'"
              + " with intersecting sites.")
        parser = ParserGenData.from_bed(intersect_sites)
        # Remove directory and extensions from the filename so it can be used as
        # a table name.
        intersect_tblname = os.path.split(intersect_sites)[1].split(".")[0]
        parser.populate_new_sqlite_table(intersect_tblname, temp_file.name)
        # Compute stats of this data file.
        stats = write_stats(parser, stats, "INTERSECT.BED", None)
        # Print these stats.
        for key, vals in stats.items():
            print("        * {}: {}".format(key, vals[-1]))

    logger.info("Finished reading input data.")

    logger.info("Executing the inner join of the database with the columns"
          + " 'chromosome' and 'pos' (finding sites shared across"
          + " all of the given files).")
    db_res, db_tblnames, con = select_shared_loci(temp_file.name)

    # Remove the "intersect_tblname" from "db_tblnames".
    if intersect_sites:
        db_tblnames = db_tblnames[:-1]
    if list(pop_labels) != list(db_tblnames):
        logger.warning("Unexpected order for the columns of the database.")
        logger.warning(pop_labels, db_tblnames)

    logger.info("Creating 'moments' SFS from the inner join of the"
          + " input polymorphism calling data.")

    # Initialize vars before loop.
    moments_dd = dict()
    skipped_triallelic = 0
    row_num = 0
    # Start iterating through sites found within the intersection (iterate
    # across sites shared across all files).
    for db_row in db_res:
        site, data_dict = from_dbrow_to_moments_dd(db_row, pop_labels)
        if not data_dict:
            logger.warning("Triallelic site in the shared sites.")
            skipped_triallelic += 1
        else:
            moments_dd[site] = data_dict
        # Avoid loading into memory all of the sites at once. When the loop
        # reaches the 100 000th site, convert this data into SFS. Repeat the
        # process and sum the vectors or matrices (1D or 2D SFS) converted at
        # each step.
        row_num = db_row[0]
        if row_num % 100000 == 0:
            logger.info(f"Reached {row_num} sites.")

            # Compute SFS for the keys (pair or single pop.) in
            # "sfs_dict.keys()".
            from_moments_dd_to_sfs_dict(moments_dd, sfs_dict, poplab_to_nchrom)

    if row_num == 0:
        # The length of 'db_res' (i.e. shared sites) is zero. The intersection
        # of ALL files is an empty set.
        raise ValueError("The intersection of ALL files is an empty set;"
                         + " there are no shared sites across all files.")

    con.close()

    # Finally, add the remaining sites in "moments_dd" if the amount of sites
    # was not exactly modulo 100 000:
    logger.info(f"Reached the end at {row_num} sites.")
    if moments_dd:
        from_moments_dd_to_sfs_dict(moments_dd, sfs_dict, poplab_to_nchrom)

    # Append stats.
    for key in stats.keys():
        if key == "poplabel":
            stats[key].append("DB-INNER-JOIN")
        elif key == "sites_total":
            stats[key].append(row_num)
        elif key == "sites_skipped_triallelic":
            stats[key].append(skipped_triallelic)
        elif key == "sites_included":
            stats[key].append(row_num - skipped_triallelic)
        else:
            stats[key].append(None)
    # Print these stats.
    for key, vals in stats.items():
        print("        * {}: {}".format(key, vals[-1]))

    # Return all of the SFS. Remember to remove the temporary list the SFS are
    # in (i.e. slice and return the first item).
    return {key: val[0] for key, val in sfs_dict.items()}, stats

