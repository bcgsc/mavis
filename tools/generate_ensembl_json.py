#!/usr/bin/env python
"""
Retrieve gene/transcript/exon/domain information from Ensembl and build a set
of annotations in custom json format for use with MAVIS. Uses 'pyensembl'
module to retrieve gene, transcript, and exon information from the Ensembl
REST-API. The pyensembl module does not have a method for retrieving domain
information, so manually retrieve domain information for each transcript from
the API, and cache it along with the cache files from pyensembl.
"""
import argparse
import csv
import os
import sys
import time
from collections import defaultdict
from functools import wraps
from glob import glob

import requests
import simplejson as json
from pyensembl import EnsemblRelease


VERSION = "1.0.0"
SCRIPT = os.path.abspath(__file__)
CACHE_DEFAULT = os.environ["HOME"] + "/.cache"

DOMAIN_CACHE = {}
DOMAIN_CACHE_PATH = None


def get_date(date_str=None):
    """
    Returns the date (str) as <weekday> <month> <day> <time hh:mm:ss> <year>
    """
    if date_str:
        return time.asctime(time.localtime(date_str))
    else:
        return time.asctime(time.localtime(time.time()))


def parse_best_file(path):
    """
    Method parse a tab-seperate file into a set of best transcripts.

    Args:
        path (str): path to a tab-seperated list with column: transcript_id

    Returns:
        set: a set of 'best transcript' ids
    """
    best = set()
    with open(path, "r") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            best.add(row["transcript_id"])

    return best


def parse_hugo_file(path):
    """
    Method parse specific columns in the tab-seperated hugo file into a dict.

    Args:
        path (str): path to a tab-seperated list with columns: ensid, hugo

    Returns:
        dict: a dict of ensembl gene ids and associated hugo aliases for the gene
    """
    hugo = defaultdict(set)
    with open(path, "r") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            hugo[row["ensid"]].add(row["hugo"])

    return hugo


def parse_cached_domains():
    """
    Method to read the domain cache into a dict.
    """
    if not os.path.isfile(DOMAIN_CACHE_PATH):
        print("Creating domain cache:", DOMAIN_CACHE_PATH, file=sys.stderr)
    else:
        with open(DOMAIN_CACHE_PATH, "r") as fh:
            for line in fh:
                protein_id, encoded = line.split("\t")
                decode = json.loads(encoded)
                DOMAIN_CACHE[protein_id] = decode
        last_mod = get_date(os.path.getmtime(DOMAIN_CACHE_PATH))
        print("Loaded domain dictionary from:", DOMAIN_CACHE_PATH, file=sys.stderr)
        print("Domain dictionary last modified:", last_mod, file=sys.stderr)


def cached_domains(func):
    """
    Decorator to check the cache for domain information first, before requesting from Ensembl
    """

    @wraps(func)
    def check_cache(*args):
        protein_id = args[1]
        if protein_id in DOMAIN_CACHE:
            return DOMAIN_CACHE[protein_id]
        result = func(*args)
        DOMAIN_CACHE[protein_id] = result
        with open(DOMAIN_CACHE_PATH, "a") as fh:
            print(protein_id, json.dumps(result), sep="\t", file=fh)

        return result

    return check_cache


# rate_limited() adapted from: https://gist.github.com/gregburek/1441055
def rate_limited(max_per_second):
    """
    Decorator to limit the number of function calls per second
    """
    min_interval = 1.0 / float(max_per_second)

    def decorate(func):
        last_time_called = [0.0]

        @wraps(func)
        def rate_limited_function(*args, **kwargs):
            elapsed = time.process_time() - last_time_called[0]
            left_to_wait = min_interval - elapsed

            if left_to_wait > 0:
                time.sleep(left_to_wait)

            ret = func(*args, **kwargs)
            last_time_called[0] = time.process_time()
            return ret

        return rate_limited_function

    return decorate


@rate_limited(0.1)
def request_ensembl_protein(protein_id):
    """
    Method request domain info for a given protein ID from the Ensembl-REST APi.
    The 'pysembl' module used in the rest of this script has no method of
    retrieving domain info, so we have to do it manually.

    Args:
        protein_id (str): Ensembl protein ID

    Returns:
        dict: a dictionary containing domain info for the given protein ID
    """
    url = "http://rest.ensembl.org/overlap/translation/{}?db_type=core".format(protein_id)
    headers = {"Content-Type": "application/json"}
    response = requests.get(url=url, headers=headers)

    # more on Ensembl response codes: https://github.com/Ensembl/ensembl-rest/wiki/HTTP-Response-Codes
    if not response.ok:
        rcode = response.status_code
        if rcode == 400:
            print("No Ensembl entry matching: {}".format(protein_id), file=sys.stderr)
            return []
        else:
            if 400 < rcode <= 499:
                print("Ensembl API refusing requests (error {})".format(rcode), file=sys.stderr)
            elif 500 <= rcode <= 599:
                print(
                    "Ensembl REST API is currently unavailable (error {}). Try again later.".format(
                        rcode
                    ),
                    file=sys.stderr,
                )
            else:
                print(
                    "Request for {} failed. Unexpected response (error {}).".format(
                        protein_id, rcode
                    ),
                    file=sys.stderr,
                )
            # response.raise_for_status()

    # 'X-RateLimit-Remaining' is the number of requests we can make before getting cut off
    # this value resets to 'X-RateLimit-Limit' after 'X-RateLimit-Reset' seconds
    try:
        headers = response.headers
        remain = int(headers["X-RateLimit-Remaining"])
        reset = int(headers["X-RateLimit-Reset"]) + 1
        if remain < 1:
            print(
                "Ensembl API refusing requests. Pausing for {}s...".format(reset), file=sys.stderr
            )
            time.sleep(reset)
    except KeyError:
        print("Could not parse API response headers", file=sys.stderr)

    decoded = response.json()
    return decoded


class EnsemblAnnotation(object):
    """
    Class for building an annotation file for MAVIS in json format.
    Args:
        species (str): species of interest
        release (int): Ensembl release to use
        output (str): path to output file
        best_file (str): path to file of "best transcripts"
        hugo_file (str): path to file with transcript aliases
    """

    def __init__(self, release, species, output, best_file=None, hugo_file=None, custom_cache=None):

        self.annotation = {}

        self.custom_cache = custom_cache
        self.cache_prefix = None
        self.gen_time = get_date()
        self.release = release
        self.species = species
        self.output = output

        self.best_file = best_file
        self.hugo_file = hugo_file

        if self.best_file:
            self.best = parse_best_file(self.best_file)
        else:
            self.best = set()

        if self.hugo_file:
            self.hugo = parse_hugo_file(self.hugo_file)
        else:
            self.hugo = defaultdict(set)

        self.data = EnsemblRelease(release, species)
        self.download_pyensembl_cache()
        self.get_domain_cache()
        self.build_json()

    def download_pyensembl_cache(self):
        """
        Method download the pyensembl cache files for this release if not already there.
        Args:
            data (EnsemblRelease): pyensembl object for the release info
            custom_cache (str): path to cirectory to cache pyensembl files
        """
        if self.custom_cache:
            os.environ["PYENSEMBL_CACHE_DIR"] = self.custom_cache
        self.data.download()
        self.data.index()
        self.cache_prefix = self.data.gtf_path.split("gtf.gz")[0]

    def get_domain_cache(self):
        global DOMAIN_CACHE_PATH
        DOMAIN_CACHE_PATH = self.cache_prefix + "domain.tsv"
        parse_cached_domains()

    def get_genes(self, eid):
        """
        Method parse gene info in the EnsemblRelease into json format.
        Args:
            eid (str): Ensembl gene ID
        Returns:
            dict: gene info formatted for json
        """
        gene = self.data.gene_by_id(eid)
        result = {
            "name": str(gene.gene_id),
            "chr": str(gene.contig),
            "start": int(gene.start),
            "end": int(gene.end),
            "strand": str(gene.strand),
            "aliases": [str(gene.gene_name)] + list(self.hugo[gene.gene_id]),
            "transcripts": [],
        }

        return result

    def get_transcripts(self, eid):
        """
        Method parse transcript info in the EnsemblRelease into json format.
               Ignore non-coding transcripts.
        Args:
            eid (str): Ensembl transcript ID
        Returns:
            dict: transcript info formatted for json
        """
        transcript = self.data.transcript_by_id(eid)
        protein_id = transcript.protein_id
        if not protein_id:
            return None

        result = {
            "name": str(transcript.transcript_id),
            "start": int(transcript.start),
            "end": int(transcript.end),
            "aliases": [str(transcript.transcript_name)],
            "is_best_transcript": str(transcript.transcript_id) in self.best,
            "protein_id": transcript.protein_id,
            "exons": [],
            "domains": [],
        }

        # start/end are absolute genomic positions, so calculate positions relative to the mRNA start
        cpos = transcript.coding_sequence_position_ranges
        if transcript.strand in ("+", "1"):
            result["cdna_coding_start"] = transcript.spliced_offset(cpos[0][0]) + 1
            result["cdna_coding_end"] = transcript.spliced_offset(cpos[-1][1]) + 1
        elif transcript.strand in ("-", "-1"):
            result["cdna_coding_start"] = transcript.spliced_offset(cpos[0][1]) + 1
            result["cdna_coding_end"] = transcript.spliced_offset(cpos[-1][0]) + 1

        return result

    def get_exons(self, eid):
        """
        Method parse exon info in the EnsemblRelease into json format.
        Args:
            eid (str): Ensembl exon ID
        Returns:
            dict: exon info formatted for json
        """
        exon = self.data.exon_by_id(eid)
        result = {"name": str(exon.exon_id), "start": int(exon.start), "end": int(exon.end)}

        return result

    @cached_domains
    def get_domains(self, eid):
        """
        Method request domain info from Ensembl and parse into json format.
        Args:
            eid (str): Ensembl protein ID
        Returns:
            list: a list of domains formatted for json
        """
        temp = {}

        protein = request_ensembl_protein(eid)
        for domain in protein:
            name = str(domain["id"])
            desc = (
                str(domain["description"]).replace('"', "").replace("'", "")
            )  # quotes causing errors when mavis loads json
            region = {"start": int(domain["start"]), "end": int(domain["end"])}
            if desc == "":
                continue
            if name in temp:
                temp[name]["regions"].append(region)
            else:
                temp[name] = {"name": name, "desc": desc, "regions": [region]}

        domain_list = list(temp.values())
        return domain_list

    def build_json(self):
        """
        Method compile a json object for MAVIS of all protein coding genes and
               associated info for the indicated species.
        Returns:
            dict: a json-formatted set of annotations for use with MAVIS
        """
        count = {"gene": 0, "transcript": 0, "non_coding": 0}

        self.annotation["script"] = SCRIPT
        self.annotation["script_version"] = VERSION
        self.annotation["hugo_mapping_file"] = self.hugo_file
        self.annotation["best_transcript_file"] = self.best_file
        self.annotation["ensembl_version"] = self.release
        self.annotation["generation_time"] = self.gen_time
        self.annotation["genes"] = []

        gene_ids = self.data.gene_ids()

        for index, gid in enumerate(gene_ids):
            print("{}/{} genes".format(index, len(gene_ids)))
            gened = self.get_genes(gid)
            count["gene"] += 1
            for tid in self.data.transcript_ids_of_gene_id(gid):
                transd = self.get_transcripts(tid)
                count["transcript"] += 1
                if transd:
                    for eid in self.data.exon_ids_of_transcript_id(tid):
                        exond = self.get_exons(eid)
                        transd["exons"].append(exond)
                    domains = self.get_domains(transd["protein_id"])
                    transd["domains"] = domains
                    gened["transcripts"].append(transd)
                else:
                    count["non_coding"] += 1
            if gened["transcripts"] != []:
                self.annotation["genes"].append(gened)

        print(
            "{gene:,} genes, {transcript:,} transcripts ({non_coding:,} non-coding transcripts, ignored)".format(
                **count
            )
        )
        return self.annotation

    def dump_json(self):
        """
        Method dump the annotations in json-format to the specified output file.
        """
        with open(self.output, "w") as fh:
            json.dump(self.annotation, fh)

    def delete_cache(self):
        """
        Method delete both the pyensembl and domain cache files.
        """
        for cache_file in glob(self.cache_prefix + "*"):
            print("Removing cache file", cache_file)
            os.remove(cache_file)


def main():
    """
    Main loop.
    """
    parser = argparse.ArgumentParser(
        description="Generate an annotation file in json format for use with MAVIS by pulling info from Ensembl.",
        add_help=False,
    )

    req_parser = parser.add_argument_group("required arguments")
    req_parser.add_argument(
        "-s", "--species", required=True, help="species of interest (human, mouse, rat, etc.)"
    )
    req_parser.add_argument(
        "-r", "--release", required=True, type=int, help="Ensembl release to use"
    )
    req_parser.add_argument("-o", "--output", required=True, help="path to output json file")

    opt_parser = parser.add_argument_group("optional arguments")
    opt_parser.add_argument(
        "-b",
        "--best-transcript-file",
        help="tab-seperated list of genes and corresponding transcripts",
    )
    opt_parser.add_argument("-m", "--hugo-mapping-file", help="one-to-one hugo gene name mapping")
    opt_parser.add_argument(
        "-c", "--custom-cache", help="use a non-default path to cache ensembl data"
    )
    opt_parser.add_argument(
        "-d",
        "--delete-cache",
        action="store_true",
        help="delete cache files when the script finishes successfully",
    )
    opt_parser.add_argument("-h", "--help", action="help", help="show this help message and exit")
    args = parser.parse_args()

    print("If this is the first time running, may take a while!")
    print("Start: {}".format(get_date()))
    annotation = EnsemblAnnotation(
        release=args.release,
        species=args.species,
        output=args.output,
        best_file=args.best_transcript_file,
        hugo_file=args.hugo_mapping_file,
        custom_cache=args.custom_cache,
    )
    annotation.dump_json()
    if args.delete_cache:
        annotation.delete_cache()
    print("End: {}".format(get_date()))


if __name__ == "__main__":
    main()
