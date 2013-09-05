#!/usr/bin/env python

__description__ = ""
__author__ = "Konrad Foerstner <konrad@foerstner.org>"
__copyright__ = "2013 by Konrad Foerstner <konrad@foerstner.org>"
__license__ = "ISC license"
__email__ = "konrad@foerstner.org"
__version__ = ""

import argparse
import csv
import sys
sys.path.append(".")
from kufpybio.gff3 import Gff3Parser, Gff3Entry
from kufpybio.gene import Gene
from kufpybio.igrfinder import IGRFinder

parser = argparse.ArgumentParser(description=__description__)
parser.add_argument("gff_file", type=argparse.FileType("r"))
parser.add_argument("output_file", type=argparse.FileType("w"))
args = parser.parse_args()
# Build gene list
gene_list = []
gff_parser = Gff3Parser()
region_entry = None
for entry in gff_parser.entries(args.gff_file):
    if entry.feature == "region":
        region_entry = entry
        continue
    gene_list.append(Gene(
        entry.seq_id, "", "", entry.start, entry.end,
        entry.strand))
# Find IGRs and generate GFF file
igr_finder = IGRFinder()
args.output_file.write("##gff-version 3\n")
for start, end in igr_finder.find_igrs(gene_list, region_entry.end):
    for strand in ["+", "-"]:
        gff3_entry = Gff3Entry({
            "seq_id" : region_entry.seq_id,
            "source" : "IGR",
            "feature" : "IGR",
            "start" : start,
            "end" : end,
            "score" : ".",
            "strand" : strand,
            "phase" : ".",
            "attributes" : "ID=IGR_%s_%s_to_%s" % (
                    region_entry.seq_id, start, end)})
        args.output_file.write(str(gff3_entry) + "\n")
