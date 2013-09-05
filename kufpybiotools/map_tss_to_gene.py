#!/usr/bin/env python

__description__ = "Allocates TSS to given gene annotations."
__author__ = "Konrad Foerstner <konrad@foerstner.org>"
__copyright__ = "2013 by Konrad Foerstner <konrad@foerstner.org>"
__license__ = "ISC license"
__email__ = "konrad@foerstner.org"
__version__ = ""

import argparse
import csv
import sys
sys.path.append(".")
from kufpybio.gff3 import Gff3Parser
from kufpybio.tss import TSS
from kufpybio.gene import Gene
import kufpybio.tssgenemapper as tssgenemapper

def main():
    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument("tss_list_file", type=argparse.FileType("r"))
    parser.add_argument("gff_file", type=argparse.FileType("r"))
    parser.add_argument("output_file", type=argparse.FileType("w"))
    parser.add_argument("--min_dist_to_gene_end", default=100)
    parser.add_argument("--orphan_distance", default=300)
    args = parser.parse_args()
    mapper = Mapper(
        args.tss_list_file, args.gff_file, args.output_file, 
        args.min_dist_to_gene_end, args.orphan_distance)
    mapper.create_tss_list()
    mapper.create_gene_list()
    mapper.map_tss()
    mapper.write_output()

class Mapper(object):

    def __init__(self, tss_list_fh, gff_fh, output_fh, min_dist_to_gene_end, 
                 orphan_distance):
        self.tss_list_fh = tss_list_fh
        self.gff_fh = gff_fh
        self.output_fh = output_fh
        self.min_dist_to_gene_end = min_dist_to_gene_end
        self.orphan_distance = orphan_distance
    
    def create_tss_list(self):
        self.tss_list = []
        for row in csv.reader(self.tss_list_fh, delimiter="\t"):
            seq_id = None
            if len(row) == 3:
                seq_id == row[2]
            try:
                self.tss_list.append(TSS(seq_id, row[0], row[1]))
            except:
                sys.stderr.write("Skipping TSS table line: \"%s\"\n" % (
                        "\t".join(row)))

    def create_gene_list(self):
        self.gene_list = []
        gff_parser = Gff3Parser()
        for entry in gff_parser.entries(self.gff_fh):
            if entry.feature != "gene":
                continue
            self.gene_list.append(Gene(
                    entry.seq_id, entry.attributes["locus_tag"], 
                    entry.attributes["Name"], entry.start, entry.end, 
                    entry.strand))
            
    def map_tss(self):
        tss_gene_mapper = tssgenemapper.TSSGeneMapper()
        self.tss_and_genes = tss_gene_mapper.map_tss(
            self.tss_list, self.gene_list)

    def write_output(self):
        self.tss_gene_formatter = tssgenemapper.TSSGeneFormatter()
        self._write_header()
        for tss_pos, tss in sorted([(tss.pos, tss) 
                                    for tss in self.tss_and_genes.keys()]):
            if self.tss_and_genes[tss] == tssgenemapper.orphan_str:
                self._write_orphan(tss)
            else:
                for gene_start, gene in sorted(
                    [(gene.start, gene) 
                     for gene in self.tss_and_genes[tss].keys()]):
                    self._write_tss_with_gene(tss, gene)
                
    def _write_tss_with_gene(self, tss, gene):
        tss_features = self.tss_and_genes[tss][gene]
        bin_features = self.tss_gene_formatter.tss_features_binary_format(
            tss_features)
        feature_string = self.tss_gene_formatter.tss_features_string(
            tss_features)
        # For some constellation there are no binary foramt of the
        # features as they will and are not being considered.
        if not bin_features:
            return
        utr_len = "-"
        if (tss_features["location"] == tssgenemapper.loc_5_prime_str):
            utr_len = str(self._utr_length(tss, gene))
        self.output_fh.write(
            "\t".join([str(tss.pos), tss.strand] + [
                    gene.gene_id, gene.name, str(gene.start), str(gene.end), 
                    gene.strand, str(gene.end-gene.start+1)] + 
                      bin_features + [feature_string] + [utr_len]) + "\n")

    def _write_orphan(self, tss):
        bin_features = self.tss_gene_formatter.tss_features_binary_format(
            tssgenemapper.orphan_str)
        self.output_fh.write(
            "\t".join([str(tss.pos), tss.strand] + ["-"] * 6 + 
                       bin_features + ["orphan", "-"]) + "\n")

    def _write_header(self):
        self.output_fh.write("\t".join(
                ["TSS pos", "TSS strand", "Gene id", "Gene name", 
                 "Gene start", "Gene end", "Gene strand", "Gene length"] +
                self.tss_gene_formatter.binary_format_header + 
                ["Status", "UTR length"]) + "\n")

    def _utr_length(self, tss, gene):
        if tss.strand == "+":
            return(gene.start - tss.pos)
        else:
            return(tss.pos - gene.end)
                
if __name__ == "__main__": 
   main()
