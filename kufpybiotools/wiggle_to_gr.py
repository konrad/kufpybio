#!/usr/bin/env python

from kufpybio.wiggle import WiggleParser

"""
Copyright (c) 2012, Konrad Foerstner <konrad@foerstner.org>

Permission to use, copy, modify, and/or distribute this software for
any purpose with or without fee is hereby granted, provided that the
above copyright notice and this permission notice appear in all
copies.

THE SOFTWARE IS PROVIDED 'AS IS' AND THE AUTHOR DISCLAIMS ALL
WARRANTIES WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE
AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL
DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER
TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
PERFORMANCE OF THIS SOFTWARE.
         
"""
__description__ = ""
__author__ = "Konrad Foerstner <konrad@foerstner.org>"
__copyright__ = "2012 by Konrad Foerstner <konrad@foerstner.org>"
__license__ = "ISC license"
__email__ = "konrad@foerstner.org"
__version__ = ""

import argparse

def main():
    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument("input_file")
    parser.add_argument("--output_prefix", default=None)
    args = parser.parse_args()
    wiggle_to_gr_converter = WiggleToGrConverter()
    wiggle_to_gr_converter.convert(
        args.input_file, args.output_prefix)

class WiggleToGrConverter(object):

    def convert(self, input_file, output_prefix):
        wiggle_parser = WiggleParser()
        for entry in wiggle_parser.entries(open(input_file)):
            output_fh = open(self._output_file(
                    input_file, output_prefix, entry), "w")
            for pos_value_pair in entry.pos_value_pairs:
                output_fh.write("\t".join([
                           str(val) for val in pos_value_pair]) + "\n")

    def _output_file(self, input_file, output_prefix, entry):
        if output_prefix:
            return("%s_%s_in_%s.gr" % (
                    output_prefix, entry.track_name, entry.chrom_name))
        else:
            return("%s_%s_in_%s.gr" % (
                    input_file, entry.track_name, entry.chrom_name))

if __name__ == "__main__": 
   main()
