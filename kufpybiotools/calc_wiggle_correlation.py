#!/usr/bin/env python
"""
Copyright (c) 2013, Konrad Foerstner <konrad@foerstner.org>

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
__description__ = ("Calculate the Pearson correlation coefficient "
                   "for the coverages of two wiggle files.")
__author__ = "Konrad Foerstner <konrad@foerstner.org>"
__copyright__ = "2013 by Konrad Foerstner <konrad@foerstner.org>"
__license__ = "ISC license"
__email__ = "konrad@foerstner.org"
__version__ = ""

import argparse

from wiggle import WiggleParser
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument("wiggle_file_1", type=argparse.FileType("r"))
    parser.add_argument("wiggle_file_2", type=argparse.FileType("r"))
    args = parser.parse_args()
    wiggel_correlator = WiggleCorrelator()
    wiggel_correlator.correlate(args.wiggle_file_1, args.wiggle_file_2)

class WiggleCorrelator(object):

    def __init__(self):
        self._wiggle_parser = WiggleParser()

    def correlate(self, wiggle_file_1, wiggle_file_2):
        self._chrom_and_pos_value_pairs = {}
        print("Replicon: Pearson correlation coefficient (p-value)")
        for entry_1, entry_2 in zip(
                self._wiggle_parser.entries(wiggle_file_1),
                self._wiggle_parser.entries(wiggle_file_2)):
                assert(entry_1.replicon == entry_2.replicon)
                pos_value_pairs_1 = dict(entry_1.pos_value_pairs)
                pos_value_pairs_2 = dict(entry_2.pos_value_pairs)
                if len(pos_value_pairs_1) == 0 or len(pos_value_pairs_2) == 0:
                    print("%s: At least one replicon has no coverage for "
                          "this libs." % (entry_2.replicon))
                    continue
                non_redu_pos = set(
                    pos_value_pairs_1.keys() + pos_value_pairs_2.keys())
                values_1 = np.array(
                    [pos_value_pairs_1.get(pos, 0.0) for pos in non_redu_pos])
                values_2 = np.array(
                    [pos_value_pairs_2.get(pos, 0.0) for pos in non_redu_pos])
                pearson, pvalue = stats.pearsonr(values_1, values_2)
                print("%s: %s (%s)" % (entry_1.replicon, pearson, pvalue))

if __name__ == "__main__":
   main()


