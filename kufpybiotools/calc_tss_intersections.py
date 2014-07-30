#!/usr/bin/env python
"""Read a TSSpredator master table and caclulates the sizes of set
   intersections

"""

__description__ = ""
__author__ = "Konrad Foerstner <konrad@foerstner.org>"
__copyright__ = "2014 by Konrad Foerstner <konrad@foerstner.org>"
__license__ = "ISC license"
__email__ = "konrad@foerstner.org"
__version__ = ""

from collections import defaultdict
import itertools
import argparse
from tsspredator import TSSPredatorReader

parser = argparse.ArgumentParser(description=__description__)
parser.add_argument("ttspredator_matertable", nargs="+")
args = parser.parse_args()

libs_and_tsss = defaultdict(set)

for master_table_file in args.ttspredator_matertable:
    tss_predator_reader = TSSPredatorReader()
    for entry in tss_predator_reader.entries(open(master_table_file)):
        # Just require detected not enriched!
        if not (entry.is_detected):
            continue
        libs_and_tsss[entry.genome].add("%s-%s-%s" % (
            master_table_file, entry.pos, entry.strand))

print("Total numbers of TSS:")
for lib, tsss in libs_and_tsss.items():
    print("- %s: %s" % (lib, len(tsss)))

print("Size of intersections:")
for combo_size in reversed(range(2, len(libs_and_tsss.keys())+1)):
    for combination in itertools.combinations(libs_and_tsss.keys(), combo_size):
        curr_intersection_set = libs_and_tsss[combination[0]]
        for lib in combination[1:]:
            curr_intersection_set = curr_intersection_set.intersection(libs_and_tsss[lib])
        print("- %s: %s" % (", ".join(combination), len(curr_intersection_set)))

for lib, tsss in libs_and_tsss.items():
    with open("%s.csv" % (lib), "w") as fh:
        fh.write("\n".join(tsss) + "\n")
