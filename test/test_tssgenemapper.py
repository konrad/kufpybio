import unittest
import sys
from io import StringIO
import kufpybio.tssgenemapper as tssgenemapper
from kufpybio.gene import Gene
from kufpybio.tss import TSS

class TestTSSGeneMapper(unittest.TestCase):

    def setUp(self):
        self.tss_gene_mapper = tssgenemapper.TSSGeneMapper()
        self.tss_gene_mapper.tss_and_hit_genes = {}
        self.tss_gene_mapper.genes_and_tss = {}

    def test_5_prime_dist_1(self):
        tss = TSS("genomeX", 5, "+")
        gene = Gene("genomeX", "g", "g" , 10, 100, "+")
        self.assertEqual(self.tss_gene_mapper._5_prime_dist(tss, gene), 5)

    def test_5_prime_dist_2(self):
        tss = TSS("genomeX", 105, "+")
        gene = Gene("genomeX", "g", "g" , 10, 100, "+")
        self.assertEqual(self.tss_gene_mapper._5_prime_dist(tss, gene), -95)

    def test_5_prime_dist_3(self):
        tss = TSS("genomeX", 2, "-")
        gene = Gene("genomeX", "g", "g" , 10, 100, "-")
        self.assertEqual(self.tss_gene_mapper._5_prime_dist(tss, gene), -98)

    def test_5_prime_dist_4(self):
        tss = TSS("genomeX", 120, "-")
        gene = Gene("genomeX", "g", "g" , 10, 100, "-")
        self.assertEqual(self.tss_gene_mapper._5_prime_dist(tss, gene), 20)

    def test_has_5_prime_association_1(self):
        """True - In 5' region"""
        tss = TSS("genomeX", 3, "+")
        gene = Gene("genomeX", "g", "g" , 10, 100, "+")
        self.assertEqual(
                self.tss_gene_mapper._has_5_prime_association(tss, gene), True)

    def test_has_5_prime_association_2(self):
        """None - in 3' region"""
        tss = TSS("genomeX", 200, "+")
        gene = Gene("genomeX", "g", "g" , 10, 100, "+")
        self.assertEqual(
                self.tss_gene_mapper._has_5_prime_association(tss, gene), None)

    def test_has_5_prime_association_3(self):
        """True - Leaderless TSS"""
        tss = TSS("genomeX", 10, "+")
        gene = Gene("genomeX", "g", "g" , 10, 100, "+")
        self.assertEqual(
                self.tss_gene_mapper._has_5_prime_association(tss, gene), True)

    def test_has_5_prime_association_4(self):
        """None - in 5' but too far away """
        tss = TSS("genomeX", 10, "+")
        gene = Gene("genomeX", "g", "g" , 1000, 1100, "+")
        self.assertEqual(
                self.tss_gene_mapper._has_5_prime_association(tss, gene), None)

    def test_has_5_prime_association_5(self):
        """None - antisense"""
        tss = TSS("genomeX", 110, "+")
        gene = Gene("genomeX", "g", "g" , 10, 100, "-")
        self.assertEqual(
                self.tss_gene_mapper._has_5_prime_association(tss, gene), None)

    def test_has_5_prime_association_6(self):
        """None - antisense"""
        tss = TSS("genomeX", 7, "+")
        gene = Gene("genomeX", "g", "g" , 10, 100, "-")
        self.assertEqual(
                self.tss_gene_mapper._has_5_prime_association(tss, gene), None)

    def test_has_5_prime_association_7(self):
        """None - antisense"""
        tss = TSS("genomeX", 100, "+")
        gene = Gene("genomeX", "g", "g" , 10, 100, "-")
        self.assertEqual(
                self.tss_gene_mapper._has_5_prime_association(tss, gene), None)

    def test_has_5_prime_association_8(self):
        """None - antisense"""
        tss = TSS("genomeX", 600, "+")
        gene = Gene("genomeX", "g", "g" , 10, 100, "-")
        self.assertEqual(
                self.tss_gene_mapper._has_5_prime_association(tss, gene), None)
        
    def test_has_5_prime_association_9(self):
        """True - In 5' region"""
        tss = TSS("genomeX", 105, "-")
        gene = Gene("genomeX", "g", "g" , 10, 100, "-")
        self.assertEqual(
                self.tss_gene_mapper._has_5_prime_association(tss, gene), True)

    def test_has_5_prime_association_10(self):
        """None - in 3' region"""
        tss = TSS("genomeX", 5, "-")
        gene = Gene("genomeX", "g", "g" , 10, 100, "-")
        self.assertEqual(
                self.tss_gene_mapper._has_5_prime_association(tss, gene), None)

    def test_has_5_prime_association_11(self):
        """True - Leaderless TSS"""
        tss = TSS("genomeX", 100, "-")
        gene = Gene("genomeX", "g", "g" , 10, 100, "-")
        self.assertEqual(
                self.tss_gene_mapper._has_5_prime_association(tss, gene), True)

    def test_has_5_prime_association_12(self):
        """None - in 5' but too far away """
        tss = TSS("genomeX", 1500, "+")
        gene = Gene("genomeX", "g", "g" , 100, 100, "+")
        self.assertEqual(
                self.tss_gene_mapper._has_5_prime_association(tss, gene), None)

    def test_has_antisense_association_1(self):
        """True - antisense and inside the range"""
        tss = TSS("genomeX", 50, "+")
        gene = Gene("genomeX", "g", "g" , 10, 100, "-")
        self.assertEqual(
                self.tss_gene_mapper._has_antisense_association(tss, gene), True)

    def test_has_antisense_association_2(self):
        """True - antisense and inside the range"""        
        tss = TSS("genomeX", 900, "+")
        gene = Gene("genomeX", "g", "g" , 1000, 1100, "-")
        self.assertEqual(
                self.tss_gene_mapper._has_antisense_association(tss, gene), True)

    def test_has_antisense_association_3(self):
        """True - antisense"""
        tss = TSS("genomeX", 200, "+")
        gene = Gene("genomeX", "g", "g" , 10, 100, "-")
        self.assertEqual(
                self.tss_gene_mapper._has_antisense_association(tss, gene), True)

    def test_has_antisense_association_4(self):
        """None - out of range"""
        tss = TSS("genomeX", 899, "+")
        gene = Gene("genomeX", "g", "g" , 1000, 1100, "-")
        self.assertEqual(
                self.tss_gene_mapper._has_antisense_association(tss, gene), None)

    def test_has_antisense_association_5(self):
        """None - out of range"""
        tss = TSS("genomeX", 201, "+")
        gene = Gene("genomeX", "g", "g" , 10, 100, "-")
        self.assertEqual(
                self.tss_gene_mapper._has_antisense_association(tss, gene), None)

    def test_has_antisense_association_6(self):
        """True - antisense and inside the range"""
        tss = TSS("genomeX", 50, "-")
        gene = Gene("genomeX", "g", "g" , 10, 100, "+")
        self.assertEqual(
                self.tss_gene_mapper._has_antisense_association(tss, gene), True)

    def test_has_antisense_association_7(self):
        """True - antisense and inside the range"""        
        tss = TSS("genomeX", 900, "-")
        gene = Gene("genomeX", "g", "g" , 1000, 1100, "+")
        self.assertEqual(
                self.tss_gene_mapper._has_antisense_association(tss, gene), True)

    def test_has_antisense_association_8(self):
        """True - antisense"""
        tss = TSS("genomeX", 200, "-")
        gene = Gene("genomeX", "g", "g" , 10, 100, "+")
        self.assertEqual(
                self.tss_gene_mapper._has_antisense_association(tss, gene), True)

    def test_has_antisense_association_9(self):
        """None - out of range"""
        tss = TSS("genomeX", 899, "-")
        gene = Gene("genomeX", "g", "g" , 1000, 1100, "+")
        self.assertEqual(
                self.tss_gene_mapper._has_antisense_association(tss, gene), None)

    def test_has_antisense_association_10(self):
        """None - out of range"""
        tss = TSS("genomeX", 201, "-")
        gene = Gene("genomeX", "g", "g" , 10, 100, "+")
        self.assertEqual(
                self.tss_gene_mapper._has_antisense_association(tss, gene), None)

    def test_has_internal_association_1(self):
        """True - internal"""
        tss = TSS("genomeX", 20, "+")
        gene = Gene("genomeX", "g", "g" , 10, 100, "+")
        self.assertEqual(
                self.tss_gene_mapper._has_internal_association(tss, gene), True)

    def test_has_internal_association_2(self):
        """None - 5' region"""
        tss = TSS("genomeX", 5, "+")
        gene = Gene("genomeX", "g", "g" , 10, 100, "+")
        self.assertEqual(
                self.tss_gene_mapper._has_internal_association(tss, gene), None)

    def test_has_internal_association_3(self):
        """None - 3' region """
        tss = TSS("genomeX", 150, "+")
        gene = Gene("genomeX", "g", "g" , 10, 100, "+")
        self.assertEqual(
                self.tss_gene_mapper._has_internal_association(tss, gene), None)

    def test_has_internal_association_4(self):
        """None - on the first base => leaderless 5' """
        tss = TSS("genomeX", 10, "+")
        gene = Gene("genomeX", "g", "g" , 10, 100, "+")
        self.assertEqual(
                self.tss_gene_mapper._has_internal_association(tss, gene), None)

    def test_has_internal_association_5(self):
        """True - on the last base """
        tss = TSS("genomeX", 100, "+")
        gene = Gene("genomeX", "g", "g" , 10, 100, "+")
        self.assertEqual(
                self.tss_gene_mapper._has_internal_association(tss, gene), True)

    def test_has_internal_association_6(self):
        """None - antisense"""
        tss = TSS("genomeX", 20, "+")
        gene = Gene("genomeX", "g", "g" , 10, 100, "-")
        self.assertEqual(
                self.tss_gene_mapper._has_internal_association(tss, gene), None)

    def test_has_internal_association_7(self):
        """True - on the last base """
        tss = TSS("genomeX", 10, "-")
        gene = Gene("genomeX", "g", "g" , 10, 100, "-")
        self.assertEqual(
                self.tss_gene_mapper._has_internal_association(tss, gene), True)

    def test_has_internal_association_8(self):
        """None - on the first base => 5' leaderless """
        tss = TSS("genomeX", 100, "-")
        gene = Gene("genomeX", "g", "g" , 10, 100, "-")
        self.assertEqual(
                self.tss_gene_mapper._has_internal_association(tss, gene), None)

    def test_set_type_of_5_prime_tss(self):
        self.tss_gene_mapper.genes_and_5_prime_tss = {}
        self.tss_gene_mapper.tss_and_hit_genes = {}
        tss_primary = TSS("genomeX", 40, "+")
        tss_secondary = TSS("genomeX", 30, "+")
        gene = Gene("genomeX", "g", "g" , 50, 100, "+")
        self.tss_gene_mapper.genes_and_5_prime_tss[gene] = [
            [10, tss_primary], [20, tss_secondary]]
        self.tss_gene_mapper.tss_and_hit_genes[tss_primary] = {}
        self.tss_gene_mapper.tss_and_hit_genes[tss_secondary] = {}
        self.tss_gene_mapper.tss_and_hit_genes[tss_primary][gene] = {
            "location" : tssgenemapper.loc_5_prime_str,  
            "orientation" : tssgenemapper.sense_str, 
            "distance" : 10, "tss_type" : None}
        self.tss_gene_mapper.tss_and_hit_genes[tss_secondary][gene] = {
            "location" : tssgenemapper.loc_5_prime_str,  
            "orientation" : tssgenemapper.sense_str, 
            "distance" : 20, "tss_type" : None}
        self.tss_gene_mapper._set_type_of_5_prime_tss()
        self.assertEqual(
            self.tss_gene_mapper.tss_and_hit_genes[tss_primary][
                gene]["tss_type"], tssgenemapper.primary_str)
        self.assertEqual(
            self.tss_gene_mapper.tss_and_hit_genes[tss_secondary][
                gene]["tss_type"], tssgenemapper.secondary_str)

class TestTSSGeneFormatter(unittest.TestCase):

    def setUp(self):
        self.tss_gene_formatter = tssgenemapper.TSSGeneFormatter()
        
    def test_tss_features_binary_format_5_prime_1(self):
        self.assertListEqual(
            self.tss_gene_formatter.tss_features_binary_format({
                    "location" : tssgenemapper.loc_5_prime_str, 
                    "tss_type" : tssgenemapper.primary_str}),
            ["1", "0", "0", "0"])

    def test_tss_features_binary_format_5_prime_2(self):
        self.assertListEqual(
            self.tss_gene_formatter.tss_features_binary_format({
                    "location" : tssgenemapper.loc_5_prime_str, 
                    "tss_type" : tssgenemapper.secondary_str}),
            ["0", "1", "0", "0"])

    def test_tss_features_binary_format_internal(self):
        self.assertEqual(
            self.tss_gene_formatter.tss_features_binary_format({
                    "location" : tssgenemapper.loc_internal_str}),
            ["0", "0", "1", "0"])

    def test_tss_features_binary_format_antisense(self):
        self.assertEqual(
            self.tss_gene_formatter.tss_features_binary_format({
                    "location" : tssgenemapper.loc_antisense_str}),
            ["0", "0", "0", "1"])

if __name__ == "__main__":
    unittest.main()
