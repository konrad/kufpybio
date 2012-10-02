import unittest
import sys
import kufpybio.genelistmerger as genelistmerger
from kufpybio.gene import Gene


class TestTSSGeneMapper_1(unittest.TestCase):
    """Use default setting. Any overlap sufficient"""

    def setUp(self):
        self.gene_list_merger = genelistmerger.GeneListMerger()

    def test_have_sufficient_overlap_1(self):
        gene_1 = Gene("genomeX", "g", "g" , 1, 50, "+")
        gene_2 = Gene("genomeX", "g", "g" , 200, 450, "+")
        self.assertFalse(self.gene_list_merger._have_sufficient_overlap(
            gene_1, gene_2))

    def test_have_sufficient_overlap_2(self):
        gene_1 = Gene("genomeX", "g", "g" , 1, 100, "+")
        gene_2 = Gene("genomeX", "g", "g" , 100, 150, "+")
        self.assertTrue(self.gene_list_merger._have_sufficient_overlap(
            gene_1, gene_2))

    def test_have_sufficient_overlap_3(self):
        gene_1 = Gene("genomeX", "g", "g" , 1, 100, "+")
        gene_2 = Gene("genomeX", "g", "g" , 50, 60, "+")
        self.assertTrue(self.gene_list_merger._have_sufficient_overlap(
            gene_1, gene_2))

    def test_overlap_percentage_1(self):
        gene = Gene("genomeX", "g", "g" , 1, 100, "+")
        self.assertEqual(self.gene_list_merger._overlap_percentage(gene, 50), 50)

class TestTSSGeneMapper_2(unittest.TestCase):
    """Set and test min overlap."""

    def setUp(self):
        self.gene_list_merger = genelistmerger.GeneListMerger(50)

    def test_have_sufficient_overlap_1(self):
        gene_1 = Gene("genomeX", "g", "g" , 1, 100, "+")
        gene_2 = Gene("genomeX", "g", "g" , 50, 150, "+")
        self.assertTrue(self.gene_list_merger._have_sufficient_overlap(
            gene_1, gene_2))

    def test_have_sufficient_overlap_2(self):
        gene_1 = Gene("genomeX", "g", "g" , 1, 100, "+")
        gene_2 = Gene("genomeX", "g", "g" , 51, 150, "+")
        self.assertFalse(self.gene_list_merger._have_sufficient_overlap(
            gene_1, gene_2))

    def test_have_sufficient_overlap_3(self):
        gene_1 = Gene("genomeX", "g", "g" , 1, 100, "+")
        gene_2 = Gene("genomeX", "g", "g" , 1, 10, "+")
        self.assertFalse(self.gene_list_merger._have_sufficient_overlap(
            gene_1, gene_2))

    def test_have_sufficient_overlap_4(self):
        gene_1 = Gene("genomeX", "g", "g" , 50, 150, "+")
        gene_2 = Gene("genomeX", "g", "g" , 1, 100, "+")
        self.assertTrue(self.gene_list_merger._have_sufficient_overlap(
            gene_1, gene_2))

    def test_have_sufficient_overlap_5(self):
        gene_1 = Gene("genomeX", "g", "g" , 51, 150, "+")
        gene_2 = Gene("genomeX", "g", "g" , 1, 100, "+")
        self.assertFalse(self.gene_list_merger._have_sufficient_overlap(
            gene_1, gene_2))

    def test_have_sufficient_overlap_6(self):
        gene_1 = Gene("genomeX", "g", "g" , 1, 10, "+")        
        gene_2 = Gene("genomeX", "g", "g" , 1, 100, "+")
        self.assertFalse(self.gene_list_merger._have_sufficient_overlap(
            gene_1, gene_2))

if __name__ == "__main__":
    unittest.main()
