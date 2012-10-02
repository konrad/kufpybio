import unittest
import sys
from kufpybio.gene import Gene

class TestGene(unittest.TestCase):

    def test_init_without_extra(self):
        gene = Gene("chr1", "gene_2342", "hacY", "3", "83", "+")
        self.assertEqual(gene.seq_id, "chr1")
        self.assertEqual(gene.gene_id, "gene_2342")
        self.assertEqual(gene.name, "hacY")
        self.assertEqual(gene.start, 3)
        self.assertEqual(gene.end, 83)
        self.assertEqual(gene.strand, "+")
        self.assertFalse(hasattr(gene, "extra"))

    def test_init_with_extra(self):
        gene = Gene("chr2", "gene_0005", "hacZ", "15", "30", "-", extra="mope")
        self.assertEqual(gene.seq_id, "chr2")
        self.assertEqual(gene.gene_id, "gene_0005")
        self.assertEqual(gene.name, "hacZ")
        self.assertEqual(gene.start, 15)
        self.assertEqual(gene.end, 30)
        self.assertEqual(gene.strand, "-")
        self.assertEqual(gene.extra, "mope")

    def test_init_start_end_sorting(self):
        """Test that the start and end position are ordered"""
        gene = Gene("chr1", "gene_2342", "hacY", "1000", "5", "+")
        self.assertEqual(gene.start, 5)
        self.assertEqual(gene.end, 1000)

    def test_len(self):
        # 1234567890
        # ==========
        gene1 = Gene("chr1", "a", "b", "1", "10", "+")
        self.assertEqual(gene1.len(), 10)
        # 01234567890
        # ===========
        gene2 = Gene("chr1", "a", "b", "10", "20", "+")
        self.assertEqual(gene2.len(), 11)

if __name__ == "__main__":
    unittest.main()
