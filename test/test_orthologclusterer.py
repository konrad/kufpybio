import unittest
import sys
from kufpybio.orthologclusterer import OrthologClusterer, Ortholog

class TestOrthologClusterer(unittest.TestCase):
    
    def setUp(self):
        self._ortholog_clusterer = OrthologClusterer()

    def test_process_pairs_to_clusters_3_gene_in_cluster(self):
        self._ortholog_clusterer.process_pairs_to_clusters(
            [(("gene_A1", "org_A"), ("gene_B1", "org_B")),
             (("gene_A1", "org_A"), ("gene_C1", "org_C"))])
        orthologs_and_clusters = self._ortholog_clusterer._orthologs_and_clusters
        # There must be three orthologs
        self.assertEqual(sorted(
                orthologs_and_clusters.keys()), 
                         ['gene_A1_org_A', 'gene_B1_org_B', 'gene_C1_org_C'])
        # All 3 othologos must refer to the same cluster
        self.assertEqual(orthologs_and_clusters['gene_A1_org_A'],
                         orthologs_and_clusters['gene_B1_org_B'])
        self.assertEqual(orthologs_and_clusters['gene_A1_org_A'],
                         orthologs_and_clusters['gene_C1_org_C'])
        # The set must contain 3 orthologs
        self.assertEqual(len(orthologs_and_clusters['gene_A1_org_A']), 3)

    def test_process_pairs_to_clusters_2_gene_in_cluster_1_lonely_gene(self):
        self._ortholog_clusterer.process_pairs_to_clusters(
            [(("gene_A1", "org_A"), ("gene_B1", "org_B")),
             (("gene_C1", "org_C"), None)])
        orthologs_and_clusters = self._ortholog_clusterer._orthologs_and_clusters
        # There must be three orthologs
        self.assertEqual(
            sorted(orthologs_and_clusters.keys()), 
            ['gene_A1_org_A', 'gene_B1_org_B', 'gene_C1_org_C'])        
        # A and B must refer to the same cluster
        self.assertEqual(orthologs_and_clusters['gene_A1_org_A'],
                         orthologs_and_clusters['gene_B1_org_B'])
        # C must refer to a cluster of the size 1
        self.assertEqual(len(orthologs_and_clusters[
                    'gene_C1_org_C']), 1)
        
    def test_clusters(self):
        self._ortholog_clusterer.process_pairs_to_clusters(
            [(("gene_A1", "org_A"), ("gene_B1", "org_B")),
             (("gene_A1", "org_A"), ("gene_C1", "org_C")),
             (("gene_D1", "org_D"), ("gene_E1", "org_E"))])
        clusters = list(self._ortholog_clusterer.clusters())
        assert len(clusters) == 2
        
class TestOrtholog(unittest.TestCase):

    def test_init(self):
        ortholog = Ortholog("gene_A1", "org_A")
        self.assertEqual(ortholog.accession, "gene_A1")
        self.assertEqual(ortholog.source, "org_A")
        self.assertEqual(ortholog.key, "gene_A1_org_A")
        self.assertEqual(str(ortholog), "gene_A1_org_A")

if __name__ == "__main__":
    unittest.main()
