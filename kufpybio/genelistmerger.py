import kufpybio.helpers as helpers
from kufpybio.gene import Gene
import sys

class GeneListMerger(object):

    def __init__(self, min_overlap_percentage=None):
        self.result_genes = []
        self.min_overlap_percentage = min_overlap_percentage
        self._seen_genes = {}

    def add_genes(self, genes_to_add):
        for query_gene in genes_to_add:
            if self._already_seen(query_gene):
                continue
            overlapping = False
            for result_gene in self.result_genes:
                # Skip genes that are on other genetic elements or the
                # opposite strand
                if (query_gene.strand != result_gene.strand or 
                    query_gene.seq_id != result_gene.seq_id):
                    continue
                # Calculate the overlap - if sufficient merge genes
                if self._have_sufficient_overlap(query_gene, result_gene):
                    self.result_genes.remove(result_gene)
                    merged_gene = self._merge_genes(query_gene, result_gene)
                    self.result_genes.append(merged_gene)
                    overlapping = True
                    break
                else:
                    continue
            # If not overlap with any other gene was found append the
            # gene to the list
            if not overlapping:
                self.result_genes.append(query_gene)

    def _already_seen(self, gene):
        key = ":".join([gene.seq_id, gene.gene_id, gene.name, str(gene.start), 
                        str(gene.end), gene.strand])
        if key in self._seen_genes:
            return(True)
        else:
            self._seen_genes[key] = 1
            return(False)

    def _have_sufficient_overlap(self, gene_1, gene_2):
        overlap = helpers.overlap(
            gene_1.start, gene_1.end, gene_2.start, gene_2.end)
        if overlap <= 0:
            return(False)
        if not self.min_overlap_percentage:
            return(True)
        else:
            if ((self._overlap_percentage(gene_1, overlap)
                 >= self.min_overlap_percentage) and 
                (self._overlap_percentage(gene_2, overlap)
                 >= self.min_overlap_percentage)):
                return(True)
        return(False)

    def _overlap_percentage(self, gene, overlap):
        try: 
            return(float(overlap) / float(gene.len()) * 100.0)
        except ZeroDivisionError:
            sys.stderr.write("Zero length for gene '%s/%s (start and end: %s)'\n" % (
                    gene.gene_id, gene.name, gene.start))

    def _merge_genes(self, gene_1, gene_2):
        start = min([gene_1.start, gene_2.end])
        end = max([gene_1.start, gene_2.end])
        name = "_merged_with_".join([gene_1.name, gene_2.name])
        gene_id = "_merged_with_".join([gene_1.gene_id, gene_2.gene_id])
        return(Gene(gene_1.seq_id, gene_id, name, start, end, gene_1.strand))
