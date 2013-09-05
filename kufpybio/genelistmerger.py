import os
import sys
import kufpybio.helpers as helpers
import sqlite3
from kufpybio.gene import Gene

class GeneListMerger(object):

    def __init__(self, min_overlap_percentage=None, perform_gene_merging=True):
        # This flag influeces the general behavior. It it is set to
        # True gene will be merged i.e. the minimal star position
        # value and maximum end position taken. If it is False new
        # genes that are overlapping with old ones will be discarded.
        self._perform_gene_merging = perform_gene_merging
        self._min_overlap_percentage = min_overlap_percentage
        self._tmp_db_file = 'tmp.db'
        self._conn = sqlite3.connect(self._tmp_db_file)
        self._cur = self._conn.cursor()
        self._create_table()

    def add_genes(self, genes_to_add):
        for gene in genes_to_add:
            self._test_and_add(gene)

    def merged_genes(self):
        return [self._sql_row_to_gene(row) for row 
                in self._cur.execute("SELECT * FROM genes")]

    def cleanup(self):
        self._conn.close()
        os.remove(self._tmp_db_file)

    def _create_table(self):
        try:
            self._cur.execute(
                "create table genes (row_id integer primary key, name text, "
                "chrom text, strand text, start int, end int)")
        except sqlite3.OperationalError:
            os.remove(self._tmp_db_file)

    def _test_and_add(self, gene):
        rows = list(self._cur.execute(
                "SELECT * FROM genes WHERE chrom=? AND strand=? AND "
                "start <= ? AND end >= ?",
                (gene.seq_id, gene.strand, gene.end, gene.start)))
        if len(rows) == 0:
            self._store_gene_in_db(gene)
        else:
            # If the self._merge_genes is False the new gene is discarded
            if self._perform_gene_merging is False:
                return
            elif len(rows) == 1:
                self._try_gene_merge(gene, rows[0])
            else:
                self._try_multi_gene_merge(gene, rows)

    def _try_multi_gene_merge(self, gene, rows):
        overlapping_genes = [self._sql_row_to_gene(row) for row in rows]
        genes_to_merge = []
        rows_to_remove = []
        for overlapping_gene, row in zip(overlapping_genes, rows):
            if not self._have_sufficient_overlap(
                gene, overlapping_gene) is True:
                continue
            genes_to_merge.append(overlapping_gene)
            rows_to_remove.append(row)
        for row in rows_to_remove:
            self._remove_row(row)
        start = min([gene.start] + [gene.start for gene in genes_to_merge])
        end = max([gene.end] + [gene.end for gene in genes_to_merge])
        gene_id = "_merged_with_".join(
            [gene.gene_id for gene in genes_to_merge] + [gene.gene_id])
        self._store_gene_in_db(
            Gene(gene.seq_id, gene_id, gene_id, start, end, gene.strand))

    def _try_gene_merge(self, gene, row):
        overlapping_gene = self._sql_row_to_gene(row)
        if not self._have_sufficient_overlap(
            gene, overlapping_gene) is True:
            return
        if self._genes_are_identical(gene, overlapping_gene) is True:
            return
        start = min(gene.start, overlapping_gene.start)
        end = max(gene.end, overlapping_gene.end)
        gene_id = "%s_merged_with_%s" % (
            overlapping_gene.gene_id, gene.gene_id)
        self._remove_row(row)
        self._store_gene_in_db(
            Gene(gene.seq_id, gene_id, gene_id, start, end, gene.strand))

    def _store_gene_in_db(self, gene):
        self._cur.execute("INSERT INTO genes VALUES (NULL, ?, ?, ?, ?, ?)",
                              (gene.gene_id, gene.seq_id, gene.strand, 
                               gene.start, gene.end))
        self._conn.commit()

    def _genes_are_identical(self, gene_1, gene_2):
        if (gene_1.gene_id == gene_2.gene_id and 
            gene_1.start == gene_2.end and
            gene_1.end == gene_2.end):
            return True
        return False

    def _remove_row(self, row):
        self._cur.execute("DELETE FROM genes WHERE row_id=%s" % (row[0]))
        self._conn.commit()

    def _sql_row_to_gene(self, row):
        return Gene(row[2], row[1], row[1], row[4], row[5], row[3])

    def _have_sufficient_overlap(self, gene_1, gene_2):
        overlap = helpers.overlap(
            gene_1.start, gene_1.end, gene_2.start, gene_2.end)
        if overlap <= 0:
            return False
        if not self._min_overlap_percentage:
            return True
        else:
            if ((self._overlap_percentage(gene_1, overlap)
                 >= self._min_overlap_percentage) and
                (self._overlap_percentage(gene_2, overlap)
                 >= self._min_overlap_percentage)):
                return True
        return False

    def _overlap_percentage(self, gene, overlap):
        try:
            return float(overlap) / float(gene.len()) * 100.0
        except ZeroDivisionError:
            sys.stderr.write("Zero length for gene '%s/%s (start and end: %s)'\n" % (
                    gene.gene_id, gene.name, gene.start))

    def _merge_genes(self, gene_1, gene_2):
        start = min([gene_1.start, gene_2.end])
        end = max([gene_1.start, gene_2.end])
        name = "_merged_with_".join([gene_1.name, gene_2.name])
        gene_id = "_merged_with_".join([gene_1.gene_id, gene_2.gene_id])
        return Gene(gene_1.seq_id, gene_id, name, start, end, gene_1.strand)
