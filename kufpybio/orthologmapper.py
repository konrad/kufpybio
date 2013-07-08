from subprocess import call
import os.path
import sys
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML

class OrthologMapper(object):
    """Generate one-to-one ortholog mappings based on BLAST results.

    Currently this is performed only in a very strict one-to-one
    fashion.

    Example usage:
    orthology_mapper = OrthologMapper(
        "input/NC_009641.fa", "input/NC_010079.fa", "nucl", "output", 
        makeblastbin="bin/makeblastdb", blastnbin="bin/blastn")
    orthology_mapper.blast()
    for entry_1, entry_2 in orthology_mapper.ortholog_pairs():
        print("\t".join([entry_1, entry_2]))
    """
    
    def __init__(self, fasta_file_1, fasta_file_2, seq_type, output_folder,
                 max_evalue=0.1, min_rel_length=0.6, makeblastbin="makeblastdb", 
                 blastpbin="blastp", blastnbin="blastn", override=False):
        self._fasta_file_1 = fasta_file_1
        self._fasta_file_2 = fasta_file_2
        self._seq_type = seq_type
        self._output_folder = output_folder
        self._max_evalue = float(max_evalue)
        self._min_rel_length = float(min_rel_length)
        self._makeblastbin = makeblastbin
        self._blastpbin = blastpbin
        self._blastnbin = blastnbin
        self._override = override
        self._no_ortholog_str = "no_ortholog"
        self._check_args_validity()
        
    def _check_args_validity(self):
        if self._seq_type not in ["nucl", "prot"]:
            sys.stderr.write("Sequence type must be 'nucl' or 'prot'\n")
            sys.exit(1)

    def blast(self):
        self._generate_blast_db(self._fasta_file_1)
        self._generate_blast_db(self._fasta_file_2)
        self._blast(self._fasta_file_1, self._fasta_file_2)
        self._blast(self._fasta_file_2, self._fasta_file_1)

    def ortholog_pairs(self):
        best_hits_1_in_2 = self._read_blast_file(
            self._blast_result_file(self._fasta_file_1, self._fasta_file_2))
        best_hits_2_in_1 = self._read_blast_file(
            self._blast_result_file(self._fasta_file_2, self._fasta_file_1))
        seen_list_2_entries = {}
        # Return all pairs by given the entries from fasta file 1
        # and their associated orthologs
        for query_id, hit_id in best_hits_1_in_2.items():
            # Test if hits are reciprocal best hits
            if query_id == best_hits_2_in_1.get(hit_id):
                seen_list_2_entries[hit_id] = 1
                yield((query_id, hit_id))
            else:
                yield((query_id, self._no_ortholog_str))
        # Now return all the entries from fasta file 2 that have not
        # been returned before
        for seq_id in best_hits_2_in_1.keys():
            if seq_id not in seen_list_2_entries:
                yield((self._no_ortholog_str, seq_id))
        
    def _generate_blast_db(self, fasta_file):
        call([self._makeblastbin, "-in", fasta_file, "-dbtype", self._seq_type])

    def _blast(self, fasta_db, fasta_query):
        blast_result_file = self._blast_result_file(fasta_query, fasta_db)
        if self._override is False and os.path.exists(blast_result_file):
            sys.stderr.write("Blast result file \"%s\" exists. "
                             "Skip blasting.\n" % blast_result_file)
            return 
        fasta_db_basename = os.path.basename(fasta_db)
        fasta_query_basename = os.path.basename(fasta_query)
        if self._seq_type == "prot":
            blast_bin = self._blastpbin
        elif self._seq_type == "nucl":
            blast_bin = self._blastnbin
        blast_cli = NcbiblastpCommandline(
            cmd="%s" % blast_bin,
            out=blast_result_file,
            outfmt=5,
            query=fasta_query,
            db=fasta_db)
        print("Running \"%s\"." %blast_cli)
        stdout, stderr = blast_cli()

    def _blast_result_file(self, fasta_query, fasta_db):
        return "%s/Search_%s_in_%s.xml" % (
            self._output_folder, self._basename(fasta_query),
            self._basename(fasta_db))

    def _basename(self, fasta_file):
        return os.path.basename(fasta_file)

    def _read_blast_file(self, blast_result_file):
        """Read the hits for each query. 

        Accept only the one best hit (based on bit score) for each
        query.
        
        """
        seq_ids_and_hits = {}
        for result in NCBIXML.parse(open(blast_result_file)):
            max_alignment_length = 0
            seq_ids_and_hits.setdefault(result.query, [])
            best_hit = self._no_ortholog_str
            best_hit_bit_score = 0
            for alignment in result.alignments:
                for hsp in alignment.hsps:
                    rel_alignment_len = self._rel_alignment_len(
                        result.query_length, hsp.align_length)
                    if (hsp.expect <= self._max_evalue and 
                        rel_alignment_len >=  self._min_rel_length and
                        hsp.bits > best_hit_bit_score):
                        best_hit = alignment.hit_def
                        best_hit_bit_score = hsp.bits
            seq_ids_and_hits[result.query] = best_hit
        return(seq_ids_and_hits)

    def _rel_alignment_len(self, query_length, align_length):
        return(float(align_length) / float(query_length))
