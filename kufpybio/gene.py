class Gene(object):

    def __init__(self, seq_id, gene_id, name, start, end, strand, extra=None):
        """A Gene
        
        seq_id - identifier of the harboring chromosome, plasmid, etc.
        gene_id - a unique identifier of the gene
        name - name of the gene
        start - start position of the gene
        end - end postion of the gene
        strand - the strand (+ or -)
        extra - any other information that should be associated

        There is no assumptions / restrictions if the gene coordinates
        are in a 0-based or 1-based system.
        
        """
        self.seq_id = seq_id        
        self.gene_id = gene_id
        self.name = name
        self.start, self.end = sorted([int(start), int(end)])
        self.strand = strand
        if extra:
            self.extra = extra

    def len(self):
        return(self.end - self.start)
