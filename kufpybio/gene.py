class Gene(object):

    def __init__(self, seq_id, gene_id, name, start, end, strand, 
                 feature=None, extra=None):
        """A Gene
        
        Gene is meant as an abstraction for gene, transcript, exon or
        any other entity.

        seq_id - identifier of the harboring chromosome, plasmid, etc.
        gene_id - a unique identifier of the gene
        name - name of the gene
        start - start position of the gene
        end - end postion of the gene
        strand - the strand (+ or -)
        feature - e.g. exon, gene, transcript
        extra - any other information that should be associated

        There is no assumptions / restrictions if the gene coordinates
        are in a 0-based or 1-based system.
        
        """
        self.seq_id = seq_id        
        self.gene_id = gene_id
        self.name = name
        self.start, self.end = sorted([int(start), int(end)])
        self.strand = strand
        if not feature is None:
            self.feature = feature
        if not extra is None:
            self.extra = extra

    def len(self):
        return self.end - self.start + 1
    
    def __repr__(self):
        return("Gene Id: %s\n" 
               "Gene name: %s\n"
               "Sequence id: %s\n" 
               "Strand: %s\n" 
               "Start: %s\n" 
               "End: %s" % (self.gene_id, self.name, self.seq_id, 
                            self.strand, self.start, self.end))
