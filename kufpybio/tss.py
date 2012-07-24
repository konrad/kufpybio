class TSS(object):

    def __init__(self, seq_id, pos, strand, extra=None):
        """A transcription start site
        
        seq_id - identifier of the harboring chromosome, plasmid, etc.
        pos - position of the
        strand - the strand (+ or -)
        extra - any other information that should be associated

        There is no assumptions / restrictions if the TSS position is
        in a 0-based or 1-based coordinate system.
        """

        self.seq_id = seq_id
        self.pos = int(pos)
        self.strand = strand
        if extra:
            self.extra = extra
