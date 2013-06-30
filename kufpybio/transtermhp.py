class TranstermHPBAGReader(object):

    def entries(self, input_fh):
        for line in input_fh:
            if "NONE" in line:
                continue
            yield TranstermHPBAGEntry(line)

class TranstermHPBAGEntry(object):

    def __init__(self, line):
        row = line.strip().split()
        assert(len(row) == 14)
        self.gene_name = row[0]
        self.term_start = int(row[1])
        self.term_end = int(row[3])
        self.strand = row[4]
        self.hairpin_score = float(row[5])
        self.tail_score = float(row[6])
        self.term_seq = " ".join(row[7:12])
        self.term_confidence = int(row[12])
        self.approx_dist = int(row[13])
