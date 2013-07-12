from Bio import Seq
from Bio.Alphabet import DNAAlphabet

class SubSeqExtractor(object):
    """Extract subsequences from a given DNA sequence"""

    def __init__(self, seq, coordinate_start=1, include_pos=True):
        """
        seq - sequence string
        coordinate_start - 0 for 0-based system; 1 for 1-based system
        include_pos - If True then the given position is given;
                      otherwise only the sequence before or after the
                      position are returned. Is later down- and
                      upstream sequence are requeste this argument has
                      has no influence and the query base will be
                      included.
        """
        self._coordinate_start = coordinate_start
        self._include_pos = include_pos
        self._seq = Seq.Seq(seq, DNAAlphabet())

    def extract(self, pos, upstream=0, downstream=0, rev_strand=False):
        """
        pos - postiion of the nucleotide
        upstream - length of the sequence used upstream
        downstream - length of the sequence used downstream
        rev_strand - True or False - If the feature is located on the
                     reverse stand. This implies that the upstream and
                     downstream selection is reversed and that the
                     reveser complement sequence is returned.
        """
        assert self._coordinate_start == 0 or self._coordinate_start == 1
        assert pos >= 0
        if self._coordinate_start == 1:
            assert pos >= 1
        if rev_strand is True:
            upstream, downstream = downstream, upstream
        start = pos - upstream - 1
        end = pos + downstream
        if self._include_pos is False and downstream != 0 and upstream != 0:
            pass
        elif self._include_pos is False and downstream == 0 and upstream == 0:
            start += 1
        elif self._include_pos is False and upstream != 0:
            end -= 1
        elif self._include_pos is False and downstream != 0:
            start += 1
        if self._coordinate_start == 0:
            start += 1
            end += 1
        assert start >= 0
        assert end-1 < len(self._seq)
        if rev_strand is True:
            return str(self._seq[start:end].reverse_complement())
        return str(self._seq[start:end])
