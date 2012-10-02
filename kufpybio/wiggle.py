import csv

class WiggleParser(object):
    """

    Warning - this does not implement the full specification!

    """

    def entries(self, input_fh):
        track_name = None
        chrom_name = None
        span = None
        pos_value_pairs = []
        for line in input_fh:
            row = line[:-1].split()
            if len(row) == 0:
                continue
            if row[0].startswith("track"):
                track_name = self._track_name(row)
            elif row[0].startswith("variableStep"):
                if chrom_name:
                    prev_chrom_name = chrom_name
                    prev_span = span
                    prev_pos_value_pairs = pos_value_pairs
                    chrom_name = self._chrom_name(row)
                    span = None
                    pos_value_pairs = []
                    yield(WiggleEntry(
                            track_name, prev_chrom_name, prev_span, prev_pos_value_pairs))
                else:
                    chrom_name = self._chrom_name(row)
            else:
                pos_value_pairs.append([int(row[0]), float(row[1])])

    def _chrom_name(self, row):
        return(self._attrs_and_values(row)["chrom"])
                
    def _track_name(self, row):
        return(self._attrs_and_values(row)["name"])

    def _attrs_and_values(self, row):
        attrs_and_values = {}
        for attr_and_value in row:
            if not "=" in attr_and_value:
                continue
            attr, value = attr_and_value.split("=")
            value = value.replace("\"", "")
            attrs_and_values[attr] = value
        return(attrs_and_values)

class WiggleEntry(object):
    
    def __init__(self, track_name, chrom_name, span, pos_value_pairs):
        self.track_name = track_name
        self.chrom_name = chrom_name
        self.span = span
        self.pos_value_pairs = pos_value_pairs
