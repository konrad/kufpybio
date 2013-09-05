class IGRFinder(object):

    # TODO 0- or 1-based
    def __init__(self, coordinate_system_base=0):
        self.coordinate_system_base = coordinate_system_base
        self._in_gene = "g"
        self._in_irg = "i"

    def find_igrs(self, gene_list, element_length):
        """ """
        self._init_nucleotide_list(element_length)
        self._mark_nucleotides_in_genes(gene_list)
        self._determine_igr_positions()
        return self.igr_positions

    def _init_nucleotide_list(self, element_length):
        """ Initiate the nucleotide list - for each position it is
        assumed it is intergenic..

        """
        self._nucleotides = [self._in_irg] * element_length

    def _mark_nucleotides_in_genes(self, gene_list):
        """Mark the nucleotides that are located in genes."""
        for gene in gene_list:
            self._nucleotides[gene.start:gene.end] = [self._in_gene] * (
                gene.end - gene.start)

    def _determine_igr_positions(self):
        self.igr_positions = []
        curr_igr_start = None
        # If the nucleotide list starts with an IGR the start position
        # needs to be set to 0. If this is not the case it will get a
        # value later.
        if self._nucleotides[0] == self._in_irg:
            curr_igr_start = 0
        prev_nucl_val = self._nucleotides[0]
        for index, nucl_val in enumerate(self._nucleotides):
            # Only transitions from IGR to gene and gene to IGR
            # require an action.
            if prev_nucl_val == self._in_irg and nucl_val == self._in_gene:
                self.igr_positions.append((curr_igr_start+1, index-1))
            elif prev_nucl_val == self._in_gene and nucl_val == self._in_irg:
                curr_igr_start = index
            elif prev_nucl_val == 0 and index == len(self._nucleotides-1):
                self.igr_positions.append((0,index))
            prev_nucl_val = nucl_val
        if prev_nucl_val == self._in_irg:
            self.igr_positions.append((curr_igr_start+1, len(self._nucleotides)-1))
