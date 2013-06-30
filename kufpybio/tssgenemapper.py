"""Classes for associating TSS and genes."""

# Some commonly used strings:
# Strand
plus_str = "+"
minus_str = "-"
# TSS types
orphan_str = "orphan"
primary_str = "primary"
secondary_str = "secondary"
# TSS locations/orientations
sense_str = "sense"
antisense_str = "antisense"
loc_5_prime_str = "5' region"
loc_internal_str = "internal"
loc_antisense_str = "antisense"

class TSSGeneMapper(object):
    """Find the associated genes for TSS (Transcription start sites).

    The following types of TSS-gene constellation are discriminated:


        secondary   primary  internal
        +--->       +--->    +---> 
        |           |        |
    ----+-----------+---=====+=====>------
                             |       
                         <---+
                         antisense

    5' TSS must have the same orientation as the associated genes and
    must be located in a certain range upstream of the gene (see value
    of max_dist_5_prime) including the first based of the gene. This
    means that leaderless genes will have an UTR length of 0. Internal
    genes must have the same orientation as the associated gene and
    must be located inside of the gene. Anti-sense TSS must have the
    opposite orientation as the associated gene and must be located
    inside the gene or upstream or downstream of the gene in a certain
    range (see value of max_dist_antisense).

    5' TSS are classified in primary and secondary TSS. The TSS that
    is in closest proximity to the gene is the primary TSS, the
    remaining ones are considered as secondary TSS.

    """

    def __init__(self, max_dist_5_prime=300, max_dist_antisense=100):

        """
        Arguments:
        - max_dist_5_prime: The maximaum distance a 5' TSS can have
          from its associated gene.
        - max_dist_antisense: The range upstream and downstream of a
          gene in which a anti-sense TSS can be located.

        """
        self._max_dist_5_prime = max_dist_5_prime
        self._max_dist_antisense = max_dist_antisense

    def map_tss(self, tss_list, gene_list):
        """Perform the TSS to gene mapping 

        Returns a list of TSS, their associated genes, and the type of TSS.

        Arguments:
        - tss_list: must be a list of TSS (transcriptin start sites)
          with the two features tss.pos (position as integer) and
          tss.strand ("+" or "-").
        - gene_list: must be a list of genes with the following features:
          - gene.seq_id
          - gene.gene_id
          - gene.name
          - gene.start 
          - gene.end
          - gene.strand
          
        For both, plus and minus strand entries, it is assummed that
        gene.start < gene.end.

        A comment regarding performance: The searches are not
        optimized for speed - every gene is checked.

        """
        self.tss_and_hit_genes = {}
        self.genes_and_5_prime_tss = {}
        gene_list = gene_list
        for tss in tss_list:
            for gene in gene_list:
                self._check_tss_gene_associations(tss, gene)
            if tss not in self.tss_and_hit_genes:
                self.tss_and_hit_genes[tss] = orphan_str
        self._set_type_of_5_prime_tss()
        self._remove_multiple_associations()
        return self.tss_and_hit_genes

    def _set_type_of_5_prime_tss(self):
        for gene, distances_and_tss in self.genes_and_5_prime_tss.items():
            distance_sorted_tss = [
                tss for distance, tss in sorted(distances_and_tss)]
            self.tss_and_hit_genes[
                distance_sorted_tss[0]][gene]["tss_type"] = primary_str
            if len(distance_sorted_tss) > 1:
                for secondary_tss in distance_sorted_tss[1:]:
                    self.tss_and_hit_genes[secondary_tss][
                        gene]["tss_type"] = secondary_str

    def _check_tss_gene_associations(self, tss, gene):
        location = None
        distance = None
        if self._has_5_prime_association(tss, gene):
            location = loc_5_prime_str
            distance = self._5_prime_dist(tss, gene)
        elif self._has_internal_association(tss, gene):
            location = loc_internal_str
        elif self._has_antisense_association(tss, gene):
            location = loc_antisense_str
        if location:
            self.tss_and_hit_genes.setdefault(tss, {})
            self.tss_and_hit_genes[tss][gene] = {
                "location" : location, "distance" : distance, 
                "tss_type" : None}
            # If this is a TSS in the 5' region record the
            # distance. This will later be used to determine the type
            # of TSS.
            if type(distance) == int: # can be 0 or any positive integer:
                self.genes_and_5_prime_tss.setdefault(gene, [])
                self.genes_and_5_prime_tss[gene].append([distance, tss])
            
    def _has_5_prime_association(self, tss, gene):
        """Test for 5' prime association

        Any TSS including the first nucleotide (leaderless) in 5'
        should return True.
        
        """
        if tss.strand != gene.strand:
            return
        dist = self._5_prime_dist(tss, gene)
        if (tss.strand == gene.strand and 
            dist >= 0 and dist <= self._max_dist_5_prime):
            return True

    def _has_internal_association(self, tss, gene):
        """Test for internal location

        Any TSS from the second to the last base should return
        True. The first base is 5_prime_associated (leaderless TSS).

        """
        if tss.strand != gene.strand:
            return
        if (gene.strand == plus_str and tss.pos > gene.start and 
            tss.pos <= gene.end):
            return True
        if (gene.strand == minus_str and tss.pos >= gene.start and 
            tss.pos < gene.end):
            return True

    def _has_antisense_association(self, tss, gene):
        """Test for antisense association

        """
        if (tss.strand != gene.strand and 
            tss.pos >= gene.start - self._max_dist_antisense and 
            tss.pos <= gene.end + self._max_dist_antisense):
            return True

    def _5_prime_dist(self, tss, gene):
        """Calculate the TSS to gene start distance."""
        if gene.strand == plus_str:
            return gene.start - tss.pos
        elif gene.strand == minus_str:
            return tss.pos - gene.end

    def _remove_multiple_associations(self):
        """Remove multple 5' associations

        Until now each TSS-gene association was done without
        considering other genes. Here certain associations are
        removed. If a TSS is associated with a gene that is short than
        self._max_dist_5_prime is can be theoretically connected to a
        gene that is located downstream of the first gene. In such a
        case the second gene association does not make sense and needs
        to be removed. This mean only the 5' TSS with the lowest
        distance value is kept.

        """
        for tss, hit_genes in self.tss_and_hit_genes.items():
            if len(hit_genes) == 1 or hit_genes == orphan_str:
                continue
            min_dist = None
            closest_gene = None
            filtered_hit_genes = {}
            for hit_gene in hit_genes.keys():
                if hit_genes[hit_gene]["location"] == loc_5_prime_str:
                    if (not min_dist) or (
                        hit_genes[hit_gene]["distance"] < min_dist):
                        closest_gene = hit_gene
                        min_dist = hit_genes[hit_gene]["distance"]
                else:
                    filtered_hit_genes[hit_gene] = hit_genes[hit_gene]
            if closest_gene:
                filtered_hit_genes[closest_gene] = hit_genes[closest_gene]
            self.tss_and_hit_genes[tss] = filtered_hit_genes

class TSSGeneFormatter(object):

    binary_format_header = ["Primary", "Secondary", "Internal", "Antisense"]

    def tss_features_binary_format(self, tss_features):
        """Returns the TSS information in a binary format compatible
        with previous publications.

        """
        # Orphan
        if tss_features == orphan_str:
            return ["0", "0", "0", "0"]
        # 5' ...
        elif (tss_features["location"] == loc_5_prime_str):
            # ... primary
            if  tss_features["tss_type"] == primary_str:
                return ["1", "0", "0", "0"]
            # ... secondary
            elif tss_features["tss_type"] == secondary_str:
                return ["0", "1", "0", "0"]
        # Internal
        elif (tss_features["location"] == loc_internal_str):
            return ["0", "0", "1", "0"]
        # Antisense 
        elif (tss_features["location"] == loc_antisense_str):
            return ["0", "0", "0", "1"]

    def tss_features_string(self, tss_features):
        if tss_features == orphan_str:
            return orphan_str
        elif tss_features["tss_type"]:
            return "%s - %s" % (
                    tss_features["location"], tss_features["tss_type"])
        else:
            return tss_features["location"]
