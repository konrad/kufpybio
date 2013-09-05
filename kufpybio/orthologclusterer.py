class OrthologClusterer(object):
    """Translates a list of ortholog pairs into ortholog clusters"""

    def __init__(self):
        self._orthologs_and_clusters = {}

    def process_pairs_to_clusters(self, ortholog_pairs):
        """Translate the ortholog pair information into clusters

        ortholog_pairs must be a list of tuple (the two orthologs) of
        tuples (the ortholog accession and the source name). I a
        gene/protein does not have a ortholog the tuple of the partner
        must be set to None.
        """
        for ortholog_tuple_1, ortholog_tuple_2 in ortholog_pairs:
            ortholog_1 = self._tuple_to_ortholog(ortholog_tuple_1)
            ortholog_2 = self._tuple_to_ortholog(ortholog_tuple_2)
            combined_cluster = self._combine_clusters(ortholog_1, ortholog_2)
            # Update the cluster information for all orthologs in the cluster
            for ortholog in combined_cluster:
                self._orthologs_and_clusters[ortholog.key] = combined_cluster

    def _tuple_to_ortholog(self, ortholog_tuple):
        if not ortholog_tuple == None:
            return Ortholog(ortholog_tuple[0], ortholog_tuple[1])
        else:
            # If instead of the tuple a None is given no ortholog is created
            return None

    def _combine_clusters(self, ortholog_1, ortholog_2):
        cluster_of_ortholog_1 = self._get_or_generate_cluster(ortholog_1)
        cluster_of_ortholog_2 = self._get_or_generate_cluster(ortholog_2)
        return cluster_of_ortholog_1.union(cluster_of_ortholog_2)

    def _get_or_generate_cluster(self, ortholog):
        if ortholog is not None:
            return self._orthologs_and_clusters.get(ortholog.key, set([ortholog]))
        else:
            return set([])

    def clusters(self):
        """Return a non-redundant list of cluster"""
        seen_keys = {}
        for cluster in self._orthologs_and_clusters.values():
            yield_cluster = True
            for ortholog in cluster:
                if ortholog.key in seen_keys:
                    yield_cluster = False
                    break
                seen_keys[ortholog.key] = None
            if yield_cluster is False:
                continue
            yield cluster

class Ortholog(object):
    
    def __init__(self, accession, source):
        self.accession = accession
        self.source = source
        self.key = self.accession + "_" + self.source

    def __repr__(self):
        return(self.key)
