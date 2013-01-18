# http://www.kegg.jp/kegg/docs/keggapi.html

import os
import csv

import restapi

class KEGGRESTAPI(restapi.RESTAPI):

    def __init__(self, download_folder="kegg_files"):
        self._download_folder = download_folder
        self._base_url = "http://rest.kegg.jp/"
        if not os.path.exists(self._download_folder):
            os.makedirs(self._download_folder)

    def reference_pathway(self, map_number):
        """e.g. 00630"""
        data = self._get_data(
            "%s/map%s.kegg", "get/pathway:map%s", map_number)
        return(data)

    def pathway(self, map_number):
        """e.g. sax00270"""
        data = self._get_data(
            "%s/%s.kegg", "get/pathway:%s", map_number)
        return(data)

    def linked_pathways(self, entity_id):
        """entity_id can be e.g. a gene id, a KO id or a reaction id"""
        data = self._get_data(
            "%s/%s_linked_pathways.kegg", "link/pathway/%s", entity_id)
        return(self._clean_list(self._second_column_to_list(data), "path:"))

    def linked_reactions(self, entity_id):
        """entity_id can be e.g. a gene id, a KO id or as pathway id"""
        data = self._get_data(
            "%s/%s_linked_reactions.kegg", "link/reaction/%s", entity_id)
        return(self._clean_list(self._second_column_to_list(data), "rn:"))

    def _clean_list(self, id_list, prefix_to_remove):
        return([item.replace(prefix_to_remove, "") for item in id_list])

    def _second_column_to_list(self, data):
        return(
        [line.strip().split("\t")[1] for line in
         filter(lambda line: len(line) > 0, data.split("\n"))])
