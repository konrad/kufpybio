# http://www.kegg.jp/kegg/docs/keggapi.html

import os
import urllib.request
import csv

class KEGGRESTAPI(object):

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
        return(self._second_column_to_list(data))

    def linked_reactions(self, entity_id):
        """entity_id can be e.g. a gene id, a KO id or as pathway id"""
        data = self._get_data(
            "%s/%s_linked_reactions.kegg", "link/reaction/%s", entity_id)
        return(self._second_column_to_list(data))

    def _get_data(self, path_template, url_template, entity_id):
        file_path = self._file_path(path_template, entity_id)
        if not os.path.exists(file_path):
            self._retrive_data(url_template, entity_id, file_path)
        return(open(file_path).read())

    def _retrive_data(self, url_template, entity_id, file_path):
        print(self._base_url + url_template % (entity_id))
        data = urllib.request.urlopen(
            self._base_url + url_template % (entity_id)).read()
        data_fh = open(file_path, "wb")
        data_fh.write(data)
        data_fh.close()

    def _file_path(self, path_template, entity_id):
        return(path_template % (self._download_folder, entity_id))

    def _rest_url(self, url_template, entity_id):
        return(url_template % (self._base_url, entity_id))

    def _second_column_to_list(self, data):
        return(
        [line.strip().split("\t")[1] for line in
         filter(lambda line: len(line) > 0, data.split("\n"))])
