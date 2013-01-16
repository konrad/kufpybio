# http://www.kegg.jp/kegg/docs/keggapi.html

import os
import urllib.request
import csv

class KEGGAPI(object):

    def __init__(self, download_folder="kegg_files"):
        self._download_folder = download_folder
        self._base_url = "http://rest.kegg.jp/"
        if not os.path.exists(self._download_folder):
            os.makedirs(self._download_folder)

    def pathway(self, map_number):
        map_file = self._map_file_path(map_number)
        if not os.path.exists(map_file):
            self._retrieve_map(map_number)
        return(self._file_to_dict(map_file))

    def _map_file_path(self, map_number):
        return("%s/map%s.kegg" % (self._download_folder, map_number))

    def _retrieve_map(self, map_number):
        map_file = self._map_file_path(map_number)
        data = urllib.request.urlopen(
            self._base_url + "/get/pathway:map%s" % map_number).read()
        map_fh = open(map_file, "wb")
        map_fh.write(data)
        map_fh.close()

    def _file_to_dict(self, kegg_file):
        data_dict = {}
        for line in open(kegg_file):
            if line.startswith("///"):
                continue
            split_line = line.strip().split()
            data_dict[split_line[0]] = " ".join(split_line[1:])
        print(data_dict)
