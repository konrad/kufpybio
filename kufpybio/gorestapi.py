# https://www.ebi.ac.uk/QuickGO/WebServices.html

import os
import csv

import restapi

class GORESTAPI(restapi.RESTAPI):

    def __init__(self, download_folder="go_files"):
        self._download_folder = download_folder
        self._base_url = "https://www.ebi.ac.uk/QuickGO/GTerm?"
        self._create_download_folder()

    def go_term_information_xml(self, go_id):
        """e.g. GO:0003824"""
        data = self._get_data(
            "%s/%s.xml", "id=%s&format=oboxml", go_id)
        return(data)
