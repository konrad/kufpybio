# http://www.uniprot.org/faq/28

"""This serves just the download of files.

Parsing can be done using Biopython:

from Bio import SeqIO
uniprot_entry = SeqIO.read(open("uniprot_files/A8Z556.xml"), "uniprot-xml")

"""

import os
import urllib.request

from restapi import RESTAPI

class UniprotRESTAPI(RESTAPI):

    def __init__(self, download_folder="uniprot_files"):
        self._download_folder = download_folder
        self._base_url = "http://www.uniprot.org/uniprot/"
        self._create_download_folder()

    def uniprot_protein_xml(self, uniprot_id):
        """e.g. Q5FJ41"""
        return(self._get_data("%s/%s.xml", "%s.xml", uniprot_id))
