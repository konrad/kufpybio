# https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch

import os
from xml.etree.ElementTree import ElementTree
import urllib.request

class NCBIEfetch(object):

    def __init__(self, email_address="", download_folder="ncbi_files"):
        self._email_address = email_address
        self._base_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
        self._download_folder = download_folder
        if not os.path.exists(self._download_folder):
            os.makedirs(self._download_folder)

    def gene(self, gene_id):
        gene_file = self._file_path(gene_id)
        if not os.path.exists(gene_file):
            self._retrieve_gene_data(gene_id)
        return(self._gene_file_to_dict(gene_id))

    def _retrieve_gene_data(self, gene_id):
        data = urllib.request.urlopen(self._gene_xml_url(gene_id)).read()
        gene_fh = open(self._file_path(gene_id), "wb")
        gene_fh.write(data)
        gene_fh.close()

    def _file_path(self, gene_id):
        return("%s/%s.xml" % (self._download_folder, gene_id))

    def _gene_xml_url(self, gene_id):
                return("%sdb=gene&id=%s&rettype=xml&email=%s" % (
                    self._base_url, gene_id, self._email_address))

    def _gene_file_to_dict(self, gene_id):
        gene_dict = {}
        tree = ElementTree()
        tree.parse(self._file_path(gene_id))
        gene_dict["uniprot_id"] = None
        gene_dict["kegg_id"] = None
        for db_tag in tree.findall(".//Dbtag"):
            dbtag_db = db_tag.find(".//Dbtag_db")
            object_id_str = db_tag.find(".//Object-id_str")
            if dbtag_db.text == "UniProtKB/TrEMBL" and object_id_str != None:
                gene_dict["uniprot_id"] = object_id_str.text
            elif dbtag_db.text == "KEGG" and object_id_str != None:
                gene_dict["kegg_id"] = object_id_str.text
        return(gene_dict)
