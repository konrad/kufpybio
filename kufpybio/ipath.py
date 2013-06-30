from httplib2 import Http
from urllib.parse import urlencode
import re
import urllib.request

class IPathDownload(object):

    """

    More info at
    http://pathways.embl.de/mapping.cgi
    http://pathways.embl.de/help.html

    Example:
    from ipath import IPathDownload
    kegg_ids = ["R01827 #009900 W15 0.8", "R01281"]
    ipath_download = IPathDownload(kegg_ids, output_file="my_pathway_viz.svg")
    ipath_download.get_ipath()

    """

    def __init__(self, selection=[], output_file="ipath.svg"):
        self._url = "http://pathways.embl.de/mapping.cgi"
        self._form_values = {
            "selection" : "",
            "opacity" : 1,
            "default_width" : 3,
            "default_radius" : 7,
            "keep_colors" : "",
            "default_color" : "#666666",
            "background_color" : "#ffffff",
            "whole_pathways" : 0,
            "query_reactions" : "",
            "tax_filter" : "",
            "map_type" : "svg", # svg, png, interactive
            "include_metabolic" : 1,
            "include_regulatory" : 0,
            "include_secondary" : 0,
            "png_dpi" : 72 # only needed for png output
            }
        self._output_file = output_file
        self._selection = selection

    def set_form_values(self, opacity=None, default_width=None,
                        default_radius=None, keep_colors=None,
                        default_color=None, background_color=None,
                        whole_pathways=None, query_reactions=None,
                        tax_filter=None, map_type=None,
                        include_metabolic=None):
        for parameter, value in locals().items():
            if parameter == "self":
                continue
            if value != None:
                self._form_values[parameter] = value

    def get_ipath(self):
        self._process_selection()
        response, content = self._sent_form()
        svg_url = self._extract_svg_url(content)
        self._download_svg_file(svg_url)

    def _process_selection(self):
        assert len(self._selection) != 0
        assert type(self._selection) == list
        self._form_values["selection"] = "\n".join(self._selection)

    def _extract_svg_url(self, content):
        url_patter = re.compile('href="(.+)">here<')
        return re.search(url_patter, content.decode('utf-8')).groups(0)[0]

    def _sent_form(self):
        http = Http()
        data = urlencode(self._form_values)
        return http.request(self._url, "POST", data)

    def _download_svg_file(self, svg_url):
        svg_content = urllib.request.urlopen(svg_url).read()
        fh = open(self._output_file, "w")
        fh.write(svg_content.decode('utf-8'))
        fh.close()
