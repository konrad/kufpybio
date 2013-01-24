import unittest
import sys
from kufpybio.ipath import IPathDownload

class TestIPathDownload(unittest.TestCase):

    def setUp(self):
        self._ipath_download = IPathDownload(selection=[
            "R01827 #009900 W15 0.8", "R01281"])
        self._initial_values = {
            "selection" : "", "opacity" : 1, "default_width" : 3,
            "default_radius" : 7,"keep_colors" : "",
            "default_color" : "#666666", "background_color" : "#ffffff",
            "whole_pathways" : 0, "query_reactions" : "", "tax_filter" : "",
            "map_type" : "svg", "include_metabolic" : 1,
            "include_regulatory" : 0, "include_secondary" : 0, "png_dpi" : 72}

    def test_initiation(self):
        self.assertEqual(
            self._ipath_download._url, "http://pathways.embl.de/mapping.cgi")
        self.assertDictEqual(
            self._ipath_download._form_values, self._initial_values)

    def test_process_selection(self):
        self._ipath_download._process_selection()
        self.assertEqual(self._ipath_download._form_values["selection"],
                         "R01827 #009900 W15 0.8\nR01281")

    def test_set_form_values(self):
        self._ipath_download.set_form_values(opacity=5)
        expected_values = self._initial_values.copy()
        expected_values["opacity"] = 5
        self.assertDictEqual(self._ipath_download._form_values, expected_values)

if __name__ == "__main__":
    unittest.main()
