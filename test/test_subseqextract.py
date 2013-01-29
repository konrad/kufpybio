import unittest
import sys
from kufpybio.subseqextract import SubSeqExtractor

class TestSubSeqExtractor1(unittest.TestCase):

    """Testing with
       - a 1-based coordinate system
       - the query base included in the extracted sequence

    """
    def setUp(self):
        self.sub_seq_extractor = SubSeqExtractor(
            "AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTATTT",
            coordinate_start=1,
            include_pos=True)

    def test_extract_single_pos(self):
         self.assertEqual(self.sub_seq_extractor.extract(10), "A")
         self.assertEqual(self.sub_seq_extractor.extract(11), "C")
         self.assertEqual(self.sub_seq_extractor.extract(20), "C")
         self.assertEqual(self.sub_seq_extractor.extract(21), "G")

    def test_extract_upstream_seq(self):
        self.assertEqual(self.sub_seq_extractor.extract(10, upstream=9),
                         "AAAAAAAAAA")
        self.assertEqual(
            len(self.sub_seq_extractor.extract(10, upstream=9)), 10)
        self.assertEqual(self.sub_seq_extractor.extract(20, upstream=9),
                         "CCCCCCCCCC")
        self.assertEqual(self.sub_seq_extractor.extract(25, upstream=9),
                         "CCCCCGGGGG")

    def test_extract_downstream_seq(self):
        self.assertEqual(
            self.sub_seq_extractor.extract(
                11, downstream=9),
            "CCCCCCCCCC")
        self.assertEqual(
            self.sub_seq_extractor.extract(21, downstream=9),
            "GGGGGGGGGG")
        self.assertEqual(
            self.sub_seq_extractor.extract(30, downstream=5),
            "GTTTAT")

    def test_extract_down_and_upstream_seq(self):
        self.assertEqual(self.sub_seq_extractor.extract(
            34, upstream=3, downstream=3), "TTTATTT")
        self.assertEqual(self.sub_seq_extractor.extract(
                10, upstream=5, downstream=5), "AAAAAACCCCC")

    def test_extract_single_pos_revese(self):
         self.assertEqual(
             self.sub_seq_extractor.extract(10, rev_strand=True), "T")
         self.assertEqual(
             self.sub_seq_extractor.extract(11, rev_strand=True), "G")
         self.assertEqual(
             self.sub_seq_extractor.extract(20, rev_strand=True), "G")
         self.assertEqual(
             self.sub_seq_extractor.extract(21, rev_strand=True), "C")

    def test_extract_upstream_seq_reverse(self):
        self.assertEqual(
            self.sub_seq_extractor.extract(
                11, upstream=9, rev_strand=True),
                "GGGGGGGGGG")
        self.assertEqual(
            self.sub_seq_extractor.extract(
                10, upstream=9, rev_strand=True),
                "GGGGGGGGGT")

    def test_extract_downstream_seq_reverse(self):
        self.assertEqual(
            self.sub_seq_extractor.extract(
                20, downstream=9, rev_strand=True), "GGGGGGGGGG")
        self.assertEqual(
            self.sub_seq_extractor.extract(
                11, downstream=9, rev_strand=True), "GTTTTTTTTT")
        self.assertEqual(
            self.sub_seq_extractor.extract(
                37, downstream=9, rev_strand=True), "AAATAAACCC")

    def test_extract_down_and_upstream_seq_reverse(self):
        self.assertEqual(self.sub_seq_extractor.extract(
            34, upstream=3, downstream=3, rev_strand=True),
            "AAATAAA")
        self.assertEqual(self.sub_seq_extractor.extract(
            10, upstream=5, downstream=5, rev_strand=True),
            "GGGGGTTTTTT")

class TestSubSeqExtractor2(unittest.TestCase):

    """Testing with
       - a 0-based coordinate system
       - the query base included in the extracted sequence

    """
    def setUp(self):
        self.sub_seq_extractor = SubSeqExtractor(
            "AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTATTT",
            coordinate_start=0,
            include_pos=True)

    def test_extract_single_pos(self):
         self.assertEqual(self.sub_seq_extractor.extract(9), "A")
         self.assertEqual(self.sub_seq_extractor.extract(10), "C")
         self.assertEqual(self.sub_seq_extractor.extract(19), "C")
         self.assertEqual(self.sub_seq_extractor.extract(20), "G")

    def test_extract_upstream_seq(self):
        self.assertEqual(self.sub_seq_extractor.extract(9, upstream=9),
                         "AAAAAAAAAA")
        self.assertEqual(
            len(self.sub_seq_extractor.extract(9, upstream=9)), 10)
        self.assertEqual(self.sub_seq_extractor.extract(19, upstream=9),
                         "CCCCCCCCCC")
        self.assertEqual(self.sub_seq_extractor.extract(24, upstream=9),
                         "CCCCCGGGGG")

    def test_extract_downstream_seq(self):
        self.assertEqual(
            self.sub_seq_extractor.extract(10, downstream=9),
            "CCCCCCCCCC")
        self.assertEqual(
            self.sub_seq_extractor.extract(20, downstream=9),
            "GGGGGGGGGG")

    def test_extract_down_and_upstream_seq(self):
        self.assertEqual(self.sub_seq_extractor.extract(
            33, upstream=3, downstream=3), "TTTATTT")
        self.assertEqual(self.sub_seq_extractor.extract(
                9, upstream=5, downstream=5), "AAAAAACCCCC")

    def test_extract_single_pos_revese(self):
         self.assertEqual(
             self.sub_seq_extractor.extract(9, rev_strand=True), "T")
         self.assertEqual(
             self.sub_seq_extractor.extract(10, rev_strand=True), "G")
         self.assertEqual(
             self.sub_seq_extractor.extract(19, rev_strand=True), "G")
         self.assertEqual(
             self.sub_seq_extractor.extract(20, rev_strand=True), "C")

    def test_extract_upstream_seq_reverse(self):
        self.assertEqual(
            self.sub_seq_extractor.extract(
                10, upstream=9, rev_strand=True),
                "GGGGGGGGGG")
        self.assertEqual(
            self.sub_seq_extractor.extract(
                9, upstream=9, rev_strand=True),
                "GGGGGGGGGT")

    def test_extract_downstream_seq_reverse(self):
        self.assertEqual(
            self.sub_seq_extractor.extract(
                19, downstream=9, rev_strand=True), "GGGGGGGGGG")
        self.assertEqual(
            self.sub_seq_extractor.extract(
                10, downstream=9, rev_strand=True), "GTTTTTTTTT")

    def test_extract_down_and_upstream_seq_reverse(self):
        self.assertEqual(self.sub_seq_extractor.extract(
            33, upstream=3, downstream=3, rev_strand=True),
            "AAATAAA")
        self.assertEqual(self.sub_seq_extractor.extract(
            9, upstream=5, downstream=5, rev_strand=True),
            "GGGGGTTTTTT")


class TestSubSeqExtractor3(unittest.TestCase):

    """Testing with
       - a 1-based coordinate system
       - the query base NOT included in the extracted sequence

    """
    def setUp(self):
        self.sub_seq_extractor = SubSeqExtractor(
            "AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTATTT",
            coordinate_start=1,
            include_pos=False)

    def test_extract_single_pos(self):
         self.assertEqual(self.sub_seq_extractor.extract(10), "")
         self.assertEqual(self.sub_seq_extractor.extract(11), "")
         self.assertEqual(self.sub_seq_extractor.extract(20), "")
         self.assertEqual(self.sub_seq_extractor.extract(21), "")

    def test_extract_upstream_seq(self):
        self.assertEqual(self.sub_seq_extractor.extract(10, upstream=9),
                         "AAAAAAAAA")
        self.assertEqual(
            len(self.sub_seq_extractor.extract(10, upstream=9)), 9)
        self.assertEqual(self.sub_seq_extractor.extract(20, upstream=9),
                         "CCCCCCCCC")
        self.assertEqual(self.sub_seq_extractor.extract(25, upstream=9),
                         "CCCCCGGGG")

    def test_extract_downstream_seq(self):
        self.assertEqual(
            self.sub_seq_extractor.extract(
                11, downstream=9),
            "CCCCCCCCC")
        self.assertEqual(
            self.sub_seq_extractor.extract(21, downstream=9),
            "GGGGGGGGG")
        self.assertEqual(
            self.sub_seq_extractor.extract(30, downstream=5),
            "TTTAT")

    def test_extract_down_and_upstream_seq(self):
        self.assertEqual(self.sub_seq_extractor.extract(
            34, upstream=3, downstream=3), "TTTATTT")
        self.assertEqual(self.sub_seq_extractor.extract(
                10, upstream=5, downstream=5), "AAAAAACCCCC")

    def test_extract_single_pos_revese(self):
         self.assertEqual(
             self.sub_seq_extractor.extract(10, rev_strand=True), "")
         self.assertEqual(
             self.sub_seq_extractor.extract(11, rev_strand=True), "")
         self.assertEqual(
             self.sub_seq_extractor.extract(20, rev_strand=True), "")
         self.assertEqual(
             self.sub_seq_extractor.extract(21, rev_strand=True), "")

    def test_extract_upstream_seq_reverse(self):
        self.assertEqual(
            self.sub_seq_extractor.extract(
                11, upstream=9, rev_strand=True),
                "GGGGGGGGG")
        self.assertEqual(
            self.sub_seq_extractor.extract(
                11, upstream=9, rev_strand=True),
                "GGGGGGGGG")

    def test_extract_downstream_seq_reverse(self):
        self.assertEqual(
            self.sub_seq_extractor.extract(
                20, downstream=9, rev_strand=True), "GGGGGGGGG")
        self.assertEqual(
            self.sub_seq_extractor.extract(
                11, downstream=9, rev_strand=True), "TTTTTTTTT")
        self.assertEqual(
            self.sub_seq_extractor.extract(
                37, downstream=9, rev_strand=True), "AATAAACCC")

    def test_extract_down_and_upstream_seq_reverse(self):
        self.assertEqual(self.sub_seq_extractor.extract(
            34, upstream=3, downstream=3, rev_strand=True),
            "AAATAAA")
        self.assertEqual(self.sub_seq_extractor.extract(
            10, upstream=5, downstream=5, rev_strand=True),
            "GGGGGTTTTTT")


class TestSubSeqExtractor4(unittest.TestCase):

    """Testing with
       - a 0-based coordinate system
       - the query base NOT included in the extracted sequence

    """
    def setUp(self):
        self.sub_seq_extractor = SubSeqExtractor(
            "AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTATTT",
            coordinate_start=0,
            include_pos=False)

    def test_extract_single_pos(self):
         self.assertEqual(self.sub_seq_extractor.extract(9), "")
         self.assertEqual(self.sub_seq_extractor.extract(10), "")
         self.assertEqual(self.sub_seq_extractor.extract(19), "")
         self.assertEqual(self.sub_seq_extractor.extract(20), "")

    def test_extract_upstream_seq(self):
        self.assertEqual(self.sub_seq_extractor.extract(9, upstream=9),
                         "AAAAAAAAA")
        self.assertEqual(
            len(self.sub_seq_extractor.extract(9, upstream=9)), 9)
        self.assertEqual(self.sub_seq_extractor.extract(19, upstream=9),
                         "CCCCCCCCC")
        self.assertEqual(self.sub_seq_extractor.extract(24, upstream=9),
                         "CCCCCGGGG")

    def test_extract_downstream_seq(self):
        self.assertEqual(
            self.sub_seq_extractor.extract(10, downstream=9),
            "CCCCCCCCC")
        self.assertEqual(
            self.sub_seq_extractor.extract(20, downstream=9),
            "GGGGGGGGG")

    def test_extract_down_and_upstream_seq(self):
        self.assertEqual(self.sub_seq_extractor.extract(
            33, upstream=3, downstream=3), "TTTATTT")
        self.assertEqual(self.sub_seq_extractor.extract(
                9, upstream=5, downstream=5), "AAAAAACCCCC")

    def test_extract_single_pos_revese(self):
         self.assertEqual(
             self.sub_seq_extractor.extract(9, rev_strand=True), "")
         self.assertEqual(
             self.sub_seq_extractor.extract(10, rev_strand=True), "")
         self.assertEqual(
             self.sub_seq_extractor.extract(19, rev_strand=True), "")
         self.assertEqual(
             self.sub_seq_extractor.extract(20, rev_strand=True), "")

    def test_extract_upstream_seq_reverse(self):
        self.assertEqual(
            self.sub_seq_extractor.extract(
                10, upstream=9, rev_strand=True),
                "GGGGGGGGG")
        self.assertEqual(
            self.sub_seq_extractor.extract(
                8, upstream=9, rev_strand=True),
                "GGGGGGGGT")

    def test_extract_downstream_seq_reverse(self):
        self.assertEqual(
            self.sub_seq_extractor.extract(
                19, downstream=9, rev_strand=True), "GGGGGGGGG")
        self.assertEqual(
            self.sub_seq_extractor.extract(
                12, downstream=4, rev_strand=True), "GGTT")

    def test_extract_down_and_upstream_seq_reverse(self):
        self.assertEqual(self.sub_seq_extractor.extract(
            33, upstream=3, downstream=3, rev_strand=True),
            "AAATAAA")
        self.assertEqual(self.sub_seq_extractor.extract(
            9, upstream=5, downstream=5, rev_strand=True),
            "GGGGGTTTTTT")

if __name__ == "__main__":
    unittest.main()
