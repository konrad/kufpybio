import unittest
import sys
from kufpybio.tss import TSS

class TestTSS(unittest.TestCase):

    def test_init_without_extra(self):
        tss = TSS("plasmidX", "20", "-")
        self.assertEqual(tss.seq_id, "plasmidX")
        self.assertEqual(tss.pos, 20)
        self.assertEqual(tss.strand, "-")
        self.assertFalse(hasattr(tss, "extra"))

    def test_init_with_extra(self):
        tss = TSS("plasmidY", "30", "+", extra="weak")
        self.assertEqual(tss.seq_id, "plasmidY")
        self.assertEqual(tss.pos, 30)
        self.assertEqual(tss.strand, "+")
        self.assertEqual(tss.extra, "weak")

if __name__ == "__main__":
    unittest.main()
