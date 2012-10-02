import unittest
import sys
from io import StringIO
import kufpybio.helpers as helpers

class TestHelpers(unittest.TestCase):

    def test_overlap(self):
        self.assertEqual(helpers.overlap(1, 100, 50, 99), 50)
        self.assertEqual(helpers.overlap(50, 99, 1, 100), 50)
        self.assertEqual(helpers.overlap(0, 100, 100, 200), 1)
        self.assertEqual(helpers.overlap(0, 100, 101, 300), 0)
        self.assertEqual(helpers.overlap(0, 100, 200, 300), 0)

if __name__ == "__main__":
    unittest.main()
