import unittest
from your_project.api import Example


class TestExample(unittest.TestCase):

    def test___init___simple(self):
        ex = Example()
        self.assertEqual(ex.x, 1)