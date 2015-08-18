from gunfolds.tools import bfutils
from gunfolds.tools import conversions
from gunfolds.tools import unknownrate as ur
from gunfolds.tools import zickle as zkl
import os
import unittest


class TestBFUtilsFunctions(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls._G = {'1': {'1': set([(0, 1)]),
                   '2': set([(0, 1)]),
                   '3': set([(0, 1)])},
             '3': {'1': set([(0, 1)]),
                   '3': set([(0, 1)]),
                   '2': set([(0, 1)]),
                   '4': set([(0, 1)])},
             '2': {'3': set([(0, 1)]),
                   '2': set([(0, 1)]),
                   '4': set([(0, 1)])},
             '5': {'1': set([(0, 1)])},
             '4': {'1': set([(0, 1)]),
                   '3': set([(0, 1)]),
                   '2': set([(0, 1)]),
                   '5': set([(0, 1)])}}

    def test__call_undersamples(self):
        gs = [{'1': {'1': {(0, 1)}, '2': {(0, 1)}, '3': {(0, 1)}},
              '2': {'2': {(0, 1)}, '3': {(0, 1)}, '4': {(0, 1)}},
              '3': {'1': {(0, 1)}, '2': {(0, 1)}, '3': {(0, 1)}, '4': {(0, 1)}},
              '4': {'1': {(0, 1)}, '2': {(0, 1)}, '3': {(0, 1)}, '5': {(0, 1)}},
              '5': {'1': {(0, 1)}}},
              {'1': {'1': {(0, 1)},
               '2': {(0, 1), (2, 0)},
               '3': {(0, 1), (2, 0)},
               '4': {(0, 1), (2, 0)},
               '5': {(2, 0)}},
              '2': {'1': {(0, 1), (2, 0)},
               '2': {(0, 1)},
               '3': {(0, 1), (2, 0)},
               '4': {(0, 1), (2, 0)},
               '5': {(0, 1), (2, 0)}},
              '3': {'1': {(0, 1), (2, 0)},
               '2': {(0, 1), (2, 0)},
               '3': {(0, 1)},
               '4': {(0, 1), (2, 0)},
               '5': {(0, 1), (2, 0)}},
              '4': {'1': {(0, 1), (2, 0)},
               '2': {(0, 1), (2, 0)},
               '3': {(0, 1), (2, 0)},
               '4': {(0, 1)}},
              '5': {'1': {(0, 1), (2, 0)}, '2': {(0, 1), (2, 0)}, '3': {(0, 1), (2, 0)}}},
              {'1': {'1': {(0, 1)},
               '2': {(0, 1), (2, 0)},
               '3': {(0, 1), (2, 0)},
               '4': {(0, 1), (2, 0)},
               '5': {(0, 1), (2, 0)}},
              '2': {'1': {(0, 1), (2, 0)},
               '2': {(0, 1)},
               '3': {(0, 1), (2, 0)},
               '4': {(0, 1), (2, 0)},
               '5': {(0, 1), (2, 0)}},
              '3': {'1': {(0, 1), (2, 0)},
               '2': {(0, 1), (2, 0)},
               '3': {(0, 1)},
               '4': {(0, 1), (2, 0)},
               '5': {(0, 1), (2, 0)}},
              '4': {'1': {(0, 1), (2, 0)},
               '2': {(0, 1), (2, 0)},
               '3': {(0, 1), (2, 0)},
               '4': {(0, 1)},
               '5': {(0, 1), (2, 0)}},
              '5': {'1': {(0, 1), (2, 0)},
               '2': {(0, 1), (2, 0)},
               '3': {(0, 1), (2, 0)},
               '4': {(0, 1), (2, 0)}}},
              {'1': {'1': {(0, 1)},
               '2': {(0, 1), (2, 0)},
               '3': {(0, 1), (2, 0)},
               '4': {(0, 1), (2, 0)},
               '5': {(0, 1), (2, 0)}},
              '2': {'1': {(0, 1), (2, 0)},
               '2': {(0, 1)},
               '3': {(0, 1), (2, 0)},
               '4': {(0, 1), (2, 0)},
               '5': {(0, 1), (2, 0)}},
              '3': {'1': {(0, 1), (2, 0)},
               '2': {(0, 1), (2, 0)},
               '3': {(0, 1)},
               '4': {(0, 1), (2, 0)},
               '5': {(0, 1), (2, 0)}},
              '4': {'1': {(0, 1), (2, 0)},
               '2': {(0, 1), (2, 0)},
               '3': {(0, 1), (2, 0)},
               '4': {(0, 1)},
               '5': {(0, 1), (2, 0)}},
              '5': {'1': {(0, 1), (2, 0)},
               '2': {(0, 1), (2, 0)},
               '3': {(0, 1), (2, 0)},
               '4': {(0, 1), (2, 0)},
               '5': {(0, 1)}}}]

        gs_test = bfutils.call_undersamples(self._G)
        self.assertEqual(gs, gs_test)

    def test__call_undersample(self):
        u = 1
        g_u_1 = {'1': {'1': {(0, 1)},
                       '2': {(0, 1), (2, 0)},
                 '3': {(0, 1), (2, 0)},
                       '4': {(0, 1), (2, 0)},
                       '5': {(2, 0)}},
                 '2': {'1': {(0, 1), (2, 0)},
                       '2': {(0, 1)},
                       '3': {(0, 1), (2, 0)},
                       '4': {(0, 1), (2, 0)},
                       '5': {(0, 1), (2, 0)}},
                 '3': {'1': {(0, 1), (2, 0)},
                       '2': {(0, 1), (2, 0)},
                       '3': {(0, 1)},
                       '4': {(0, 1), (2, 0)},
                       '5': {(0, 1), (2, 0)}},
                 '4': {'1': {(0, 1), (2, 0)},
                       '2': {(0, 1), (2, 0)},
                       '3': {(0, 1), (2, 0)},
                       '4': {(0, 1)}},
                 '5': {'1': {(0, 1), (2, 0)}, '2': {(0, 1), (2, 0)}, '3': {(0, 1), (2, 0)}}}
        g2 = bfutils.undersample(self._G, u)
        self.assertEqual(g_u_1, g2)

        u = 2
        g_u_2 = {'1': {'1': {(0, 1)},
                       '2': {(0, 1), (2, 0)},
                 '3': {(0, 1), (2, 0)},
                       '4': {(0, 1), (2, 0)},
                       '5': {(0, 1), (2, 0)}},
                 '2': {'1': {(0, 1), (2, 0)},
                       '2': {(0, 1)},
                       '3': {(0, 1), (2, 0)},
                       '4': {(0, 1), (2, 0)},
                       '5': {(0, 1), (2, 0)}},
                 '3': {'1': {(0, 1), (2, 0)},
                       '2': {(0, 1), (2, 0)},
                       '3': {(0, 1)},
                       '4': {(0, 1), (2, 0)},
                       '5': {(0, 1), (2, 0)}},
                 '4': {'1': {(0, 1), (2, 0)},
                       '2': {(0, 1), (2, 0)},
                       '3': {(0, 1), (2, 0)},
                       '4': {(0, 1)},
                       '5': {(0, 1), (2, 0)}},
                 '5': {'1': {(0, 1), (2, 0)},
                       '2': {(0, 1), (2, 0)},
                       '3': {(0, 1), (2, 0)},
                       '4': {(0, 1), (2, 0)}}}
        g2 = bfutils.undersample(self._G, u)
        self.assertEqual(g_u_2, g2)

        u = 4
        g_u_4 = {'1': {'1': {(0, 1)},
                       '2': {(0, 1), (2, 0)},
                 '3': {(0, 1), (2, 0)},
                       '4': {(0, 1), (2, 0)},
                       '5': {(0, 1), (2, 0)}},
                 '2': {'1': {(0, 1), (2, 0)},
                       '2': {(0, 1)},
                       '3': {(0, 1), (2, 0)},
                       '4': {(0, 1), (2, 0)},
                       '5': {(0, 1), (2, 0)}},
                 '3': {'1': {(0, 1), (2, 0)},
                       '2': {(0, 1), (2, 0)},
                       '3': {(0, 1)},
                       '4': {(0, 1), (2, 0)},
                       '5': {(0, 1), (2, 0)}},
                 '4': {'1': {(0, 1), (2, 0)},
                       '2': {(0, 1), (2, 0)},
                       '3': {(0, 1), (2, 0)},
                       '4': {(0, 1)},
                       '5': {(0, 1), (2, 0)}},
                 '5': {'1': {(0, 1), (2, 0)},
                       '2': {(0, 1), (2, 0)},
                       '3': {(0, 1), (2, 0)},
                       '4': {(0, 1), (2, 0)},
                       '5': {(0, 1)}}}
        g2 = bfutils.undersample(self._G, u)
        self.assertEqual(g_u_4, g2)

class TestTraversalFunctions(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # Read in pickle file with results
        cls._G1 = {1: {1: 1, 2: 1, 3: 1},
                   2: {1: 1, 3: 1, 4: 1},
                   3: {1: 1, 2: 1, 3: 1, 4: 1},
                   4: {1: 1, 2: 1, 3: 1, 5: 1},
                   5: {2: 1}}
        cls._G2 = {1: {1: 1, 2: 3, 3: 3, 4: 3, 5: 2},
                   2: {1: 3, 2: 1, 3: 3, 4: 3, 5: 3},
                   3: {1: 3, 2: 3, 3: 1, 4: 3, 5: 3},
                   4: {1: 3, 2: 3, 3: 3, 4: 1},
                   5: {1: 3, 2: 2, 3: 3, 4: 1}}
        DIR_NAME = os.path.dirname(__file__)
        cls._ABS_PATH = os.path.abspath(os.path.join(DIR_NAME))

    def setUp(self):
        # Copy these before each test incase a test modifies them
        self.G1 = deepcopy(self._G1)
        self.G2 = deepcopy(self._G2)

    def test__v2g22g1(self):
        expected = zkl.load("{}/v2g22g1_output.zkl".format(self._ABS_PATH))
        self.assertEqual(expected, traversal.v2g22g1(self.G2))


class TestConversionFunctions(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls._G = {'1': {'1': set([(0, 1)]),
                   '2': set([(0, 1)]),
                   '3': set([(0, 1)])},
             '3': {'1': set([(0, 1)]),
                   '3': set([(0, 1)]),
                   '2': set([(0, 1)]),
                   '4': set([(0, 1)])},
             '2': {'3': set([(0, 1)]),
                   '2': set([(0, 1)]),
                   '4': set([(0, 1)])},
             '5': {'1': set([(0, 1)])},
             '4': {'1': set([(0, 1)]),
                   '3': set([(0, 1)]),
                   '2': set([(0, 1)]),
                   '5': set([(0, 1)])}}


    def test__dict_format_converter(self):
        expected = {1: {1: 1, 2: 1, 3: 1},
                    2: {2: 1, 3: 1, 4: 1},
                    3: {1: 1, 2: 1, 3: 1, 4: 1},
                    4: {1: 1, 2: 1, 3: 1, 5: 1},
                    5: {1: 1}}
        converted = conversions.dict_format_converter(self._G)
        self.assertEqual(expected, converted)




class TestUnknownRateFunctions(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # Read in pickle file with results
        cls._G = {'1': {'1': set([(0, 1)]),
                   '2': set([(0, 1)]),
                   '3': set([(0, 1)])},
             '3': {'1': set([(0, 1)]),
                   '3': set([(0, 1)]),
                   '2': set([(0, 1)]),
                   '4': set([(0, 1)])},
             '2': {'3': set([(0, 1)]),
                   '2': set([(0, 1)]),
                   '4': set([(0, 1)])},
             '5': {'1': set([(0, 1)])},
             '4': {'1': set([(0, 1)]),
                   '3': set([(0, 1)]),
                   '2': set([(0, 1)]),
                   '5': set([(0, 1)])}}
        DIR_NAME = os.path.dirname(__file__)
        cls._ABS_PATH = os.path.abspath(os.path.join(DIR_NAME))

    def test__liteqclass(self):
        u = 1
        set_of_u_1 = zkl.load("{}/liteqclass_output_n_5_u_1.zkl".format(self._ABS_PATH))
        g2 = bfutils.undersample(self._G, u)
        s = ur.liteqclass(g2, verbose=False, capsize=1000)
        self.assertEqual(s, set_of_u_1)

        u = 2
        set_of_u_2 = zkl.load("{}/liteqclass_output_n_5_u_2.zkl".format(self._ABS_PATH))
        g2 = bfutils.undersample(self._G, u)
        s = ur.liteqclass(g2, verbose=False, capsize=1000)
        self.assertEqual(s, set_of_u_2)


if __name__ == '__main__':
    try:
        from teamcity import is_running_under_teamcity
        from teamcity.unittestpy import TeamcityTestRunner

        if is_running_under_teamcity():
            runner = TeamcityTestRunner()
        else:
            runner = unittest.TextTestRunner()
    except ImportError:
        runner = unittest.TextTestRunner()

    unittest.main(testRunner=runner, verbosity=2)
