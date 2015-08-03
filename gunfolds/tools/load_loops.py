import os
import sys
import zickle as zkl

DIR_NAME = os.path.dirname(__file__)
ABS_PATH = os.path.abspath(os.path.join(DIR_NAME))

alloops = zkl.load('{}/../data/allloops.zkl'.format(ABS_PATH))
circp = zkl.load('{}/../data/circular_p.zkl'.format(ABS_PATH))