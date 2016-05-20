import sys

sys.path.append('./tools/')

import bclique as bq
import latents as lt
import pathtreetools as ptt

def bpts(bc):
    """
    Given a b-clique returns a set of PathTrees (one for each edge in bclique) such that loops are identified across edges in bclique
    :param bc: b-clique
    :return:
    """

def getagenerator(g):
    bclqs = bq.bcliques(g)

