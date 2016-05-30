import sys

sys.path.append('./tools/')

import bclique as bq
import latents as lt
import pathtreetools as ptt
from pathtree import PathTree


def can_add_loop(pt, num, elements):
    r = False
    if ptt.isptelement(PathTree(num - pt.preset, pre=pt.preset), num):
        r = True
    for e in pt.loopset:
        if type(e) is int:
            can_add_loop(PathTree(pt.preset))

def learn_path_tree(pt):
    elements = ptt.pt2seq(pt, 1)
    newpt = PathTree(set(), pre={elements[0]})

    def rpath(elements, npt):
        if not ptt.isptelement(npt, element[0]):


    return newpt

def bpts(bc):
    """
    Given a b-clique returns a set of PathTrees (one for each edge in bclique) such that loops are identified across edges in bclique
    :param bc: b-clique
    :return:
    """

def getagenerator(g):
    bclqs = bq.bcliques(g)
    for clq in bclqs:
        bpts(clq)
    return bclqs
