import sys, os
import sys, os

sys.path.append('./tools/')
import traversal as trv
import bfutils as bfu
import graphkit as gk
import numpy as np
from copy import deepcopy
import ipdb

def g2lg(g):
    """
    Convert a data structure encoding the MSL-type graph into a structure encoding latents graph
    :return: a graph with integer indices and sets as weights
    """
    edge_type = {(0, 1): 1, (2, 0): 2}
    edge_weight = {(0, 1): 1, (2, 0): 0}
    lg = {int(e): {int(c): {edge_type[w]:{edge_weight[w]}
                            for w in g[e][c]}
                   for c in g[e]}
          for e in g}
    return lg


def osumnum(s, num):
    return set(num + x for x in s)


def osumset(s1, s2):
    s = set()
    for e in s1:
        s.update(osumnum(s2, e))
    return s


class LoopSet:
    def __init__(self, lset, pre=None):
        assert (type(lset) == set)
        self.loopset = lset
        if pre is None:
            self.preset = set([0])
        else:
            assert (type(pre) == set)
            self.preset = pre

    def __add__(self, other):
        if type(other) is int:
            return LoopSet(self.loopset, pre=osumnum(self.preset, other))
        if type(other) is set:
            return LoopSet(self.loopset, pre=osumset(self.preset, other))
        return LoopSet(self.loopset.union(other.loopset), pre=osumset(self.preset, other.preset))

    def __radd__(self, other):
        return self.__add__(other)

    def __str__(self):
        s = ''
        if not self.preset == set([0]):
            s += '{' + ', '.join(map(str, self.preset)) + '} '
        if not self.loopset == set([0]):
            s += '<{ ' + ', '.join(map(str, self.loopset)) + ' }>'
        return s


def gtranspose(G):  # Transpose (rev. edges of) G
    GT = {u: {} for u in G}
    for u in G:
        for v in G[u]:
            if 1 in G[u][v]:
                GT[v][u] = {1: set([1])}  # Add all reverse edges
    return GT


def parents(g, N):
    t = gtranspose(g)  # Transpose the graph
    return {n: g[n][N][1]
            for n in t[N] if n != N}  # Return all children


def children(g, N):
    return {n: g[N][n][1]
            for n in g[N] if n != N}


def remove_node(g, N):
    del g[N]
    for V in g:
        if N in g[V]:
            del g[V][N]


def merge_weightsets(ab, ah, hb, hh):
    ws = osumset(ah, hb)
    if hh:
        ws = ws + hh
    if not type(ws) is set:
        ws = set([ws])
    return ws.union(ab)


def hide_node(g, H):
    """
    Removes a node from a graph taking care of the edges and weights as if the node was not observed
    :param g: input graph
    :param H: node to hide
    :return: the graph
    """

    gg = deepcopy(g)

    if not H in g:
        raise KeyError
    ch = children(g, H)
    pa = parents(g, H)
    if H in g[H]:
        sl = g[H][H][1]
    else:
        sl = set()
    remove_node(gg, H)

    for p in pa:
        for c in ch:
            if c in gg[p]:
                ab = gg[p][c]
            else:
                ab = set()
            gg[p][c] = {1: merge_weightsets(ab, pa[p], ch[c], sl)}

    return gg


def test_osumnum():
    assert osumnum(set(range(5)), 1) == set(range(1, 5 + 1))
