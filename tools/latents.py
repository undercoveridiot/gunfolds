import sys, os
import sys, os

sys.path.append('./tools/')
import traversal as trv
import bfutils as bfu
import graphkit as gk
import numpy as np
import ecj
from copy import deepcopy
import ipdb


def g2lg(g):
    """
    Convert a data structure encoding the MSL-type graph into a structure encoding latents graph
    :return: a graph with integer indices and sets as weights
    """
    edge_type = {(0, 1): 1, (2, 0): 2}
    edge_weight = {(0, 1): 1, (2, 0): 0}
    lg = {int(e): {int(c): {edge_type[w]: {edge_weight[w]}
                            for w in g[e][c]}
                   for c in g[e]}
          for e in g}
    for e in lg:
        if e in lg[e]:
            lg[e][e] = {1: {LoopSet({1})}}
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
        self.loopset = lset
        if pre is None:
            self.preset = {[0]}
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
        s = '('
        comma = False
        if not self.preset == {0}:
            s += '{' + ', '.join(map(str, self.preset)) + '} '
            comma = True
        if not self.loopset == {0}:
            if comma:
                s += ', '
            s += '<{ ' + ', '.join(map(str, self.loopset)) + ' }>'
        s += ')'
        return s


def gtranspose(G):  # Transpose (rev. edges of) G
    GT = {u: {} for u in G}
    for u in G:
        for v in G[u]:
            if 1 in G[u][v]:
                GT[v][u] = {1: {1}}  # Add all reverse edges
    return GT


def parents(g, N): # an inefficient way to do it
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
        ws = osumset(ws, hh)
    if not type(ws) is set:
        ws = {ws}
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
                ab = gg[p][c][1]  # self loop
            else:
                ab = set()
            w = merge_weightsets(ab, pa[p], ch[c], sl)  # new weight set
            if c == p:  # a new loop is forming
                w = {LoopSet(w)}
            gg[p][c] = {1: w}

    return gg


def degrees(nodes, g):
    return [len(parents(g, v)) + len(children(g, v)) for v in nodes]


def sortbydegree(nodes, g):
    idx = np.argsort(degrees(nodes, g))
    return list(np.asarray(nodes)[idx])


def hide_nodes(g, nodelist, dosort=True):
    nodeset = set()  # make sure not to delete a node twice
    if dosort: nodelist = sortbydegree(nodelist, g)
    gg = deepcopy(g)
    for n in nodelist:
        if n in nodeset: continue
        gg = hide_node(gg, n)
        nodeset.add(n)
    return gg


def print_ws(ws):
    print '{',
    for e in ws:
        print e, ', ',
    print '}'


def test_osumnum():
    assert osumnum(set(range(5)), 1) == set(range(1, 5 + 1))


def testcase(n):
    g1 = {1: {2: {1: {1}}, 4: {1: {1}}},
          2: {3: {1: {1}}, 7: {1: {1}}},
          3: {},
          4: {5: {1: {1}}},
          5: {3: {1: {1}}, 6: {1: {1}}},
          6: {5: {1: {1}}},
          7: {8: {1: {1}}},
          8: {2: {1: {1}}}}

    g2 = {1: {2: {1: {1}}},
          2: {3: {1: {1}}, 6: {1: {1}}},
          3: {4: {1: {1}}},
          4: {5: {1: {1}}},
          5: {},
          6: {7: {1: {1}}},
          7: {3: {1: {1}}}}

    g3 = {1: {2: {1: {1}}},
          2: {3: {1: {1}}, 6: {1: {1}}},
          3: {4: {1: {1}}, 8: {1: {1}}},
          4: {5: {1: {1}}},
          5: {},
          6: {7: {1: {1}}},
          7: {2: {1: {1}}},
          8: {9: {1: {1}}},
          9: {10: {1: {1}}},
          10: {11: {1: {1}}},
          11: {3: {1: {1}}}}

    g4 = {1: {2: {1: {1}}},
          2: {3: {1: {1}}, 4: {1: {1}}},
          4: {5: {1: {1}}},
          5: {6: {1: {1}}, 7: {1: {1}}},
          6: {2: {1: {1}}},
          7: {8: {1: {1}}},
          8: {5: {1: {1}}},
          3: {}}

    g5 = {1: {2: {1: {1}}},
          2: {3: {1: {1}}, 4: {1: {1}}},
          4: {5: {1: {1}}},
          5: {6: {1: {1}}, 7: {1: {1}}},
          6: {2: {1: {1}}},
          7: {8: {1: {1}}},
          8: {5: {1: {1}}},
          3: {9: {1: {1}}},
          9: {10: {1: {1}}, 11: {1: {1}}},
          11: {9: {1: {1}}},
          10: {}}

    cases = [g1, g2, g3, g4, g5]

    return cases[n]
