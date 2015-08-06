""" This module contains graph format conversion functions """
import gmpy as gmp
import networkx as nx
from numpy import unravel_index



##++++ CONVERTED ++++##





def g2num(g):
    """ Convert a graph into a binary format """
    n = len(g)
    n2 = n ** 2 + n
    num = 0
    for v in range(1, n + 1):
        for w in g[v]:
            num |= (1 << (n2 - v * n - w))
    return num


def ug2num(g):
    """ Convert non-empty edges into a tuple of (directed, bidriected) in binary format """
    n = len(g)
    n2 = n ** 2 + n
    num = 0
    mask = 0
    num2 = 0
    for v in g:
        for w in g[v]:
            if g[v][w] in (1, 3):
                mask = (1 << (n2 - v * n - w))
                num |= mask
            if g[v][w] in (2, 3):
                num2 |= mask
    return num, num2


def bg2num(g):
    """ Convert bidirected edges into a binary format """
    n = len(g)
    n2 = n ** 2 + n
    num = 0
    for v in g:
        for w in g[v]:
            if g[v][w] in (2, 3):
                num = num | (1 << (n2 - v * n - w))
    return num


def graph2nx(G):
    g = nx.DiGraph()
    for v in G:
        g.add_edges_from([(v, x) for x in G[v] if G[v][x] in (1, 3)])
    return g


def nx2graph(G):
    g = {n: {} for n in G}
    for n in G:
        g[n] = {x: 1 for x in G[n]}
    return g


def num2CG(num, n):
    """num2CG - converts a number  whose binary representaion encodes edge
    presence/absence into a compressed graph representaion

    """
    n2 = n * n
    G = {i + 1: {} for i in xrange(n)}
    if num == 0:
        return G
    bl = gmp.bit_length(num)
    idx = [n2 - i - 1 for i in xrange(bl) if num & (1 << i)]
    idx = unravel_index(idx, (n, n))
    x = idx[0] + 1
    y = idx[1] + 1
    for i in xrange(len(x)):
        G[x[i]][y[i]] = 1
    return G


def dict_format_converter(H):
    """ Convert a graph from the set style dictionary format to the integer style

        >>> test = {'1': {'1': {(0, 1)},
        ...   '2': {(0, 1), (2, 0)},
        ...   '3': {(0, 1), (2, 0)},
        ...   '4': {(2, 0)},
        ...   '5': {(0, 1)}},
        ...  '2': {'1': {(2, 0)}, '2': {(0, 1)}, '5': {(0, 1), (2, 0)}},
        ...  '3': {'1': {(0, 1), (2, 0)}, '2': {(0, 1)}, '5': {(0, 1)}},
        ...  '4': {'1': {(2, 0)},
        ...   '2': {(0, 1)},
        ...   '3': {(0, 1)},
        ...   '4': {(0, 1)},
        ...   '5': {(0, 1)}},
        ...  '5': {'1': {(0, 1)}, '2': {(0, 1), (2, 0)}, '5': {(0, 1)}}}
        >>> dict_format_converter(test)
        {1: {1: 1, 2: 3, 3: 3, 4: 2, 5: 1}, 2: {1: 2, 2: 1, 5: 3}, 3: {1: 3, 2: 1, 5: 1}, 4: {1: 2, 2: 1, 3: 1, 4: 1, 5: 1}, 5: {1: 1, 2: 3, 5: 1}}
        >>>
    """
    H_new = {}
    for vert_a in H:
        H_new[int(vert_a)] = {}
        for vert_b in H[vert_a]:
            edge_val = 0
            if (0, 1) in H[vert_a][vert_b]:
                edge_val = 1
            if (2, 0) in H[vert_a][vert_b]:
                edge_val = 2 if edge_val == 0 else 3
            if edge_val:
                H_new[int(vert_a)][int(vert_b)] = edge_val
    return H_new
