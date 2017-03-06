# tools to construct (random) graphs
from gunfolds.tools.conversions import nx2graph
from gunfolds.tools import ecj
import igraph
import networkx as nx
import numpy as np
from numpy.random import randint
import random as std_random
import scipy


def edgelist(g):  # directed
    '''
    return a list of tuples for edges of g
    '''
    l = []
    for n in g:
        l.extend([(n, e) for e in g[n] if (0, 1) in g[n][e]])
    return l


def inedgelist(g):  # missing directed iterator
    '''
    iterate over the list of tuples for edges of g
    '''
    n = len(g)
    for v in g:
        for i in xrange(1, n + 1):
            w = str(i)
            if not w in g[v]:
                yield (v, w)
            elif not (0, 1) in g[v][w]:
                yield (v, w)


def inbedgelist(g):  # missing bidirected iterator
    '''
    iterate over the list of tuples for edges of g
    '''
    for v in g:
        for w in g:
            if v != w:
                if not w in g[v]:
                    yield (v, w)
                elif not (2, 0) in g[v][w]:
                    yield (v, w)


def bedgelist(g):  # bidirected edge list with flips
    l = []
    for n in g:
        l.extend([tuple(sorted((n, e))) for e in g[n] if (2, 0) in g[n][e]])
    l = list(set(l))
    l = l + map(lambda x: (x[1], x[0]), l)
    return l


def CG2adj(G):
    n = len(G)
    A = [[0 for i in range(0, n)] for j in range(0, n)]
    for v in G:
        if G[v]:
            directed = [w for w in G[v] if (0, 1) in G[v][w]]
            for w in directed:
                A[int(w) - 1][int(v) - 1] = 1
    A = np.double(np.asarray(A))
    return A


def g2ig(g):
    """
    Converts our graph represenataion to an igraph for plotting
    """
    t = scipy.where(CG2adj(g) == 1)
    l = zip(t[0], t[1])
    ig = igraph.Graph(l, directed=True)
    ig.vs["name"] = scipy.sort([u for u in g])
    ig.vs["label"] = ig.vs["name"]
    return ig


def superclique(n):
    g = {}
    for i in range(n):
        g[str(i + 1)] = {str(j + 1): set([(0, 1), (2, 0)])
                         for j in range(n) if j != i}
        g[str(i + 1)][str(i + 1)] = set([(0, 1)])
    return g


def complement(g):
    n = len(g)
    sq = superclique(n)
    for v in g:
        for w in g[v]:
            sq[v][w].difference_update(g[v][w])
            if not sq[v][w]:
                sq[v].pop(w)
    return sq


def gtranspose(G):                      # Transpose (rev. edges of) G
    GT = {u: {} for u in G}
    for u in G:
        for v in G[u]:
            if (0, 1) in G[u][v]:
                GT[v][u] = set([(0, 1)])        # Add all reverse edges
    return GT


def scale_free(n, alpha=0.7, beta=0.25,
               delta_in=0.2, delta_out=0.2):
    g = nx.scale_free_graph(n, alpha=alpha,
                            beta=beta,
                            delta_in=delta_in, delta_out=delta_out)
    g = nx2graph(g)
    g = gtranspose(g)
    addAring(g)
    return g


def ring(n, permute=False):
    g = {}
    names = [str(x+1) for x in range(n)]
    if permute: names = np.random.permutation(names) 
    for i in range(n-1):
        g[names[i]] = {names[i+1]: set([(0,1)])}
    g[names[n-1]] = {names[0]: set([(0,1)])}
    return g


def addAring(g):
    for i in range(1, len(g)):
        if str(i + 1) in g[str(i)]:
            g[str(i)][str(i + 1)].add((0, 1))
        else:
            g[str(i)][str(i + 1)] = set([(0, 1)])
    if '1' in g[str(len(g))]:
        g[str(len(g))]['1'].add((0, 1))
    else:
        g[str(len(g))]['1'] = set([(0, 1)])


def upairs(n, k):
    '''
    n unique nonsequential pairs
    '''
    s = set()
    for p in randint(n, size=(3 * k, 2)):
        if p[1] - p[0] == 1:
            continue
        s.add(tuple(p))
    l = [e for e in s]
    return l[:k]


def ringarcs(g, n):
    for edge in upairs(len(g), n):
        g[str(edge[0] + 1)][str(edge[1] + 1)] = set([(0, 1)])
    return g


def ringmore(n, m, permute=False):
    return ringarcs(ring(n,permute=permute), m)


def digonly(H):
    """returns a subgraph of H contatining all directed edges of H

    Arguments:
    - `H`: undersampled graph
    """
    g = {n: {} for n in H}
    for v in g:
        g[v] = {w: set([(0, 1)]) for w in H[v] if not H[v][w] == set([(2, 0)])}
    return g


def OCE(g1, g2):
    '''
    omission/commision error of g1 referenced to g2
    '''
    s1 = set(edgelist(g1))
    s2 = set(edgelist(g2))
    omitted = len(s2 - s1)
    comitted = len(s1 - s2)

    s1 = set(bedgelist(g1))
    s2 = set(bedgelist(g2))
    bomitted = len(s2 - s1)
    bcomitted = len(s1 - s2)

    return {'directed': (omitted, comitted),
            'bidirected': (bomitted, bcomitted)}


def clean_leaf_nodes(g):
    for v in g:
        g[v] = {w: g[v][w] for w in g[v] if g[v][w]}


def cerror(d):
    return d['OCE']['directed'][1] / np.double(len(d['gt']['graph']) ** 2 - len(edgelist(d['gt']['graph'])))


def oerror(d):
    return d['OCE']['directed'][0] / np.double(len(edgelist(d['gt']['graph'])))


def bidirected_no_fork(g):
    be = bedgelist(g)
    T = gtranspose(g)
    for e in be:
        if not set(T[e[0]].keys()) & set(T[e[1]].keys()):
            return True
    return False


def no_parents(g):
    T = gtranspose(g)
    for n in T:
        if not T[n]:
            return True
    return False


def no_children(g):
    for n in g:
        if not g[n]:
            return True
    return False


def scc_unreachable(g):
    if bidirected_no_fork(g):
        return True
    if no_parents(g):
        return True
    if no_children(g):
        return True
    return False

# unlike functions from traversal package these do no checking


def addanedge(g, e):
    g[e[0]][e[1]] = set([(0, 1)])


def delanedge(g, e):
    g[e[0]].pop(e[1], None)


def addedges(g, es):
    for e in es:
        addanedge(g, e)


def deledges(g, es):
    for e in es:
        delanedge(g, e)


def isdedgesubset(g2star, g2):
    '''
    check if g2star directed edges are a subset of those of g2
    '''
    for n in g2star:
        for h in g2star[n]:
            if h in g2[n]:
                if not (0, 1) in g2[n][h]:
                    return False
            else:
                return False
    return True


def isedgesubset(g2star, g2):
    '''
    check if all g2star edges are a subset of those of g2
    '''
    for n in g2star:
        for h in g2star[n]:
            if h in g2[n]:
                if not g2star[n][h].issubset(g2[n][h]):
                    return False
            else:
                return False
    return True
