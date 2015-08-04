from gunfolds.tools.conversions import num2CG
from gunfolds.tools import ecj
import gmpy as gmp
import networkx as nx
import numpy as np
import operator
import scipy

np.random.RandomState()


def hasSelfLoops(G):
    for u in G:
        if G[u].has_key(u):
            return True
    return False


def randSCC(n):
    G = num2CG(scipy.random.randint(2 ** (n ** 2)), n)
    while (len(ecj.scc(G)) > 1) or gcd4scc(G) > 1:
        G = num2CG(scipy.random.randint(2 ** (n ** 2)), n)
    return G


def SM_fixed(Gstar, G, iter=5):
    compat = []
    for j in range(0, iter):
        if Gstar == ecj.undersample(G, j):
            compat.append(j)
    return compat


def SM_converging(Gstar, G):
    """Gstar is the undersampled reference graph, while G is the starting
    graph. The  code searches  over all undersampled  version of  G to
    find all matches with Gstar
    """
    compat = []
    GG = G
    Gprev = G
    if G == Gstar:
        return [0]
    j = 1
    G = ecj.undersample(GG, j)
    while not (G == Gprev):
        if Gstar == G:
            compat.append(j)
        j += 1
        Gprev = G
        G = ecj.undersample(GG, j)
    return compat


def searchMatch(Gstar, G, iter=5):
    if gcd4scc(G) > 1:
        return SM_fixed(Gstar, G, iter=iter)
    return SM_converging(Gstar, G)


def hasSink(G):
    return not reduce(operator.and_, [bool(G[n]) for n in G], True)


def hasRoot(G):
    return hasSink(ecj.tr(G))


def graph2nx(G):
    g = nx.DiGraph()
    for v in G:
        g.add_edges_from([(v, x) for x in G[v] if (0, 1) in G[v][x]])
    return g


def nx2graph(G):
    g = {str(n + 1): {} for n in G}
    for n in G:
        g['%i' % (n + 1)] = {'%i' % (x + 1): set([(0, 1)]) for x in G[n]}
    return g


def gcd4scc(SCC):
    g = graph2nx(SCC)
    return ecj.listgcd(map(lambda x: len(x) - 1, nx.simple_cycles(g)))


def compatibleAtU(uGstar):
    compat = []
    n = len(uGstar)
    numG = 2 ** (n ** 2)
    # pbar = Percentage()
    for i in range(1, numG):
        G = num2CG(i, n)
        # pbar.update(i+1)
        if len(ecj.scc(G)) > 1:
            continue
        l = searchMatch(uGstar, G, iter=5)
        if l:
            compat.append((l, G))
    # pbar.finish()
    return compat
