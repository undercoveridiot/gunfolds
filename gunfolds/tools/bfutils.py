from gunfolds.tools.conversions import g2num, ug2num, num2CG
import gunfolds.tools.zickle as zkl
import gunfolds.tools.graphkit as gk
import itertools
from numpy.random import randint
import numpy as np
import scipy
from itertools import permutations


#### Undersampling functions

def isSclique(G):
    n = len(G)
    for v in G:
        if sum([(0, 1) in G[v][w] for w in G[v]]) < n:
            return False
        if sum([(2, 0) in G[v][w] for w in G[v]]) < n - 1:
            return False
    return True


def directed_inc(G, D):
    G_un = {}
    # directed edges
    for v in D:
        G_un[v] = {}
        for w in D[v]:
            if G[w] and (0, 1) in D[v][w]:
                for e in G[w]:
                    G_un[v][e] = set([(0, 1)])
    return G_un


def bidirected_inc(G, D):
    # bidirected edges
    for w in G:
        # transfer old bidirected edges
        for l in D[w]:
            if (2, 0) in D[w][l]:
                G[w].setdefault(l, set()).add((2, 0))
        # new bidirected edges
        l = [e for e in D[w] if (0, 1) in D[w][e]]
        for pair in permutations(l, 2):
            G[pair[0]].setdefault(pair[1], set()).add((2, 0))
    return G


def increment_u(G_star, G_u):
    # directed edges
    G_un = directed_inc(G_star, G_u)
    # bidirected edges
    G_un = bidirected_inc(G_un, G_u)
    return G_un


def pure_directed_inc(G, D):
    G_un = {}
    # directed edges
    for v in D:
        G_un[v] = {}
        for w in D[v]:
            if G[w]:
                for e in G[w]:
                    G_un[v][e] = set([(0, 1)])
    return G_un


def increment(g):
    '''
    undersample g by 2
    only works for g1 to g2 directed
    '''
    r = {n: {} for n in g}

    for n in g:
        for h in g[n]:
            for e in g[h]:
                if not e in r[n]:
                    r[n][e] = set([(0, 1)])

    for n in g:
        for pair in itertools.combinations(g[n], 2):

            if pair[1] in r[pair[0]]:
                r[pair[0]][pair[1]].add((2, 0))
            else:
                r[pair[0]][pair[1]] = set([(2, 0)])

            if pair[0] in r[pair[1]]:
                r[pair[1]][pair[0]].add((2, 0))
            else:
                r[pair[1]][pair[0]] = set([(2, 0)])

    return r


def dincrement_u(G_star, G_u):
    # directed edges
    G_un = pure_directed_inc(G_star, G_u)
    return G_un


def undersample(G, u):
    Gu = G
    for i in range(u):
        Gu = increment_u(G, Gu)
    return Gu


def all_undersamples(G_star):
    glist = [G_star]
    while True:
        g = increment_u(G_star, glist[-1])
        if isSclique(g):
            return glist  # superclique convergence
        # this will (may be) capture DAGs and oscillations
        if g in glist:
            return glist
        glist.append(g)
    return glist


def call_undersamples(G_star):
    glist = [G_star]
    while True:
        g = increment_u(G_star, glist[-1])
        if g in glist:
            return glist
        glist.append(g)
    return glist


def compact_call_undersamples(G_star, steps=None):
    glist = [ug2num(G_star)]
    lastgraph = G_star
    while True:
        g = increment_u(G_star, lastgraph)
        if ug2num(g) in glist:
            return glist
        glist.append(ug2num(g))
        lastgraph = g
    return glist


def cc_undersamples(G_star, steps=1):
    glist = [ug2num(G_star)]
    lastgraph = G_star
    for i in xrange(steps):
        g = increment_u(G_star, lastgraph)
        n = ug2num(g)
        if n in glist:
            return []
        glist.append(n)
        lastgraph = g
    return glist[-1]




#### Adjacency matrix functions

def graph2adj(G):
    n = len(G)
    A = scipy.zeros((n, n), dtype=np.int8)
    for v in G:
        A[int(v) - 1, [int(w)-1 for w in G[v] if (0, 1) in G[v][w]]] = 1
    return A


def graph2badj(G):
    n = len(G)
    A = scipy.zeros((n, n), dtype=np.int8)
    for v in G:
        A[int(v) - 1, [int(w)-1 for w in G[v] if (2, 0) in G[v][w]]] = 1
    return A


def adjs2graph(A, B):
    names = [str(i) for i in range(1, A.shape[0] + 1)]
    G = {}
    for name in names:
        G[name] = {}
    for i in range(A.shape[0]):
        for name in map(str, np.where(A[i,:] == 1)[0]+1):
            G[str(i + 1)][name] = set([(0, 1)])

    for i in range(B.shape[0]):
        for j in range(B.shape[1]):
            if B[i, j]:
                if str(j + 1) in G[str(i+1)]:
                    G[str(i + 1)][str(j+1)].add((2, 0))
                else:
                    G[str(i + 1)][str(j+1)] = set([(2, 0)])
    return G


def g2vec(g):
    A = graph2adj(g)
    B = graph2badj(g)
    return np.r_[A.flatten(), B[np.triu_indices(B.shape[0])]]


def vec2adj(v, n):
    A = np.zeros((n, n))
    B = np.zeros((n, n))
    A[:] = v[:n ** 2].reshape(n, n)
    B[np.triu_indices(n)] = v[n ** 2:]
    B = B + B.T
    return A, B


def vec2g(v, n):
    A, B = vec2adj(v, n)
    return adjs2graph(A, B)




#### Misc graph functions

def overshoot(G_star, H):
    glist = [G_star]
    while True:
        g = increment_u(G_star, glist[-1])
        if isSclique(g):
            return False
        if gk.isedgesubset(H, g):
            return True
        if g in glist:
            return False
        glist.append(g)
    return False


def forms_loop(G_star, loop):
    glist = [G_star]
    while True:
        g = increment_u(G_star, glist[-1])
        if (g2num(gk.digonly(g)) & loop) == loop:
            return True
        if g in glist:
            return False
        glist.append(g)
    return False


def call_u_conflicts_d(G_star, H, checkrate=0):
    glist = [G_star]
    while True:
        g = dincrement_u(G_star, glist[-1])
        if gk.isedgesubset(g, H):
            return False
        if g in glist:
            return True
        glist.append(g)
    return True


def call_u_conflicts(G_star, H, checkrate=0):
    glist = [G_star]
    while True:
        # g = increment_u(G_star, glist[-1])
        g = directed_inc(G_star, glist[-1])
        if gk.isdedgesubset(g, H):
            return False
        g = bidirected_inc(g, glist[-1])
        if gk.isedgesubset(g, H):
            return False
        if g in glist:
            return True
        glist.append(g)
    return True


def call_u_conflicts2(G_star, H):
    glist = [G_star]
    while True:
        g = increment_u(G_star, glist[-1])
        if gk.isedgesubset(g, H):
            return False, glist
        if g in glist:
            return True, glist
        glist.append(g)
    return True, glist


def call_u_equals2(G_star, glist, H):
    while True:
        g = increment_u(G_star, glist[-1])
        if g == H:
            return True
        if g in glist:
            return False
        glist.append(g)
    return False


def call_u_equals(G_star, H):
    glist = [G_star]
    while True:
        g = increment_u(G_star, glist[-1])
        if g == H:
            return True
        if g in glist:
            return False
        glist.append(g)
    return False


def compatible(d1, d2):
    idx = scipy.where(scipy.array([[r == l for l in d2] for r in d1]))
    return idx


def compat(G):
    n = len(G)
    num = g2num(G)
    # sample all the graph for gStar
    star_l = all_undersamples(G)
    hits = {}
    # brute force all graphs
    for i in range(0, 2 ** (n**2)):
        tmpG = num2CG(i, n)
        tmp_l = all_undersamples(tmpG)
        c = compatible(tmp_l, star_l)
        if len(sum(c)) > 0:
            hits[i] = c
    return hits


def icompat(i, nodes):
    print i
    g = num2CG(i, nodes)
    return compat(g)


def ilength(i, nodes):
    print i
    g = num2CG(i, nodes)
    return len(call_undersamples(g))


def iall(i, nodes):
    print i
    g = num2CG(i, nodes)
    return compact_call_undersamples(g)


def cc_all(i, nodes, steps):
    # print i
    g = num2CG(i, nodes)
    return cc_undersamples(g, steps=steps)


def make_rect(l):
    max_seq = max(map(len, l))
    nl = []
    for e in l:
        e += [e[-1]] * (max_seq - len(e))
        nl.append(e)
    return nl


def loadgraphs(fname):
    g = zkl.load(fname)
    return g


def savegraphs(l, fname):
    zkl.save(l, fname)


# talking about extra edges on top of the ring


def dens2edgenum(d, n=10):
    return int(d * n**2)-n


def edgenum2dens(e, n=10):
    return np.double(e + n)/n**2


def randH(n, d1, d2):
    g = gk.ringmore(n, d1)
    pairs = [x for x in itertools.combinations(g.keys(), 2)]
    for p in np.random.permutation(pairs)[:d2]:
        g[p[0]].setdefault(p[1], set()).add((2, 0))
        g[p[1]].setdefault(p[0], set()).add((2, 0))
    return g


def checkconflict(H, G_test, au=None):
    if not au:
        allundersamples = call_undersamples(G_test)
    else:
        allundersamples = au
    for graph in allundersamples:
        if isedgesubset(graph, H):
            return False
    return True


def checkconflict_(Hnum, G_test, au=None):
    if not au:
        allundersamples = call_undersamples(G_test)
    else:
        allundersamples = au
    # Hnum = ug2num(H)
    for graph in allundersamples:
        gnum = ug2num(graph)
        if gnum[0] & Hnum[0] == gnum[0] and gnum[1] & Hnum[1] == gnum[1]:
            return False
    return True


def checkequality(H, G_test, au=None):
    if not au:
        allundersamples = call_undersamples(G_test)
    else:
        allundersamples = au
    for graph in allundersamples:
        if graph == H:
            return True
    return False
