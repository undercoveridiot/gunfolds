# tools to construct (random) graphs
import sys,os
TOOLSPATH='~/soft/src/dev/craft/gunfolds/tools/'
sys.path.append(os.path.expanduser(TOOLSPATH))

import ecj
import bfutils as bfu
import traversal as trv

import random as std_random
import numpy as np
import scipy
import networkx as nx
import igraph
from numpy.random import randint

def edgelist(g): # directed
    '''
    return a list of tuples for edges of g
    '''
    l = []
    for n in g:
        l.extend([(n,e) for e in g[n] if (0,1) in g[n][e]])
    return l

def edgenumber(g):
    return sum([sum([len(g[y][x]) for x in g[y]]) for y in g])

def bedgelist(g): # bidirected edge list with flips
    l = []
    for n in g:
        l.extend([tuple(sorted((n,e))) for e in g[n] if (2,0) in g[n][e]])
    l = list(set(l))
    l = l + map(lambda x: (x[1],x[0]), l)
    return l

def rnd_edges(n):
    """ generate a random uniformly distributed mask
    """
    rnum = std_random.getrandbits(n**2)
    l = list(bin(rnum)[2:])
    l = ['0' for i in range(0,n**2 - len(l))] + l
    return l

def list2dbn(l):
    """ convert list of edge presences/absences (0,1) to a DBN graph
    """
    n = scipy.sqrt(len(l))
    l = scipy.reshape(map(int, l), [n, n])
    G = ecj.adj2DBN(l)
    return G

def list2CG(l):
    """ convert list of edge presences/absences (0,1) to a compressed
    graph (CG) representation
    """
    n = scipy.sqrt(len(l))
    l = scipy.reshape(map(int, l), [n, n])
    G = bfu.adj2graph(l)
    return G

def rnd_dbn(n): return list2dbn(rnd_edges(n))
def rnd_cg(n):  return list2CG(rnd_edges(n))

def rnd_adj(n, maxindegree=5):
    l = scipy.zeros([n,n])
    for u in range(0,n):
        cap = scipy.random.randint(min([n,maxindegree+1]))
        idx = scipy.random.randint(n,size=cap)
        l[u, idx] = 1
    return l

def sp_rnd_edges(n, maxindegree=5):
    '''
    a sparse set of random edges
    '''
    l = rnd_adj(n, maxindegree=maxdegree)
    return scipy.reshape(l, n**2)

def sp_rnd_dbn(n, maxindegree=3):
    '''
    a sparse random DBN graph
    '''
    l = sp_rnd_edges(n)
    return list2dbn(l)

def emptyG(n):
    A = [[0 for j in range(n)] for i in range(n)]
    return bfu.adj2graph(np.asarray(A))

def fullG(n):
    A = [[1 for j in range(n)] for i in range(n)]
    return bfu.adj2graph(np.asarray(A))


def CG2uCG(cg):
    """
    convert to an undirected graph
    """
    G = {}
    for u in cg:
        G[u] = cg[u].copy()
    for u in cg:
        for v in cg[u]:
            G[v][u] = cg[u][v]
    return G

def connected(cg):
    n = len(cg)
    return sum(1 for _ in ecj.traverse(CG2uCG(cg),'1')) == n

def sp_rnd_CG(n, maxindegree=3, force_connected=False):
    l = sp_rnd_edges(n, maxindegree=maxindegree)
    cg = list2CG(l)
    if force_connected:
        while not connected(cg):
            cg = list2CG(sp_rnd_edges(n, maxindegree=maxindegree))
    return cg

def CG2adj(G):
    n = len(G)
    A = [[0 for i in range(0,n)] for j in range(0,n)]
    for v in G:
        if G[v]:
            directed = [w for w in G[v] if (0,1) in G[v][w]]
            for w in directed:
                A[int(w)-1][int(v)-1] = 1
    A = np.double(np.asarray(A))
    return A

def g2ig(g):
    """
    Converts our graph represenataion to an igraph for plotting
    """
    t = scipy.where(CG2adj(g)==1)
    l = zip(t[0],t[1])
    ig = igraph.Graph(l,directed=True)
    ig.vs["name"] = scipy.sort([u for u in g])
    ig.vs["label"] = ig.vs["name"]
    return ig

def superclique(n):
    g = {}
    for i in range(n):
        g[str(i+1)] = {str(j+1):set([(0,1),(2,0)])
                       for j in range(n) if j!=i}
        g[str(i+1)][str(i+1)] = set([(0,1)])
    return g

def complement(g):
    n = len(g)
    sq = superclique(n)
    for v in g:
        for w in g[v]:
            sq[v][w].difference_update(g[v][w])
            if not sq[v][w]: sq[v].pop(w)
    return sq

def gtranspose(G):                      # Transpose (rev. edges of) G
    GT = {u:{} for u in G}
    for u in G:
        for v in G[u]:
            GT[v][u] = set([(0,1)])        # Add all reverse edges
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

def ring(n):
    g = {}
    for i in range(1,n):
        g[str(i)] = {str(i+1): set([(0,1)])}
    g[str(n)] = {'1': set([(0,1)])}
    return g

def addAring(g):
    for i in range(1,len(g)):
        if str(i+1) in g[str(i)]:
            g[str(i)][str(i+1)].add((0,1))
        else:
            g[str(i)][str(i+1)] = set([(0,1)])
    if '1' in g[str(len(g))]:
        g[str(len(g))]['1'].add((0,1))
    else:
        g[str(len(g))]['1'] = set([(0,1)])

def upairs(n,k):
    '''
    n unique nonsequential pairs
    '''
    s = set()
    for p in randint(n, size=(3*k, 2)):
        if p[1]-p[0] == 1: continue
        s.add(tuple(p))
    l = [e for e in s]
    return l[:k]

def ringarcs(g,n):
    for edge in upairs(len(g),n):
        g[str(edge[0]+1)][str(edge[1]+1)] = set([(0,1)])
    return g
def ringmore(n,m):
    return ringarcs(ring(n),m)


# Justin's ternary representation: 1 = directed edge; 2 = bidirected; 3 = both
def justin2graph(g):
    r = {}
    d = {1: set([(0,1)]),
         2: set([(2,0)]),
         3: set([(0,1),(2,0)]) }
    for head in g:
        r[head] = {}
        for tail in g[head]:
            r[head][tail] = d[g[head][tail]]
    return r

def graph2justin(g):
    r = {}
    for head in g:
        r[head] = {}
        for tail in g[head]:
            if g[head][tail] == set([(0,1)]):
                r[head][tail] = 1
            elif g[head][tail] == set([(2,0)]):
                r[head][tail] = 2
            elif g[head][tail] == set([(0,1),(2,0)]):
                r[head][tail] = 3
    return r

def OCE(g1,g2):
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
    for v in g: g[v] = {w:g[v][w] for w in g[v] if g[v][w]}        

def cerror(d):
    return d['OCE']['directed'][1]/np.double(len(d['gt']['graph'])**2-len(edgelist(d['gt']['graph'])))

def oerror(d):
    return d['OCE']['directed'][0]/np.double(len(edgelist(d['gt']['graph'])))

def bidirected_no_fork(g):
    be = bedgelist(g)
    T = gtranspose(g)
    for e in be:
        if not set(T[e[0]].keys())&set(T[e[1]].keys()):
            return True
    return False

def no_parents(g):
    T = gtranspose(g)
    for n in T:
        if not T[n]: return True
    return False

def no_children(g):
    for n in g:
        if not g[n]: return True
    return False

def scc_unreachable(g):
    if bidirected_no_fork(g): return True
    if no_parents(g): return True
    if no_children(g): return True    
    return False 
