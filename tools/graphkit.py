# tools to construct (random) graphs
import sys
sys.path.append('/home/splis/soft/src/dev/craft/gunfolds/tools/')

import ecj, bfutils
import random as std_random
import scipy
import networkx as nx
import igraph

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
    G = ecj.adj2graph(l)
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
    return ecj.adj2graph(np.asarray(A))

def fullG(n):
    A = [[1 for j in range(n)] for i in range(n)]
    return ecj.adj2graph(np.asarray(A))


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
