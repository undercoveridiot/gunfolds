import scipy
import itertools
#from progressbar import ProgressBar, Percentage
from multiprocessing import Pool,Array,Process,Manager
from numpy.random import randint
import numpy as np
#import ipdb
import networkx as nx
#local
import ecj
import zickle as zkl
from comparison import num2CG, nx2graph

def directed_inc(G,D):
    G_un = {}
    # directed edges
    for v in D:
        G_un[v] = {}
        for w in [el for el in D[v] if (0,1) in D[v][el]]:
            for e in G[w]:
                G_un[v][e] = set([(0,1)])
    return G_un
def bidirected_inc(G,D):
    # bidirected edges
    for w in G:
        # transfer old bidirected edges
        l = [e for e in D[w] if (2,0) in D[w][e]]
        for p in l:
            if p in G[w]: 
                G[w][p].add((2,0))
            else: 
                G[w][p] = set([(2,0)])
        # new bidirected edges
        l = [e for e in D[w] if (0,1) in D[w][e]]
        for pair in itertools.permutations(l,2):
            if pair[1] in G[pair[0]]:
                G[pair[0]][pair[1]].add((2,0))
            else:
                G[pair[0]][pair[1]] = set([(2,0)])
    return G
def increment(g):
    '''
    undersample g by 2
    only works for g1 to g2 directed
    '''
    r = {n:{} for n in g}

    for n in g:
        for h in g[n]:
            for e in g[h]:
                if not e in r[n]:
                    r[n][e] = set([(0,1)])

    for n in g:
        for pair in itertools.combinations(g[n],2):

            if pair[1] in r[pair[0]]:
                r[pair[0]][pair[1]].add((2,0))
            else:
                r[pair[0]][pair[1]] = set([(2,0)])

            if pair[0] in r[pair[1]]:
                r[pair[1]][pair[0]].add((2,0))
            else:
                r[pair[1]][pair[0]] = set([(2,0)])

    return r

def increment_u(G_star, G_u):
    # directed edges
    G_un = directed_inc(G_star,G_u)
    # bidirected edges
    G_un = bidirected_inc(G_un,G_u)
    return G_un
def undersample(G, u):
    Gu = G
    for i in range(u):
        Gu = increment_u(G, Gu)
    return Gu
def all_undersamples(G_star,steps=5):
    glist = [G_star]
    while True:
        g = increment_u(G_star, glist[-1])
        if ecj.isSclique(g): return glist # superclique convergence
        # this will (may be) capture DAGs and oscillations
        if g in glist: return glist
        glist.append(g)
    return glist

def graph2adj(G):
    n = len(G)
    A = scipy.zeros((n,n), dtype=np.int8)
    for v in G:
        A[int(v)-1,[int(w)-1 for w in G[v] if (0,1) in G[v][w] ]] = 1
    return A
def graph2badj(G):
    n = len(G)
    A = scipy.zeros((n,n), dtype=np.int8)
    for v in G:
        A[int(v)-1,[int(w)-1 for w in G[v] if (2,0) in G[v][w] ]] = 1
    return A

def adj2graph(A):
    names = [str(i) for i in range(1,A.shape[0]+1)]
    G = {}
    for name in names:
        G[name] = {}
    for i in range(0,A.shape[0]):
        for name in map(str, np.where(A[i,:]==1)[0]+1):
            G[str(i+1)][name]=set([(0,1)])
    return G

def adjs2graph(A,B):
    names = [str(i) for i in range(1,A.shape[0]+1)]
    G = {}
    for name in names:
        G[name] = {}
    for i in range(A.shape[0]):
        for name in map(str, np.where(A[i,:]==1)[0]+1):
            G[str(i+1)][name]=set([(0,1)])
            
    for i in range(B.shape[0]):
        for j in range(B.shape[1]):
            if B[i,j]:
                if str(j+1) in G[str(i+1)]:
                    G[str(i+1)][str(j+1)].add((2,0))
                else:
                    G[str(i+1)][str(j+1)] = set([(2,0)])                    
    return G

def g2vec(g):
    A = graph2adj(g)
    B = graph2badj(g)
    return np.r_[A.flatten(),B[np.triu_indices(B.shape[0])]]

def vec2adj(v,n):
    A = np.zeros((n,n))
    B = np.zeros((n,n))    
    A[:] = v[:n**2].reshape(n,n)
    B[np.triu_indices(n)] = v[n**2:]
    B = B+B.T
    return A,B

def vec2g(v,n):
    A,B = vec2adj(v,n)
    return adjs2graph(A,B)

# tried mutable ctypes buffer - not faster :(
def graph2str(G):
    n = len(G)
    d = {((0,1),):'1', ((2,0),):'0',((2,0),(0,1),):'0',((0,1),(2,0),):'0'}
    A = ['0']*(n*n)
    for v in G:
        for w in G[v]:
            A[n*(int(v)-1)+int(w)-1] = d[tuple(G[v][w])]
    return ''.join(A)

def graph2bstr(G):
    n = len(G)
    d = {((0,1),):'0', ((2,0),):'1',((2,0),(0,1),):'1',((0,1),(2,0),):'1'}
    A = ['0']*(n*n)
    for v in G:
        for w in G[v]:
            A[n*(int(v)-1)+int(w)-1] = d[tuple(G[v][w])]
    return ''.join(A)


def adj2num(A):
    s = reduce(lambda y,x: y+str(x), 
               A.flatten().tolist(),'')
    return int(s,2)

def g2num(G): return int(graph2str(G),2) #adj2num(graph2adj(G))

def bg2num(G): return int(graph2bstr(G),2)#adj2num(graph2badj(G))
def ug2num(G): return (g2num(G),bg2num(G))#(adj2num(graph2adj(G)),adj2num(graph2badj(G)))

def num2adj(num,n):
    l = list(bin(num)[2:])
    l = ['0' for i in range(0,n**2 - len(l))] + l
    return scipy.reshape(map(int, l),[n,n])
def add_bd_by_adj(G,adj):
    c = 0
    for e in adj:
        for v in range(len(e)):
            if e[v] == 1:
                try:
                    G[str(c+1)][str(v+1)].add((2,0))
                except KeyError:
                    G[str(c+1)][str(v+1)] = set([(2,0)])                    
        c += 1
    return G

def tuple2graph(t,n):
    g = num2CG(t[0],n)
    return add_bd_by_adj(g, num2adj(t[1],n))

def call_undersamples(G_star,steps=5):
    glist = [G_star]
    while True:
        g = increment_u(G_star, glist[-1])
        if g in glist: return glist
        glist.append(g)
    return glist

def compact_call_undersamples(G_star,steps=None):
    glist = [ug2num(G_star)]
    lastgraph = G_star
    while True:
        g = increment_u(G_star, lastgraph)
        if ug2num(g) in glist: return glist
        glist.append(ug2num(g))
        lastgraph = g
    return glist

def cc_undersamples(G_star,steps=1):
    glist = [ug2num(G_star)]
    lastgraph = G_star
    for i in xrange(steps):
        g = increment_u(G_star, lastgraph)
        n = ug2num(g)
        if n in glist: return []
        glist.append(n)
        lastgraph = g       
    return glist[-1]

def compatible(d1,d2):
    idx = scipy.where(scipy.array([[r==l for l in d2] for r in d1]))
    return idx

def compat(G):
    n = len(G)
    num = g2num(G)
    #sample all the graph for gStar
    star_l = all_undersamples(G)
    hits = {}
    #brute force all graphs
    for i in range(0, 2**(n**2)):
        tmpG = num2CG(i,n)
        tmp_l = all_undersamples(tmpG)
        c = compatible(tmp_l,star_l)
        if len(sum(c))>0: hits[i] = c
    return hits

def icompat(i,nodes):
    print i
    g = num2CG(i,nodes)
    return compat(g)

def ilength(i,nodes):
    print i
    g = num2CG(i,nodes)
    return len(call_undersamples(g))

def iall(i,nodes):
    print i
    g = num2CG(i,nodes)
    return compact_call_undersamples(g)

def cc_all(i,nodes,steps):
    #print i
    g = num2CG(i, nodes)
    return cc_undersamples(g, steps=steps)


def make_rect(l):
    max_seq = max(map(len,l))
    nl = []
    for e in l:
        e += [e[-1]] * (max_seq - len(e))
        nl.append(e)
    return nl

def uniqseq(l):
    s = []
    ltr = map(lambda *a: list(a), *l)
    for i in range(len(ltr)):
        s.append(len(np.unique(ltr[i])))
        
def loadgraphs(fname):
    g = zkl.load(fname)
    return g

def savegraphs(l,fname):
    zkl.save(l,fname)

def jason2graph(g):
    r = {}
    d = {1: set([(0,1)]),
         2: set([(2,0)]),
         3: set([(0,1),(2,0)]) }
    for head in g:
        r[head] = {}
        for tail in g[head]:
            r[head][tail] = d[g[head][tail]]
    return r

def graph2jason(g):
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

# talking about extra edges on top of the ring
def dens2edgenum(d, n=10): return int(d*n**2)-n
def edgenum2dens(e, n=10): return np.double(e+n)/n**2

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

