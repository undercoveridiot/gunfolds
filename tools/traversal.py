import sys
sys.path.append('/home/splis/soft/src/dev/craft/gunfolds/tools/')
from bfutils import increment_u, g2num, num2CG
from functools import wraps
import numpy as np
import random
#import ipdb
import ecj
import itertools, copy

def pincrement(g,g2,n1,n2,n3):
    r = copy.deepcopy(g2)
    rg= ecj.tr(g)
    for k in rg[n1]:
        if n2 in r[k]: 
            r[k][n2].add((0,1))
        else:
            r[k][n2] = set([(0,1)])        
    for k in g[n2]:
        if k in r[n1]: 
            r[n1][k].add((0,1))
        else:
            r[n1][k] = set([(0,1)])
    for k in g[n3]:
        if k in r[n2]: 
            r[n2][k].add((0,1))
        else:
            r[n2][k] = set([(0,1)])
    for k in g[n1]:
        if k != n2:
            if k in r[n2]:
                r[n2][k].add((2,0))
            else:
                r[n2][k] = set([(2,0)])
            if n2 in r[k]:
                r[k][n2].add((2,0))
            else:
                r[k][n2] = set([(2,0)])
    for k in g[n2]:
        if k != n3:
            if k in r[n3]:
                r[n3][k].add((2,0))
            else:
                r[n3][k] = set([(2,0)])
            if n3 in r[k]:
                r[k][n3].add((2,0))
            else:
                r[k][n3] = set([(2,0)])
    return r

def increment(g):
    '''
    undersample g by 2
    only works for g1 to g2 directed
    '''
    r = {n:{} for n in g}
    for n in g:
        for h in g[n]:
            for e in g[h]:
                r[n][e] = set([(0,1)])
    for n in g:
        for pair in itertools.permutations(g[n],2):
            if pair[1] in r[pair[0]]:
                r[pair[0]][pair[1]].add((2,0))
            else:
                r[pair[0]][pair[1]] = set([(2,0)])
    return r

def isedgesubset(g2star,g2):
    '''
    check if g2star edges are a subset of those of g2
    '''
    for n in g2star:
        for h in g2star[n]:
            if h in g2[n]:
                #if not (0,1) in g2[n][h]:
                if not g2star[n][h].issubset(g2[n][h]):
                    return False
            else:
                    return False
    return True

def edgelist(g):
    '''
    return a list of tuples for edges of g
    '''
    l = []
    for n in g:
        l.extend([(n,e) for e in g[n] if (0,1) in g[n][e]])
    return l

def cloneempty(g): return {n:{} for n in g} # return a graph with no edges

def add2edges(g,e,p):
    '''
    break edge e[0] -> e[1] into two pieces 
    e[0] -> p and p -> e[1]
    and add them to g
    '''
    e1 = p in g[e[0]]
    e2 = e[1] in g[p]
    g[e[0]][p] = set([(0,1)])
    g[p][e[1]] = set([(0,1)])
    return e1, e2

def deledges(g,e,p,e1,e2):
    '''
    restore the graph as it was before adding e[0]->p and p->e[1]
    '''
    if not e1: del g[e[0]][p]
    if not e2 and e[1] in g[p]: del g[p][e[1]]

def rotate(l): return l[1:] + l[:1] # rotate a list
def density(g): return len(edgelist(g))/np.double(len(g)**2)

def esig(l):
    '''
    turns edge list into a hash string
    '''
    n = map(lambda x: '.'.join(x), l)
    n.sort()
    n = ','.join(n)
    return n

def gsig(g):
    '''
    turns input graph g into a hash string using edges
    '''
    return esig(edgelist(g))

def signature(g, edges): return (gsig(g),esig(edges))

def memo(func):
    cache = {}                        # Stored subproblem solutions
    @wraps(func)                      # Make wrap look like func
    def wrap(*args):                  # The memoized wrapper
        s = signature(args[0],args[2])# Signature: g and edges
        if s not in cache:            # Not already computed?
            cache[s] = func(*args)    # Compute & cache the solution
        return cache[s]               # Return the cached solution
    return wrap

def g22g1(g2, capsize=None):
    '''
    computes all g1 that are in the equivalence class for g2
    '''
    if ecj.isSclique(g2):
        print 'Superclique - any SCC with GCD = 1 fits'
        return set([-1])

    single_cache = {}

    @memo # memoize the search
    def nodesearch(g, g2, edges, s):
        if edges:
            e = edges.pop()
            #gg = increment(g)
            ln = [n for n in g2]
            random.shuffle(ln)
            for n in ln:
                if (n,e) in single_cache: continue
                e1, e2 = add2edges(g,e,n)
                #if isedgesubset(pincrement(g,gg,e[0],n,e[1]), g2):
                if isedgesubset(increment(g), g2):
                    r = nodesearch(g,g2,edges,s)
                    if r and increment(r)==g2:
                        s.add(g2num(r))    
                        if capsize and len(s)>capsize:
                            raise ValueError('Too many elements in eqclass')
                deledges(g,e,n,e1,e2)
            edges.append(e)
        else:
            return g
    # find all directed g1's not conflicting with g2
    n = len(g2)
    edges = edgelist(g2)
    random.shuffle(edges)
    g = cloneempty(g2)

    for e in edges:
        for n in g2:        
            e1, e2 = add2edges(g,e,n)
            if not isedgesubset(increment(g), g2):
                single_cache[(n,e)] = False
            deledges(g,e,n,e1,e2)

    s = set()
    try:
        nodesearch(g,g2,edges,s)
    except ValueError:
        s.add(0)
    return s
