import sys
sys.path.append('/home/splis/soft/src/dev/craft/gunfolds/tools/')
from bfutils import increment_u, g2num, num2CG
from functools import wraps
import numpy as np
import ipdb
import itertools

def increment(g):
    '''
    undersample g by 2
    only works for g1 to g2 directed
    '''
    r = {n:{} for n in g}
    for n in g:
        for h in g[n]:
            r[n] = {e:set([(0,1)]) for e in g[h]}
        for pair in itertools.permutations(g[n],2):
            try:
                r[pair[0]][pair[1]].add((2,0))
            except KeyError:
                r[pair[0]][pair[1]] = set([(2,0)])
    return r

def isedgesubset(g2star,g2):
    '''
    check if g2star edges are a subset of those of g2
    '''
    for n in g2star:
        for h in g2star[n]:
            try:
                #if not (0,1) in g2[n][h]:
                if not g2star[n][h].issubset(g2[n][h]):
                    return False
            except KeyError:
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

def gsig(g):
    n = map(lambda x: ''.join(x), edgelist(g))
    n.sort()
    n = ','.join(n)
    return n

def signature(g, edges):
    n = map(lambda x: ''.join(x), edgelist(g))
    n.sort()
    n = ','.join(n)
    l = map(lambda x: ''.join(x), edges)
    l.sort()
    e = ','.join(l)
    return (n,e)

def memo(func):
    cache = {}                        # Stored subproblem solutions
    @wraps(func)                      # Make wrap look like func
    def wrap(*args):                  # The memoized wrapper
        s = signature(args[0],args[2])# Signature: g and edges
        if s not in cache:            # Not already computed?
            cache[s] = func(*args)    # Compute & cache the solution
        return cache[s]               # Return the cached solution
    return wrap

def eqc(g2):
    '''
    computes all g1 that are in the equivalence class for g2
    '''
    @memo # memoize the search
    def nodesearch(g, g2, edges, s):
        if edges:
            e = edges.pop()
            for n in g:
                e1, e2 = add2edges(g,e,n)
                if isedgesubset(increment(g), g2):
                    r = nodesearch(g,g2,edges,s)
                    #if r: s.add(g2num(r))
                    if r and increment_u(r,r)==g2:
                        s.add(g2num(r))
                deledges(g,e,n,e1,e2)
            edges.append(e)
        else:
            return g
    # find all directed g1's not conflicting with g2
    n = len(g2)
    edges = edgelist(g2)
    g = cloneempty(g2)
    s = set()
    nodesearch(g,g2,edges,s)
    # prune those that are not consistent with the bidirected edge structure
    # now we  filter out those  g1's that are compatible  via directed
    # edges but not via bidirected
    #f = lambda x: increment_u(num2CG(x,n),num2CG(x,n))
    #l = [e for e in s if f(e) == g2]
    return s
