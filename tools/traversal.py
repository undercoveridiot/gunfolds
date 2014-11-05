import sys
sys.path.append('/home/splis/soft/src/dev/craft/gunfolds/tools/')
from bfutils import increment_u, g2num, num2CG
from functools import wraps
import numpy as np
import random
import ipdb
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

def bedgelist(g): # bidirected edge list with flips
    l = []
    for n in g:
        l.extend([tuple(sorted((n,e))) for e in g[n] if (2,0) in g[n][e]])
    l = list(set(l))
    l = l + map(lambda x: (x[1],x[0]), l)
    return l

def chunks(l, n):
    """ Yield successive n-sized chunks from l.
    """
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

def vedgelist(g):
    '''
    return a list of tuples for edges of g and forks
    '''
    l = []
    for n in g: # just all children
        c = [e for e in g[n] if (0,1) in g[n][e]]
        random.shuffle(c)
        if len(c) == 1:
            l.append((n,c[0]))
        elif len(c) > 1:
            b = [tuple([n]+i) for i in chunks(c,2)]
            l.extend(b)
    return l

def checkbedges(v,bel,g2):
    for e in bel:
        if e[0] in v and not e[0] in g2[e[0]]: 
            bel.remove(e)
        elif e[1] in v and not e[1] in g2[e[1]]:
            bel.remove(e)
    return bel

def checkedge(e, g2):
    l = [n for n in g2 if not n in e]
    for n in e: 
        if n in g2[n]: l.append(n)
    return l

def checkvedge(v, g2):
    bl = checkbedges(v,bedgelist(g2),g2)
    l = bl + [(n,n) for n in g2 if not n in v]
    for n in v: 
        if n in g2[n]: l.append((n,n))
    return list(set(l))
    
def isvedge(v): return len(v) == 3

def checkable(g2):
    d = {}
    vlist = vedgelist(g2)
    for v in vlist:
        if isvedge(v): 
            d[v] = checkvedge(v,g2)
        else:
            d[v] = checkedge(v,g2)            
    return d

def cloneempty(g): return {n:{} for n in g} # return a graph with no edges

def add2edges(g,e,p):
    '''
    break edge e[0] -> e[1] into two pieces
    e[0] -> p and p -> e[1]
    and add them to g
    '''
    mask = [p in g[e[0]], e[1] in g[p]]
    g[e[0]][p] = set([(0,1)])
    g[p][e[1]] = set([(0,1)])
    return mask

def del2edges(g,e,p,mask):
    '''
    restore the graph as it was before adding e[0]->p and p->e[1]
    '''
    if not mask[0]: g[e[0]].pop(p, None)
    if not mask[1] and e[1] in g[p]: g[p].pop(e[1], None)

def addavedge(g,v,b):
    mask = [b[0] in g[v[0]], b[1] in g[v[0]],
            v[1] in g[b[0]], v[2] in g[b[1]]]

    if mask[0]: g[v[0]][b[0]].add((0,1))
    if mask[1]: g[v[0]][b[1]].add((0,1))
    if mask[2]: g[b[0]][v[1]].add((0,1))
    if mask[3]: g[b[1]][v[2]].add((0,1))

    if not mask[0]: g[v[0]][b[0]] = set([(0,1)])
    if not mask[1]: g[v[0]][b[1]] = set([(0,1)])
    if not mask[2]: g[b[0]][v[1]] = set([(0,1)])
    if not mask[3]: g[b[1]][v[2]] = set([(0,1)])

    return mask

def delavedge(g,v,b,mask):
    if not mask[0]: g[v[0]].pop(b[0], None)
    if not mask[1]: g[v[0]].pop(b[1], None)
    if not mask[2]: g[b[0]].pop(v[1], None)
    if not mask[3]: g[b[1]].pop(v[2], None)

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
            #random.shuffle(ln)
            for n in ln:
                if (n,e) in single_cache: continue
                mask = add2edges(g,e,n)
                #if isedgesubset(pincrement(g,gg,e[0],n,e[1]), g2):
                if isedgesubset(increment(g), g2):
                    r = nodesearch(g,g2,edges,s)
                    if r and increment(r)==g2:
                        s.add(g2num(r))
                        if capsize and len(s)>capsize:
                            raise ValueError('Too many elements in eqclass')
                del2edges(g,e,n,mask)
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
            mask = add2edges(g,e,n)
            if not isedgesubset(increment(g), g2):
                single_cache[(n,e)] = False
            del2edges(g,e,n,mask)

    s = set()
    try:
        nodesearch(g,g2,edges,s)
    except ValueError:
        s.add(0)
    return s

def vg22g1(g2, capsize=None):
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
            key, checklist = edges.popitem()
            if len(key) == 2:
                adder = add2edges
                remover = del2edges
            else:
                adder = addavedge
                remover = delavedge
            for n in checklist:
                mask = adder(g,key,n)
                if isedgesubset(increment(g), g2):
                    r = nodesearch(g,g2,edges,s)
                    if r and increment(r)==g2:
                        s.add(g2num(r))
                        if capsize and len(s)>capsize:
                            raise ValueError('Too many elements')
                remover(g,key,n,mask)
            edges[key] = checklist
        else:
            return g

    # find all directed g1's not conflicting with g2
    n = len(g2)
    chlist = checkable(g2)
    g = cloneempty(g2)
    s = set()
    try:
        nodesearch(g,g2,chlist,s)
    except ValueError:
        s.add(0)
    return s
