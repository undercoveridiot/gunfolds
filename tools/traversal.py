import sys
sys.path.append('/home/splis/soft/src/dev/craft/gunfolds/tools/')
from bfutils import increment_u, g2num, num2CG
#import ipdb

def compare(g2star,g2):
    #gse = set(edgelist(g2star))
    #return gse.issubset(edgelist(g2))
    for n in g2star:
        for h in g2star[n]:
            try:
                if not (0,1) in g2[n][h]:
                    return False
            except KeyError:
                    return False
    return True

def edgelist(g):
    l = []
    for n in g:
        l.extend([(n,e) for e in g[n] if (0,1) in g[n][e]])
    return l

def cloneempty(g):
    e = {}
    for n in g: e[n] = {}
    return e

def add2edges(g,e,p):
    try:
        e1 = g[e[0]][p]
    except KeyError:
        e1 = False
    try:
        e2 = g[p][e[1]]
    except KeyError:
        e2 = False
    try:
        g[e[0]][p].add((0,1))
    except KeyError:
        g[e[0]][p] = set([(0,1)])
    try:
        g[p][e[1]].add((0,1))
    except KeyError:
        g[p][e[1]] = set([(0,1)])
    return e1, e2

def deledges(g,e,p,e1,e2):
    if e[0] == e[1] == p:
        if e1:
            g[e[0]][p] = e1
        else:
            del g[e[0]][p]
    else:
        if e1:
            g[e[0]][p] = e1
        else:
            del g[e[0]][p]
        if e2:
            g[p][e[1]] = e2
        else:
            del g[p][e[1]]

def rotate(l): return l[1:] + l[:1]

def nodesearch(g, g2, edges, s):
    if edges:
        e = edges.pop()
        for n in g:
            e1, e2 = add2edges(g, e, n)
            if compare(increment_u(g,g), g2):
                nodesearch(g,g2,edges,s)
            deledges(g,e,n,e1,e2)
        edges.insert(0,e)
    else:
        s.add(g2num(g))                


def edgesearch(g, g2, edges, s):
    print ''.join(['-' for x in edges])
    for e in edges:
        print edges
        nodesearch(g, g2, edges, s)
        edges = rotate(edges)

def eqc(g2):
    '''
    computes all g1 that are in the equivalence class for g2
    '''
    # find all directed g1's not conflicting with g2
    edges = edgelist(g2)
    g = cloneempty(g2)
    s = set()
    nodesearch(g,g2,edges,s)
    return s
    # prune those that are not consistent with the bidirected edge structure


