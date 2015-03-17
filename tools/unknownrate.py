#BFS implementation of subgraph and supergraph
#Gu to G1 algorithm
from itertools import combinations
from functools import wraps
import copy
import time
import sys,os

TOOLSPATH='./tools/'
sys.path.append(os.path.expanduser(TOOLSPATH))

import bfutils as bfu
import traversal as trv
import graphkit as gk
import comparison as cmp

def memo(func):
    cache = {}                        # Stored subproblem solutions
    @wraps(func)                      # Make wrap look like func
    def wrap(*args):                  # The memoized wrapper
        s = trv.gsig(args[0])        # Signature: just the g
        #s = tool.signature(args[0],args[2])# Signature: g and edges
        if s not in cache:            # Not already computed?
            cache[s] = func(*args)    # Compute & cache the solution
        return cache[s]               # Return the cached solution
    return wrap

def prune_conflicts(H, g, elist):
    """checks if adding an edge from the list to graph g causes a
    conflict with respect to H and if it does removes the edge
    from the list

    Arguments:
    - `H`: the undersampled graph
    - `g`: a graph under construction
    - `elist`: list of edges to check
    """
    masks  = []
    for e in elist:
        gk.addanedge(g,e)
        if gk.checkconflict(H,g):
            masks.append(False)
        else:
            masks.append(True)
        gk.delanedge(g,e)
    return [elist[i] for i in range(len(elist)) if masks[i]]

def eqclass(H):
    '''
    Find all graphs in the same equivalence class with respect to
    graph H and any undesampling rate.
    '''
    g = {n:{} for n in H}
    s = set()

    @memo
    def addedges(g,H,edges):
        if edges:
            nedges = prune_conflicts(H, g, edges)
            n = len(nedges)

            if n == 0: return None

            for i in range(n):
                gk.addanedge(g,nedges[i])
                if gk.checkequality(H,g): return trv.gsig(g)
                s.add(addedges(g,H,nedges[:i]+nedges[i+1:]))
                gk.delanedge(g,nedges[i])

    edges = gk.edgelist(gk.complement(g))
    addedges(g,H,edges)
    return s-set([None])

# these two functions come fromt his answer:
# http://stackoverflow.com/a/12174125
def set_bit(value, bit): return value | (1<<bit)
def clear_bit(value, bit): return value & ~(1<<bit)

def add2set_(ds, H):
    n = len(H)
    
    dsr = {}
    s = set()
    ss = set()

    for gnum in ds:#gset:
        g = bfu.num2CG(gnum, n)
        glist = []
        elist = []
        for e in ds[gnum]:
            if not e[1] in g[e[0]]:
                gk.addanedge(g,e)
                num = bfu.g2num(g)
                if not num in s:
                    au = bfu.call_undersamples(g)
                    if not gk.checkconflict(H, g, au=au):
                        glist.append(num)
                        elist.append(e)
                        s.add(num)
                        if gk.checkequality(H, g, au=au): ss.add(num)
                gk.delanedge(g,e)
        for gn in glist: dsr[gn] = elist

    return dsr, ss

def iteqclass(H):
    '''
    Find all graphs in the same equivalence class with respect to
    graph H and any undesampling rate.
    '''
    if cmp.isSclique(H):
        print 'not running on superclique'
        return None
    g = {n:{} for n in H}
    s = set()

    edges = gk.edgelist(gk.complement(g))

    ds = {bfu.g2num(g): edges}

    for i in range(len(H)**2):
        print i, len(ds)
        ds, ss = add2set_(ds, H)
        s = s | ss
        if not ds: break

    return s

def add2set(gset, elist, H):
    n = len(H)

    s = set()
    ss = set()

    eremove = {e: True for e in elist}

    for gnum in gset:
        g = bfu.num2CG(gnum, n)
        for e in elist:
            if not e[1] in g[e[0]]:
                gk.addanedge(g,e)
                num = bfu.g2num(g)
                if not num in s:
                    au = bfu.call_undersamples(g)
                    if not gk.checkconflict(H, g, au=au):
                        eremove[e] = False
                        s.add(num)
                        if gk.checkequality(H, g, au=au): ss.add(num)
                gk.delanedge(g,e)

    for e in eremove:
        if eremove[e]: elist.remove(e)

    return s, ss, elist

def eqclass_list(H):
    '''
    Find all graphs in the same equivalence class with respect to
    graph H and any undesampling rate.
    '''
    g = {n:{} for n in H}
    s = set()

    edges = gk.edgelist(gk.complement(g))
    #edges = prune_conflicts(H, g, edges)

    gset = set([bfu.g2num(g)])
    for i in range(len(H)**2):
        print i
        gset, ss, edges = add2set(gset, edges, H)
        s = s | ss
        if not edges: break

    return s

def main():
    g = bfu.ringmore(4,2);
    H = bfu.undersample(g,2);
    ss = eqclass(H)
    print ss

if __name__ == "__main__":
    main()
