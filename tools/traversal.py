import sys
sys.path.append('/home/splis/soft/src/dev/craft/gunfolds/tools/')
from bfutils import increment_u, g2num, num2CG
from functools import wraps
import numpy as np
import random
import ipdb
import ecj
import itertools, copy
import munkres
from comparison import graph2nx
from networkx import strongly_connected_components

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

def edgelist(g): # directed
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

def purgepath(path, l):
    for i in range(1,len(path)-1):
        l.remove((path[i],path[i+1]))

def next_or_none(it):
    try:
        n = it.next()
    except StopIteration:
        return None
    return n

def try_till_d_path(g,d,order=None):
    k = []
    i = 1
    while not k:
        if order:
            k = next_or_none(length_d_paths(g,order(i),d))
        else:
            k = next_or_none(length_d_paths(g,str(i),d))
        i += 1
        if i > len(g): return []
    return k

def try_till_path(g):
    d = len(g)-1
    gx = graph2nx(g)
    sccl = [x for x in strongly_connected_components(gx)]
    # take the largest
    ln = [len(x) for x in sccl]
    idx = np.argsort(ln)
    d = len(sccl[idx[-1]])-1
    sccl = [sccl[i] for i in idx[::-1]]
    order = [item for sublist in sccl for item in sublist]
    k = []
    while not k:
        k = try_till_d_path(g,d)
        d -= 1
        if d == 0: return []
    return k

def gpurgepath(g,path):
    for i in range(1,len(path)-1):
        del g[path[i]][path[i+1]]


def vedgelist(g):
    """ Return a list of tuples for edges of g and forks
    """
    l = []
    el = edgelist(g)

    gc = copy.deepcopy(g)
    for i in range(16):

        k = try_till_path(gc)
        if len(k) < 5: break
        if k:
            l.append(('2',)+tuple(k))
            purgepath(l[-1],el)
            gpurgepath(gc,l[-1])
        else:
            break

    bl = bedgelist(g)
    for n in g:
        c = [e for e in g[n] if (0,1) in g[n][e]]# all children
        if len(c) == 1:
            if (n,c[0]) in el:
                l.append((n,c[0]))
                el.remove((n,c[0]))
        elif len(c) > 1:
            r = set()
            for p in [x for x in itertools.combinations(c,2)]:
                if (not p in bl) and (n,p[0]) in el and (n,p[1]) in el:
                    l.append((n,)+p)
                    el.remove((n,p[0]))
                    el.remove((n,p[1]))
                    r.add(p[0])
                    r.add(p[1])
            for e in r: c.remove(e)
            #b = [tuple([n]+i) for i in chunks(c,2)]

            b = []
            for p in chunks(c,2):

                try:
                    if (n,p[0]) in el:
                        if (n,p[1]) in el:
                            b.append(tuple([n]+p))
                            el.remove((n,p[0]))
                            el.remove((n,p[1]))
                        else:
                            b.append((n,p[0]))
                            el.remove((n,p[0]))
                except IndexError:
                    b.append(tuple([n]+p))
                    el.remove(tuple([n]+p))

            l.extend(b)
    r = twoedges(l)
    A, singles = makechains(r)

    if singles:
        B, singles = makesinks(singles)
    else:
        B, singles = [], []

    l = longpaths(l)+threedges(l) + A + B + singles
    return l

def twoedges(l):  return [e for e in l if len(e)==2]
def threedges(l): return [e for e in l if len(e)==3]
def longpaths(l): return [e for e in l if len(e)>3 and e[0]=='2']
def makechains(l):
    """ Greedily construct 2 edge chains from edge list
    """
    ends = {e[1]:e for e in l}
    starts = {e[0]:e for e in l}
    r = []
    singles = []
    while l:
        e = l.pop()
        if e[1] in starts and e[0] != e[1] and starts[e[1]] in l:
            r.append(('0', e[0],)+starts[e[1]])
            l.remove(starts[e[1]])
        elif e[0] in ends and e[0] != e[1] and ends[e[0]] in l:
            r.append(('0',)+ends[e[0]]+(e[1],))
            l.remove(ends[e[0]])
        else:
            singles.append(e)
    return r, singles
def makesink(es): return ('1', es[0][0],) + es[1]
def makesinks(l):
    """ Greedily construct 2 edge sinks ( b->a<-c ) from edge list
    """
    sinks = {}
    for e in l:
        if e[1] in sinks:
            sinks[e[1]].append(e)
        else:
            sinks[e[1]] = [e]
    r = []
    singles = []
    for e in sinks:
        if len(sinks[e])>1:
            for es in chunks(sinks[e],2):
                if len(es)==2:
                    r.append(makesink(es))
                else:
                    singles.append(es[0])
        else:
            singles.append(sinks[e][0])
    return r, singles

def selfloop(n,g):
    return n in g[n]

def selfloops(l,g):
    return reduce(lambda x,y: x and y, map(lambda x: selfloop(x,g), l))

def checkbedges(v,bel,g2):
    r = []
    for e in bel:
        if e == tuple(v[1:]) and not selfloops(e, g2):
            r.append(e)
        if e == (v[2],v[1]) and not selfloops(e, g2):
            r.append(e)
    for e in r: bel.remove(e)
    return bel

def checkedge(e, g2):
    if e[0] == e[1]:
        l = [n for n in g2 if n in g2[n]]
        s = set()
        for v in g2[e[0]]: s=s.union(g2[e[0]][v])
        if not (2,0) in s:
            l.remove(e[0])
        return l
    else:
        return g2.keys()

def single_nodes(v,g2):
    """ Returns a list of singleton nodes allowed for merging with v
    """
    l = [(n,n) for n in g2 if not n in v and len(g2[n])>1]
    return l

def checkvedge(v, g2):
    """ Nodes to check to merge the virtual nodes of v ( b<-a->c )
    """
    l = bedgelist(g2)
    if (v[1],v[2]) in l:
        l = single_nodes(v,g2) + checkbedges(v,l,g2)
        for n in v:
            if n in g2[n]: l.append((n,n))
    else:
        l = checkbedges(v,l,g2)
    return list(set(l))

def checkAedge(v, g2):
    """ Nodes to check to merge the virtual nodes of A ( b->a<-c )
    """
    l = []
    # try all pairs but the sources
    for pair in itertools.combinations(g2,2):
        if pair == (v[1],v[2]): continue
        if pair == (v[2],v[1]): continue
        l.append(pair)
        l.append(pair[::-1])
    for n in g2:
        l.append((n,n))
    return l

def checkcedge(c, g2):
    """ Nodes to check to merge the virtual nodes of c ( a->b->c )
    """
    l = edgelist(g2)
    return list(set(l))

def checkApath(p, g2):
    sl = [x for x in g2 if selfloop(x,g2)]
    d = len(p) - 2
    l = []
    for n in g2:
        l.extend([x for x in length_d_loopy_paths(g2, n, d, p[1:])])
    k = prunepaths_1D(g2, p, l)
    return k


def isedge(v):  return len(v) == 2 # a->b
def isvedge(v): return len(v) == 3 # b<-a->c
def isCedge(v): return len(v) == 4 and v[0] == '0' # a->b->c
def isAedge(v): return len(v) == 4 and v[0] == '1'# b->a<-c
def isApath(v):  return len(v) >= 4 and v[0] == '2'# a->b->...->z

def checkable(g2):
    d = {}
    g = cloneempty(g2)

    vlist = vedgelist(g2)
    for v in vlist:
        if isvedge(v):
            d[v] = checkvedge(v,g2)
        elif isCedge(v):
            d[v] = checkcedge(v,g2)
        elif isAedge(v):
            d[v] = checkAedge(v,g2)
        elif isApath(v):
            d[v] = checkApath(v,g2)
        else:
            d[v] = checkedge(v,g2)

    # check if some of the otherwise permissible nodes still fail
    f = [(add2edges, del2edges),
         (addavedge,delavedge),
         (addacedge,delacedge),
         (addaAedge,delaAedge),
         (addapath,delapath)]
    for e in d:
        adder, remover = f[min(4,len(e))-2+min(max(3,len(e))-3,1)*int(e[0])]#f[len(e)-2]
        for n in d[e]:
            mask = adder(g,e,n)
            if not isedgesubset(increment(g), g2):
                d[e].remove(n)
            remover(g,e,n,mask)

    return d

def inorder_check2(e1, e2, j1, j2, g2):
    g = cloneempty(g2) # the graph to be used for checking
    f = [(add2edges, del2edges),
         (addavedge,delavedge),
         (addacedge,delacedge),
         (addaAedge,delaAedge),
         (addapath,delapath)]

    adder1, remover1 = f[min(4,len(e1))-2+min(max(3,len(e1))-3,1)*int(e1[0])]
    adder2, remover2 = f[min(4,len(e2))-2+min(max(3,len(e2))-3,1)*int(e2[0])]

    d = {}
    for c1 in j1: # for each connector
        mask1 = adder1(g,e1,c1)
        d[c1] = set()
        for c2 in j2:
            mask2 = adder2(g,e2,c2)
            if isedgesubset(increment(g), g2):
                d[c1].add(c2)
            remover2(g,e2,c2,mask2)
        remover1(g,e1,c1,mask1)
    return d

def inorder_check3(e1, e2, e3, j1, j2, j3, g2):
    g = cloneempty(g2) # the graph to be used for checking
    f = [(add2edges, del2edges),
         (addavedge,delavedge),
         (addacedge,delacedge)]

    adder1, remover1 = f[len(e1)-2+(max(3,len(e1))-3)*int(e1[0])]
    adder2, remover2 = f[len(e2)-2+(max(3,len(e2))-3)*int(e2[0])]
    adder3, remover3 = f[len(e3)-2+(max(3,len(e3))-3)*int(e3[0])]

    d = {}
    for c1 in j1: # for each connector
        mask1 = adder1(g,e1,c1)
        d[c1] = set()
        for c2 in j2:
            mask2 = adder2(g,e2,c2)
            if isedgesubset(increment(g), g2):
                d[c1].add(c2)
            remover2(g,e2,c2,mask2)
        remover1(g,e1,c1,mask1)
    return d

def del_empty(d):
    l = [e for e in d]
    for e in l:
        if d[e]==set(): del d[e]
    return d
def inorder_checks(g2, gg):
    ee = [e for e in gg] # to preserve the order
    #cds = conformanceDS(g2, ee)
    #oo = new_order(g2, ee, repeats=100, cds=None)
    #ee = oo[0]
    #random.shuffle(ee)
    d = {} # new datastructure
    d[ee[0]] = {('0'):gg[ee[0]]}
    for i in range(len(ee)-1):
        d[ee[i+1]] = del_empty(inorder_check2(ee[i], ee[i+1],
                                              gg[ee[i]], gg[ee[i+1]], g2))
    return ee, d

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
    if not mask[1]: g[p].pop(e[1], None)

def addavedge(g,v,b):
    mask = [b[0] in g[v[0]], b[1] in g[v[0]],
            v[1] in g[b[0]], v[2] in g[b[1]]]

    g[v[0]][b[0]] = g[v[0]][b[1]] = g[b[0]][v[1]] = g[b[1]][v[2]] = set([(0,1)])
    return mask

def delavedge(g,v,b,mask):
    if not mask[0]: g[v[0]].pop(b[0], None)
    if not mask[1]: g[v[0]].pop(b[1], None)
    if not mask[2]: g[b[0]].pop(v[1], None)
    if not mask[3]: g[b[1]].pop(v[2], None)

def addaAedge(g,v,b):
    mask = [b[0] in g[v[1]], b[1] in g[v[2]],
            v[3] in g[b[0]], v[3] in g[b[1]]]

    g[v[1]][b[0]] = g[v[2]][b[1]] = g[b[0]][v[3]] = g[b[1]][v[3]] = set([(0,1)])
    return mask

def delaAedge(g,v,b,mask):
    if not mask[0]: g[v[1]].pop(b[0], None)
    if not mask[1]: g[v[2]].pop(b[1], None)
    if not mask[2]: g[b[0]].pop(v[3], None)
    if not mask[3]: g[b[1]].pop(v[3], None)

def addapath(g,v,b):
    mask = []
    s = set([(0,1)])
    for i in range(len(b)):
        mask.append(b[i] in g[v[i+1]])
        mask.append(v[i+2] in g[b[i]])

    for i in range(len(b)):
        g[v[i+1]][b[i]] = g[b[i]][v[i+2]] = s

    return mask

def delapath(g, v, b, mask):
    for i in range(len(b)):
        if not mask[2*i]: g[v[i+1]].pop(b[i], None)
        if not mask[2*i+1]:g[b[i]].pop(v[i+2], None)

def prunepaths_1D(g2, path, conn):
    c = []
    g = cloneempty(g2)
    for p in conn:
        mask = addapath(g,path,p)
        if isedgesubset(increment(g), g2): c.append(tuple(p))
        delapath(g,path,p,mask)
    return c

def addacedge(g,v,b): # chain
    mask = [b[0] in g[v[1]], v[2] in g[b[0]],
            b[1] in g[v[2]], v[3] in g[b[1]]]

    g[v[1]][b[0]] = set([(0,1)])
    g[v[2]][b[1]] = set([(0,1)])
    g[b[0]][v[2]] = set([(0,1)])
    g[b[1]][v[3]] = set([(0,1)])

    return mask

def delacedge(g,v,b,mask):
    if not mask[0]: g[v[1]].pop(b[0], None)
    if not mask[1]: g[b[0]].pop(v[2], None)
    if not mask[2]: g[v[2]].pop(b[1], None)
    if not mask[3]: g[b[1]].pop(v[3], None)

def rotate(l): return l[1:] + l[:1] # rotate a list
def density(g): return len(edgelist(g))/np.double(len(g)**2)

def esig(l,n):
    '''
    turns edge list into a hash string
    '''
    z = len(str(n))
    n = map(lambda x: ''.join(map(lambda y: y.zfill(z),x)), l)
    n.sort()
    n = ''.join(n[::-1])
    return int('1'+n)

def gsig(g):
    '''
    turns input graph g into a hash string using edges
    '''
    return g2num(g)#esig(edgelist(g))

def signature(g, edges): return (gsig(g),esig(edges,len(g)))

def memo(func):
    cache = {}                        # Stored subproblem solutions
    @wraps(func)                      # Make wrap look like func
    def wrap(*args):                  # The memoized wrapper
        s = signature(args[0],args[2])# Signature: g and edges
        if s not in cache:            # Not already computed?
            cache[s] = func(*args)    # Compute & cache the solution
        return cache[s]               # Return the cached solution
    return wrap

def memo2(func):
    cache = {}                        # Stored subproblem solutions
    @wraps(func)                      # Make wrap look like func
    def wrap(*args):                  # The memoized wrapper
        s = signature(args[0],args[2].keys())# Signature: g and edges
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
            ln = g2.keys()
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

    f = [(add2edges, del2edges),
         (addavedge,delavedge),
         (addacedge,delacedge)]
    @memo2 # memoize the search
    def nodesearch(g, g2, edges, s):
        if edges:
            #key, checklist = edges.popitem()
            key = random.choice(edges.keys())
            checklist = edges.pop(key)
            adder, remover = f[len(key)-2]
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

def v2g22g1(g2, capsize=None):
    '''
    computes all g1 that are in the equivalence class for g2
    '''
    if ecj.isSclique(g2):
        print 'Superclique - any SCC with GCD = 1 fits'
        return set([-1])
    f = [(add2edges, del2edges),
         (addavedge,delavedge),
         (addacedge,delacedge),
         (addaAedge,delaAedge),
         (addapath,delapath)]

    @memo2 # memoize the search
    def nodesearch(g, g2, edges, inlist, order, s, cds):
        if edges:
            key = order.pop(0)
            checklist = edges.pop(key)

            if inlist[0] in checklist:
                if len(inlist) == 1:
                    tocheck = checklist[inlist[0]]
                else:
                    tocheck = conformant(cds, inlist)

                #if tocheck:
                adder, remover = f[min(4,len(key))-2+min(max(3,len(key))-3,1)*int(key[0])]
                #f[len(key)-2]
                for n in tocheck:
                    mask = adder(g,key,n)
                    if isedgesubset(increment(g), g2):
                        r = nodesearch(g,g2,edges,[n]+inlist,order,s, cds)
                        if r and increment(r)==g2:
                            s.add(g2num(r))
                            if capsize and len(s)>capsize:
                                raise ValueError('Too many elements')
                    remover(g,key,n,mask)

            order.insert(0,key)
            edges[key] = checklist

        else:
            return g

    # find all directed g1's not conflicting with g2
    chlist = checkable(g2)
    order, d = inorder_checks(g2,chlist)
    cds = conformanceDS(g2, order)
    g = cloneempty(g2)

    s = set()
    try:
        nodesearch(g,g2,d,['0'],order,s, cds)
    except ValueError:
        s.add(0)
    return s


def conformanceDS(g2, order):
    gg = checkable(g2)
    CDS = {}
    CDS[0] = set(gg[order[0]])
    for x in itertools.combinations(range(len(order)),2):
        d = del_empty(inorder_check2(order[x[0]], order[x[1]],
                                     gg[order[x[0]]], gg[order[x[1]]], g2))
        if not x[1] in CDS:
            CDS[x[1]] = {}
            CDS[x[1]][x[0]] = d
        else:
            CDS[x[1]][x[0]] = d

    return pruneCDS(CDS)

def conformant(cds, inlist):
    if inlist[len(inlist)-2] in cds[len(inlist)-1][0]:
        s = cds[len(inlist)-1][0][inlist[len(inlist)-2]]
    else:
        return set()
    for i in range(1,len(inlist)-1):
        if inlist[len(inlist)-i-2] in cds[len(inlist)-1][i]:
            s = s.intersection(cds[len(inlist)-1][i][inlist[len(inlist)-i-2]])
        else:
            return set()
    return s

def pruneCDS(cds):
    cds[0] = cds[0].intersection(cds[1][0].keys())
    #for i in range(2, len(cds)):
    return cds


def cost_matrix(g2, order, cds=None):
    gg = checkable(g2)
    if cds==None: cds = conformanceDS(g2, order)
    m = 100000.0*np.eye(len(order))
    for x in itertools.combinations(range(len(order)),2):
#        d = del_empty(inorder_check2(order[x[0]], order[x[1]],
#                                     gg[order[x[0]]], gg[order[x[1]]], g2))
        d = cds[x[1]][x[0]]
        s = sum([len(d[k]) for k in d])
        r = len(gg[order[x[0]]])*len(gg[order[x[1]]])
        m[x[0],x[1]] = np.double(s)#/r
        m[x[1],x[0]] = np.double(s)#/r
    return m, cds

def cost_path(p, W):
    x = [0]+p[:-1]
    return W[x,p]
def path_weight(p, W):
    return sum(cost_path(p,W))

def new_order(g2, order, repeats = 100, cds=None):
    M, cds = cost_matrix(g2,order,cds)
    p = range(1,M.shape[0])
    mnw = path_weight(p,M) #- 10*sum(np.diff(cost_path(p,M)))
    mnp = p[:]
    for i in range(repeats):
        random.shuffle(p)
        pw = path_weight(p,M) #- 10*sum(np.diff(cost_path(p,M)))
        if pw < mnw:
            mnw = pw
            mnp = p[:]
    return [order[i] for i in [0]+mnp], mnw, cds

def length_d_paths(G, s, d):
    """
    Iterate over nodes in G reachable from s in exactly d steps
    """
    def recurse(G, s, d, path=[]):
        if d == 0:
            yield path
            return

        for u in G[s]:
            if G[s][u] == set([(2,0)]) or u in path: continue
            for v in recurse(G, u, d-1, path+[u]):
                yield v

    for u in recurse(G, s, d, [s]):
            yield u

def subfinder(mylist, pattern):
    matches = []
    for i in range(len(mylist)):
        if mylist[i] == pattern[0] and mylist[i:i+len(pattern)] == pattern:
            matches.append((i,pattern))
    return matches

def make_triplets(path, selfloops):
    tri = {}
    for u in selfloops:
        idx = path.index(u)
        tri[u] = path[idx-1:idx+2]
    return tri

def length_d_loopy_paths1(G, s, d, path=None, ok2loop=[]):
    """
    Iterate over nodes in G reachable from s in exactly d steps
    """
    candidates = [x for x in length_d_paths(G,s,d)]
    triplets = make_triplets(path, ok2loop)
    for u in ok2loop:
        ext = []
        for p in candidates:
            sub = subfinder(p, triplets[u])
            if sub:
                idx = sub[0][0]
                ext.append(p[:idx+1]+[u,]+p[idx+1:-1])
        candidates.extend(ext)

    if path[0] in ok2loop:
        ext = []
        for p in candidates:
            if p[:2]==path[:2]:
                ext.append([p[0]]+p[:-1])
        candidates.extend(ext)

    if path[-1] in ok2loop:
        ext = []
        for p in candidates:
            if p[-2:]==path[-2:]:
                ext.append(p[1:]+[p[-1]])
        candidates.extend(ext)

    return candidates

def edge_increment_ok(s,m,e,g,g2):
    """
    s - start,
    m - middle,
    e - end
    """
    # directed edges
    for u in g[m]:
        if not u in g2[s] or not (0,1) in g2[s][u]:
            return False
    for u in g:
        if s in g[u] and (not m in g2[u] or not (0,1) in g2[u][m]):
            return False
        if m in g[u] and (not e in g2[u] or not (0,1) in g2[u][e]):
            return False

    # bidirected edges
    for c in g[s]:
        if c == s: continue
        if c == m: continue
        try:
            if not (2,0) in g2[c][m]:
                return False
        except KeyError:
            return False
    for c in g[m]:
        if c == m: continue
        if c == e: continue
        try:
            if not (2,0) in g2[c][e]:
                return False
        except KeyError:
            return False
    return True

def length_d_loopy_paths(G, s, dt, p):
    """
    Iterate over nodes in G reachable from s in exactly d steps
    """
    g1 = cloneempty(G)

    def recurse(g, g2, s, d, path=[]):

        if edge_increment_ok(p[-d-2],s,p[-d-1],g,g2):
            
            if d == 0:
                yield path
                return

            mask = add2edges(g,(p[-d-2],p[-d-1]),s)
            for u in g2[s]:
                if g2[s][u] == set([(2,0)]): continue
                for v in recurse(g, g2, u, d-1, path+[u]):
                    yield v
            del2edges(g,(p[-d-2],p[-d-1]),s, mask)
            
    for u in recurse(g1, G, s, dt-1, [s]):
            yield u
