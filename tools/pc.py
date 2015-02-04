import sys,os
TOOLSPATH='~/soft/src/dev/craft/gunfolds/tools/'
sys.path.append(os.path.expanduser(TOOLSPATH))

import traversal as trv
import bfutils as bfu
import graphkit as gk
import zickle as zkl
import numpy as np
import itertools as iter
import statsmodels.api as sm
import ipdb


def independent(y,X,condset=[], pval=0.05):
        X  = sm.add_constant(X)
        est  = sm.OLS(y,X).fit()
        return est.pvalues[1] > pval


def dpc(data, pval=0.05):
    n = data.shape[0]
    data = np.asarray(np.r_[data[:,:-1],data[:,1:]])


    def cind_(y,x, condset=[], pval=pval, shift=0):
        yd = data[n+int(y)-1,:].T
        X  = data[[shift+int(x)-1]+condset,:].T
        return independent(yd, X, condset=condset,pval=pval)

    def cindependent(y, x, counter, parents=[], pval=pval):
        for S in [j for j in iter.combinations(parents,counter)]:
            if cind_(y, x, condset=list(S), pval=pval):
                return True
        return False

    def bindependent(y, x, parents=[], pval=pval):
        return cind_(y, x, condset=parents, pval=pval, shift=n)

    def prune(elist, mask, g):
        for e in mask:
            g[e[0]][e[1]].remove((0,1))
            elist.remove(e)
        gk.clean_leaf_nodes(g)

    g  = gk.superclique(n)
    gtr= bfu.gtranspose(g)

    el = gk.edgelist(g)
    for counter in range(n):
        to_remove = []
        for e in el:
            ppp = [int(k)-1 for k in gtr[e[1]] if k != e[0]]
            ind = cindependent(e[1], e[0], counter, parents=ppp, pval=pval)
            if ind:
                to_remove.append(e)
                gtr[e[1]].pop(e[0],None)
        prune(el, to_remove, g)

    bel = [map(lambda k: str(k+1), x) for x in iter.combinations(range(n),2)]
    for e in bel:
        ppp = list(set(gtr[e[0]].keys()) | set(gtr[e[1]].keys()))
        ppp = map(lambda x: int(x)-1, ppp)
        if bindependent(e[0], e[1], parents=ppp, pval=pval):
            g[e[0]][e[1]].remove((2,0))
            g[e[1]][e[0]].remove((2,0))
    gk.clean_leaf_nodes(g)

    return g
