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

def dpc(data):
    n = data.shape[0]
    data = np.asarray(np.r_[data[:,:-1],data[:,1:]])


    def cindependent(y,x, condset=[], pval=0.05):
        yd = data[n+int(y)-1,:].T
        if condset: condset=map(lambda k: int(k)-1,condset)
        X  = data[[int(x)-1]+condset,:].T
        X  = sm.add_constant(X)
        est  = sm.OLS(yd,X).fit()
        return est.pvalues[0] > pval

    def prune(elist, mask, g):
        for e in mask:
            g[e[0]][e[1]].remove((0,1))
            elist.remove(e)
        to_remove = []
        for v in g:
            g[v] = {w for w in g[v] if g[v]}

    g  = gk.superclique(n)
    gtr= bfu.gtranspose(g)
    el = gk.edgelist(g)

    counter = 0
    to_remove = []
    while counter <= n-2:
        print el
        for e in el:
            ind = False
            for S in [j for j in iter.combinations(
                    [k for k in gtr[e[1]] if e[0] != k],counter)]:
                if cindependent(e[1], e[0], condset=list(S)):
                    ind = True
                    break
            if ind:
                to_remove.append(e)
                gtr[e[1]].pop(e[0],None)
        counter += 1
        prune(el, to_remove, g)
    return g
