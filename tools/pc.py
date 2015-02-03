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


    def not_cindependent(y,x, condset=[], pval=0.05):
        yd = data[n+int(y)-1,:].T
        if condset: condset=map(lambda k: int(k)-1,condset)
        X  = data[[int(x)-1]+condset,:].T
        X  = sm.add_constant(X)
        est  = sm.OLS(yd,X).fit()
        print est.pvalues
        return est.pvalues[0] > pval

    g  = gk.superclique(n)
    gtr= bfu.gtranspose(g)
    el = gk.edgelist(g)
    
    counter = 0
    to_remove = []
    while counter <= n-2:
        for e in el:
            ind = True
            for S in [j for j in iter.combinations(
                    [k for k in gtr[e[1]] if e[0] != k],counter)]:
                if not_cindependent(e[1], e[0], condset=list(S)):
                    ind = False
                    break
            if not ind:
                to_remove.append(e)
                el.remove(e)
        counter += 1
        
    for e in to_remove:
        g[e[0]][e[1]].remove((0,1))
    to_remove = []
    for v in g:
        for w in g[v]:
            if not g[v][w]:
                to_remove.append((v,w))
    for v in to_remove:
        g[v[0]].pop(v[1], None)
    return g
