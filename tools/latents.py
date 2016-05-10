import sys,os
import sys, os
sys.path.append('./tools/')
import traversal as trv
import bfutils as bfu
import graphkit as gk
import numpy as np

def bclique(g, n):
    gt = bfu.gtranspose(g)
    checkedp = set()
    checkedc = set()

    bc = {}
    bc['children'] = set([n])
    bc['parents'] = set(gt[n].keys())

    Qc = []
    Qp = gt[n].keys()
    checkedc.add(n)

    while Q:
        p = Q.pop()
        if not p in checkedp:
            checkedp.add(p)
            Qc = g[p].keys()
