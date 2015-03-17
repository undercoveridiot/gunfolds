# This is a use-case of the tools in the tools directory. The example defines a graph and shows how to generate a figure that shows the graph at different undersampling rates. Running the file in python (python dbnplot.py) generates a figure in figures folder: shipfig.pdf


# system packages
import os, sys
import numpy as np
from random import random
sys.path.append('tools')
import zickle as zkl
# local packages
import dbn2latex as d2l
from bfutils import jason2graph    

def ring(n):
    g = {}
    g['1'] = {'1': {(0,1)}, '2': {(0,1)}}
    for i in range(2,n):
        g[str(i)] = {str(i+1): {(0,1)}}
    g[str(n)] = {'1': {(0,1)}}
    return g

def listplot(fname, mname='JJ', stl='', width=5):
    l = zkl.load(fname)
    y = min(width,len(l))
    x = np.int(np.ceil(len(l)/float(y)))
    d2l.matrix_list(l,y,x,R=2, w_gap=1, h_gap=2, mname=mname, stl=stl)


g = {'1': {'2': set([(0, 1)]), '7': set([(0, 1)])},
     '10': {'1': set([(0, 1)]), '5': set([(0, 1)]), '9': set([(0, 1)])},
     '2': {'3': set([(0, 1)]),
           '4': set([(0, 1)]),
           '6': set([(0, 1)]),
           '7': set([(0, 1)])},
     '3': {'4': set([(0, 1)])},
     '4': {'1': set([(0, 1)]), '4': set([(0, 1)]), '5': set([(0, 1)])},
     '5': {'10': set([(0, 1)]),
           '5': set([(0, 1)]),
           '6': set([(0, 1)]),
           '8': set([(0, 1)]),
           '9': set([(0, 1)])},
     '6': {'2': set([(0, 1)]), '7': set([(0, 1)])},
     '7': {'8': set([(0, 1)])},
     '8': {'4': set([(0, 1)]),
           '7': set([(0, 1)]),
           '8': set([(0, 1)]),
           '9': set([(0, 1)])},
     '9': {'1': set([(0, 1)]),
           '10': set([(0, 1)]),
           '2': set([(0, 1)]),
           '6': set([(0, 1)]),
           '7': set([(0, 1)])}}


# output file
foo = open('figures/shipfig_figure.tex', 'wb')
sys.stdout = foo

# generation of the output
g = {'1': {'2':set([(0,1)])},
     '2': {'3':set([(0,1)]), '1':set([(0,1)])},
     '3': {'4':set([(0,1)])},
     '4': {'1':set([(0,1)])},     
}
#d2l.matrix_fold(g,2,1,R=5, w_gap=1, h_gap=2, mname='TT1')

listplot('list.zkl', width=4)

sys.stdout = sys.__stdout__              # remember to reset sys.stdout!
foo.flush()
foo.close()
PPP = os.getcwd()
os.chdir('figures')
os.system('pdflatex --shell-escape shipfig.tex 2>&1 > /dev/null')
os.chdir(PPP)
