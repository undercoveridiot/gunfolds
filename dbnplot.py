# system packages
import os, sys
import numpy as np
from random import random
import pickle as pkl
sys.path.append('tools')

# local pakcages
import dbn2latex as d2l
    
g = {
    '1': {'2': {(0, 1)}, '5': {(0, 1)}},
    '2': {'3': {(0, 1)}},
    '3': {'4': {(0, 1)}},
    '4': {'1': {(0, 1)}},
    '5': {'6': {(0, 1)}},
    '6': {'1': {(0, 1)}},
}

g = {
    '1': {'2': {(0, 1)}, '3': {(0, 1)}},
    '2': {'5': {(0, 1)}},
    '3': {'4': {(0, 1)}},
    '4': {'6': {(0, 1)}},
    '5': {'6': {(0, 1)}},
    '6': {'1': {(0, 1)}, '3': {(0, 1)}},
}

# fname = 'list.pkl'
# f = open(fname,'r')
# l = pkl.load(f)
# f.close()

#y = min(5,len(l))
#x = np.int(np.ceil(len(l)/float(y)))


# output file
foo = open('figures/shipfig_figure.tex', 'wb')
sys.stdout = foo

# generation of the output
ww = 1
d2l.gmatrix_fold(g,1,ww,R=2, w_gap=1, h_gap=2, mname='TT1')
#d2l.matrix_list(l,y,x,R=2,w_gap=2,h_gap=2)


sys.stdout = sys.__stdout__              # remember to reset sys.stdout!
foo.flush()
foo.close()
PPP = os.getcwd()
os.chdir('figures')
os.system('pdflatex --shell-escape shipfig.tex 2>&1 > /dev/null')
os.chdir(PPP)
