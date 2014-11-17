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

g = {'1': {'2':1,'16':1},
     '2': {'3':1},
     '3': {'4':1},
     '4': {'5':1,'4':1},     
     '5': {'6':1},
     '6': {'7':1,'2':1},
     '7': {'8':1},
     '8': {'9':1},
     '9': {'10':1},
     '10': {'11':1,'10':1},
     '11': {'12':1},
     '12': {'13':1},
     '13': {'14':1},
     '14': {'15':1},
     '15': {'1':1},
     '16': {'17':1},
     '17': {'18':1},
     '18': {'19':1,'18':1},
     '19': {'20':1},
     '20': {'21':1},
     '21': {'22':1},
     '22': {'23':1},
     '23': {'24':1,'45':1},
     '24': {'25':1},
     '25': {'26':1},
     '26': {'27':1},
     '27': {'28':1},
     '28': {'29':1},
     '29': {'30':1,'29':1},
     '30': {'31':1},
     '31': {'32':1},
     '32': {'33':1},
     '33': {'34':1},
     '34': {'35':1},
     '35': {'36':1,'2':1},
     '36': {'37':1},
     '37': {'38':1},
     '38': {'39':1},
     '39': {'40':1},
     '40': {'41':1},
     '41': {'42':1},
     '42': {'43':1,'17':1},
     '43': {'44':1},
     '44': {'45':1},
     '45': {'46':1},
     '46': {'47':1},
     '47': {'48':1},
     '48': {'49':1},
     '49': {'50':1},
     '50': {'51':1},
     '51': {'52':1},
     '52': {'53':1},
     '53': {'54':1},
     '54': {'55':1},
     '55': {'56':1},
     '56': {'57':1},
     '57': {'58':1},
     '58': {'59':1},
     '59': {'60':1},
     '60': {'16':1}
}

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

# fname = 'list.zkl'
# l = zkl.load(fname)

# y = min(5,len(l))
# x = np.int(np.ceil(len(l)/float(y)))

# output file
foo = open('figures/shipfig_figure.tex', 'wb')
sys.stdout = foo

# generation of the output
#ww = len(l)
g = {'1': {'2': {(0, 1)}, '8': {(0, 1)}},
     '10': {'1': {(0, 1)}},
     '2': {'3': {(0, 1)}},
     '3': {'10': {(0, 1)}, '3': {(0, 1)}, '4': {(0, 1)}},
     '4': {'1': {(0, 1)}, '5': {(0, 1)}, '7': {(0, 1)}},
     '5': {'6': {(0, 1)}},
     '6': {'10': {(0, 1)}, '7': {(0, 1)}},
     '7': {'3': {(0, 1)}, '4': {(0, 1)}, '8': {(0, 1)}},
     '8': {'4': {(0, 1)}, '7': {(0, 1)}, '9': {(0, 1)}},
     '9': {'10': {(0, 1)}}}
d2l.matrix_fold(g,2,1,R=5, w_gap=1, h_gap=2, mname='TT1')
#for i in range(3,10):
#    d2l.matrix_fold(ring(i),ww,1,R=2.5, w_gap=1, h_gap=2, mname='TT'+str(i-1), 
#                    stl=', below=5cm of TT'+str(i-2)+'.west,anchor=west')

#d2l.matrix_list(l,y,x,R=2, w_gap=1, h_gap=2)


sys.stdout = sys.__stdout__              # remember to reset sys.stdout!
foo.flush()
foo.close()
PPP = os.getcwd()
os.chdir('figures')
os.system('pdflatex --shell-escape shipfig.tex 2>&1 > /dev/null')
os.chdir(PPP)
