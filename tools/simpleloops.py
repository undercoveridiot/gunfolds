import networkx
from bfutils import g2num, num2CG
from comparison import graph2nx
from ecj import undersample


def simple_loops(g, u):
    """
    iterator over the list of simple loops of graph g at the undersample rate u
    """
    gx = graph2nx(num2CG(g2num(undersample(g,u)), len(g)))
    for l in networkx.simple_cycles(gx):
        yield l

def print_loops(g, u):
    for l in simple_loops(g,u):
        print l

def ul(l):
    """
    returns list elements that are present only once
    """
    u, r = set(), set()
    for e in l:
        if e not in u:
            u.add(e)
        else:
            r.add(e)
    return u.difference(r)
