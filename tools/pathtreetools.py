import sys

sys.path.append('./tools/')
from pathtree import PathTree
from ortools.constraint_solver import pywrapcp


def ptloopnum(pt):
    """
    Given a PathTree object returns the number of loop in it
    :param pt: PathTree object
    :return: number of loops (n)
    """
    def ptn(pt,n=0):
        for e in pt.loopset:
            if type(e) is int:
                n += 1
                continue
            n += ptn(e, n=1)
        return n
    return ptn(pt)


def ptelement(pt, w):
    n = next(iter(pt.preset))

    def sumloops(pt, w):
        n = 0
        ls = list(pt.loopset)
        for i in range(len(ls)):
            if type(ls[i]) is int:
                n += w[i] * ls[i]
                continue
            n += w[i]*next(iter(ls[i].preset)) \
                 + min(1, w[i]) * sumloops(ls[i], w[len(ls):])
        return n
    return n + sumloops(pt, w)

def ptelements(pt, seqlen=100, verbose=False, maxloop=100):
    """
    Generate first `seqlen` elements from a pathtree
    :param pt: a path tree object from pathtree.py
    :param seqlen: number of elements to generate in ascending order
    :param verbose: whether to print debugging information
    :return: a list of elements
    """
    solver = pywrapcp.Solver("pt-elements")

    # declare variables
    weights = []
    N = ptloopnum(pt)
    for i in range(N):
        weights.append(solver.IntVar(0, maxloop, "w[%04i]" % i))

    # declare constraints
    #solver.Add()

    # run the solver
    solution = solver.Assignment()
    solution.Add(weights)
    db = solver.Phase(weights,
                      solver.CHOOSE_FIRST_UNBOUND,
                      solver.ASSIGN_MIN_VALUE)
    solver.NewSearch(db)

    num_solutions = 0
    els = set()
    while solver.NextSolution():
        w = [x.Value() for x in weights]
        num_solutions += 1
        els.add(ptelement(pt, w))
        if len(els) == seqlen:
            break
    solver.EndSearch()

    # output solutions
    if verbose:
        print "num_solutions:", num_solutions
        print "failures:", solver.Failures()
        print "branches:", solver.Branches()
        print "WallTime:", solver.WallTime()

    return list(els)

def isptelement(pt, element, verbose=False, maxloop=100):
    """
    Check if an integer element is in the weight set represented by the path tree
    :param pt: a path tree object from pathtree.py
    :param element: an integer to check for presence in the weight
    :param verbose: whether to print debugging information
    :return: True or False
    """
    solver = pywrapcp.Solver("pt-elements5")

    # declare variables
    weights = []
    N = ptloopnum(pt)
    for i in range(N):
        weights.append(solver.IntVar(0, maxloop, "w[%04i]" % i))

    # declare constraints
    solver.Add(element == ptelement(pt, weights))

    # run the solver
    solution = solver.Assignment()
    solution.Add(weights)
    db = solver.Phase(weights,
                      solver.CHOOSE_FIRST_UNBOUND,
                      solver.ASSIGN_MIN_VALUE)
    solver.NewSearch(db)

    solution_exists = False
    while solver.NextSolution():
        solution_exists = True
        break
    solver.EndSearch()

    # output solutions
    if verbose:
        print "failures:", solver.Failures()
        print "branches:", solver.Branches()
        print "WallTime:", solver.WallTime()

    return solution_exists


def pt_getelements(pt, num):
    i = 0
    s = set()
    while len(s) < num:
        if isptelement(pt, i, maxloop=10*num):
            s.add(i)
        i += 1
    return list(s)


def s2spt(s): # convert edge set to pt
    ss = set()
    for e in s:
        if type(e) is int:
            ss.add(PathTree({0}, pre={e}))
            continue
        ss.add(e)
    return ss


def spt_elements(spt, num):
    """
    Generate numbers from a set of PathTrees
    :param spt: set of PathTrees
    :param num: number of elements (from the first) to generate
    :return: list of num numbers
    """
    i = 0
    s = set()
    while len(s) < num:
        if issptelement(spt, i):
            s.add(i)
        i += 1
    return list(s)


def issptelement(spt, element):
    a = False
    for pt in s2spt(spt):
        a = a or isptelement(pt, element)
    return a