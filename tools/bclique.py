import sys

sys.path.append('./tools/')
import traversal as trv
import bfutils as bfu
import numpy as np
from ortools.constraint_solver import pywrapcp


dens = 0.4
N = 5
k = bfu.dens2edgenum(dens, N)

solver = pywrapcp.Solver("b-clique")

# generate a random graph
g = bfu.ringmore(N,k)

def edge_exists(i, j, g):
    return str(j+1) in g[str(i+1)]

# declare variables
edges = []
for i in range(N):
    e = []
    for j in range(N):
        #if edge_exists(i, j, g):
        e.append(solver.IntVar(0,1,"%i -> %i" % (i,j)))
    edges.append(e)

parents  = []
children = []
for i in range(N):
    parents.append(solver.IntVar(0,1,"%i" % (i)))
    children.append(solver.IntVar(0,1,"%i" % (i)))


def allkids(i, c = children, e = edges, s = solver):
    n = len(c)
    el = [e[i][k] for k in range(n)]
    dl = [e[i][k] * c[k] for k in range(n)]
    return s.Sum(el) == s.Sum(c) == s.Sum(dl)


def allpars(j, p = parents, e = edges, s = solver):
    n = len(p)
    el = [e[k][j] for k in range(n)]
    dl = [e[k][j] * p[k] for k in range(n)]
    return s.Sum(el) == s.Sum(p) == s.Sum(dl)


# declare constraints
for i in range(N):
    for j in range(N):
        # edge existance constraints
        if not edge_exists(i, j, g):
            solver.Add( 0 == edges[i][j] )
        else:
            solver.Add(parents[i] == edges[i][j])
            solver.Add(children[j] == edges[i][j])

            solver.Add(allkids(i) == edges[i][j])
            solver.Add(edges[i][j] == allpars(j))
#stop

# run the solver
solution = solver.Assignment()
solution.Add([edges[i][j] for i in range(N) for j in range(N)])
solution.Add(parents)
solution.Add(children)
collector = solver.AllSolutionCollector(solution)
solver.Solve(solver.Phase([edges[i][j] for i in range(N) for j in range(N)] + children + parents,
                          solver.CHOOSE_FIRST_UNBOUND,
                          solver.ASSIGN_MIN_VALUE),
                          [collector])
num_solutions = collector.SolutionCount()

# output solutions
print "num_solutions:", num_solutions
print "failures:", solver.Failures()
print "branches:", solver.Branches()
print "WallTime:", solver.WallTime()

if num_solutions > 0 and num_solutions < 5:
    for s in range(num_solutions):
        qval = [collector.Value(s, edges[i][j]) for i in range(N) for j in range(N)]
        for i in range(len(qval)):
            if qval[i]:
                e = np.unravel_index(i,[N,N])
                print e[0],"->",e[1]
        print
        print
