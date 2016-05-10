import sys

sys.path.append('./tools/')
import traversal as trv
import bfutils as bfu
import numpy as np
from ortools.constraint_solver import pywrapcp


N = 4
k = bfu.dens2edgenum(0.2, N)

solver = pywrapcp.Solver("b-clique")

# generate a random graph
g = bfu.ringmore(N,k)

# declare variables
edges = []
for i in range(N):
    e = []
    for j in range(N):
        e.append(solver.IntVar(0,1,"%i -> %i" % (i,j)))
    edges.append(e)


# run the solver
solution = solver.Assignment()
solution.Add([edges[i][j] for i in range(N) for j in range(N)])
collector = solver.AllSolutionCollector(solution)
solver.Solve(solver.Phase([edges[i][j] for i in range(N) for j in range(N)],
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