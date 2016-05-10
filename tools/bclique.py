import sys

sys.path.append('./tools/')
import traversal as trv
import bfutils as bfu
import numpy as np
from ortools.constraint_solver import pywrapcp

dens = 0.3
N = 5
k = bfu.dens2edgenum(dens, N)

solver = pywrapcp.Solver("b-clique")

# generate a random graph
g = bfu.ringmore(N, k)


def edge_exists(i, j, g):
    return str(j + 1) in g[str(i + 1)]


# declare variables
edges = []
for i in range(N):
    e = []
    for j in range(N):
        # if edge_exists(i, j, g):
        e.append(solver.IntVar(0, 1, "%i -> %i" % (i, j)))
    edges.append(e)

parents = []
children = []
for i in range(N):
    parents.append(solver.IntVar(0, 1, "%i" % (i)))
    children.append(solver.IntVar(0, 1, "%i" % (i)))


def numkids(i, e=edges, s=solver):
    return s.Sum([e[i][k] for k in range(N)])


def numparents(i, e=edges, s=solver):
    return s.Sum([e[k][i] for k in range(N)])


def allkids(i, c=children, e=edges, s=solver):
    n = len(c)
    el = [e[i][k] for k in range(n)]
    dl = [e[i][k] * c[k] for k in range(n)]
    return s.Sum(c) == s.Sum(dl)


def allpars(j, p=parents, e=edges, s=solver):
    n = len(p)
    el = [e[k][j] for k in range(n)]
    dl = [e[k][j] * p[k] for k in range(n)]
    return s.Sum(p) == s.Sum(dl)


# declare constraints
solver.Add(solver.Sum(parents) > 0)
solver.Add(solver.Sum(children) > 0)

for i in range(N):
    solver.Add((parents[i] == 1) == (solver.Sum([edges[i][k] for k in range(N)]) == numkids(i)))
    solver.Add((children[i] == 1) == (solver.Sum([edges[k][i] for k in range(N)]) == numparents(i)))

    for j in range(N):
        # edge existence constraints
        if not edge_exists(i, j, g):
            solver.Add(0 == edges[i][j])
        else:
            solver.Add(edges[i][j] == parents[i]*children[j])
#            solver.Add((edges[i][j] == 1) == allpars(j))
# stop

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


for s in range(num_solutions):
    #if s>5:
    #    break
    print " ----------------- ", s
    qval = [collector.Value(s, edges[i][j]) for i in range(N) for j in range(N)]
    pval = [collector.Value(s, parents[i]) for i in range(N)]
    cval = [collector.Value(s, children[i]) for i in range(N)]
    print np.asarray([pval, cval])
    for i in range(len(qval)):
        if qval[i]:
            e = np.unravel_index(i, [N, N])
            print e[0], "->", e[1]
    print "---"
    print
