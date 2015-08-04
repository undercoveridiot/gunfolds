""" This is where functions from old experiments and previous data formats have been moved
    to in order to clarify the current code.  These remain for reference and in case they
    become useful again. """



""" from bfutils.py """

# tried mutable ctypes buffer - not faster :(
def graph2str(G):
    n = len(G)
    d = {((0, 1),): '1', ((2, 0),): '0', ((2, 0), (0, 1),): '1', ((0, 1), (2, 0),): '1'}
    A = ['0'] * (n * n)
    for v in G:
        for w in G[v]:
            A[n * (int(v, 10) - 1) + int(w, 10) - 1] = d[tuple(G[v][w])]
    return ''.join(A)


def graph2bstr(G):
    n = len(G)
    d = {((0, 1),): '0', ((2, 0),): '1', ((2, 0), (0, 1),): '1', ((0, 1), (2, 0),): '1'}
    A = ['0'] * (n * n)
    for v in G:
        for w in G[v]:
            A[n * (int(v, 10) - 1) + int(w, 10) - 1] = d[tuple(G[v][w])]
    return ''.join(A)


def adj2num(A):
    s = reduce(lambda y, x: y + str(x),
               A.flatten().tolist(), '')
    return int(s, 2)


def num2adj(num, n):
    l = list(bin(num)[2:])
    l = ['0' for i in range(0, n ** 2 - len(l))] + l
    return scipy.reshape(map(int, l), [n, n])


def add_bd_by_adj(G, adj):
    c = 0
    for e in adj:
        for v in range(len(e)):
            if e[v] == 1:
                try:
                    G[str(c + 1)][str(v + 1)].add((2, 0))
                except KeyError:
                    G[str(c + 1)][str(v + 1)] = set([(2, 0)])
        c += 1
    return G


def tuple2graph(t, n):
    g = num2CG(t[0], n)
    return add_bd_by_adj(g, num2adj(t[1], n))


def uniqseq(l):
    s = []
    ltr = map(lambda *a: list(a), *l)
    for i in range(len(ltr)):
        s.append(len(np.unique(ltr[i])))


def jason2graph(g):
    r = {}
    d = {1: set([(0, 1)]),
         2: set([(2, 0)]),
         3: set([(0, 1), (2, 0)])}
    for head in g:
        r[head] = {}
        for tail in g[head]:
            r[head][tail] = d[g[head][tail]]
    return r


def graph2jason(g):
    r = {}
    for head in g:
        r[head] = {}
        for tail in g[head]:
            if g[head][tail] == set([(0, 1)]):
                r[head][tail] = 1
            elif g[head][tail] == set([(2, 0)]):
                r[head][tail] = 2
            elif g[head][tail] == set([(0, 1), (2, 0)]):
                r[head][tail] = 3
    return r




