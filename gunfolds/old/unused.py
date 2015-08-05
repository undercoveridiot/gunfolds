""" This is where functions from old experiments and previous data formats have been moved
    to in order to clarify the current code.  These remain for reference and in case they
    become useful again.  No garuntee all imports are here to use these functions."""
import numpy as np
import scipy



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










""" from ecj.py """

def bfs_print_tree(tree, r):
    """
    A modified single list solution
    """
    Q = []
    idx = 1
    print str(idx) + ': ' + r
    Q.extend(tree[r])
    while Q:
        idx += 1
        print str(idx) + ':',
        for u in range(0, len(Q)):
            e = Q.pop(0)
            print e,
            Q.extend(tree[e])
        print ''


def bfs_dict(tree, r):
    """
    Bob's suggested dictionary based solution
    """
    D = {}
    idx = 1
    D[idx] = [r]
    while D[idx]:
        idx += 1
        D[idx] = []
        for u in D[idx - 1]:
            D[idx].extend(tree[u])
    D.pop(idx)  # the last dictionary element is empty - must go
    for idx in D:
        print str(idx) + ': ' + ' '.join(D[idx])


def clrbi(G):
    for v in G:
        d = []
        for u in G[v]:
            try:
                G[v][u].remove((edge_type['bidirected'], 0))
                if len(G[v][u]) == 0:
                    d.append(u)
            except KeyError:
                pass
        for e in d:
            G[v].pop(e)


def ecj(G, s, sccs=set()):           # elementary circuit by Johnson
    blocked = {v: False for v in G}  # unblock all
    B = {v: [] for v in G}
    stack = []

    def unblock(u):
        blocked[u] = False
        for w in B[u]:
            B[u].remove(w)
            if blocked[w]:
                unblock(w)

    def circuit(v, stack):
        f = False
        stack.append(v)
        blocked[v] = True
        for u in G[v]:
            if u == s:
                f = True
                # print stack
                sccs.add(len(stack))
            elif not blocked[u]:
                if circuit(u, stack):
                    f = True
        if f:
            unblock(v)
        else:
            for w in G[v]:
                if v not in B[w]:
                    B[w].append(v)
        stack.pop()
        return f
    circuit(s, stack)
    return sccs


def ecj_loops(G, s, sl=set()):           # elementary circuit by Johnson
    blocked = {v: False for v in G}  # unblock all
    B = {v: [] for v in G}
    stack = []

    def unblock(u):
        blocked[u] = False
        for w in B[u]:
            B[u].remove(w)
            if blocked[w]:
                unblock(w)

    def circuit(v, stack):
        f = False
        stack.append(v)
        blocked[v] = True
        for u in G[v]:
            if u == s:
                f = True
                # print scipy.sort(stack)
                sl.add(tuple(scipy.sort(stack)))
            elif not blocked[u]:
                if circuit(u, stack):
                    f = True
        if f:
            unblock(v)
        else:
            for w in G[v]:
                if v not in B[w]:
                    B[w].append(v)
        stack.pop()
        return f
    circuit(s, stack)
    return sl


def lcm(a, b):
    return a * b / gcd(a, b)


def chmatch(n, m, delta):
    m, n = scipy.sort([n, m])
    sq = scipy.mod(range(n, lcm(n, m) + 1, n), m)
    return scipy.mod(delta, m) in sq


def allpaths(G, s, g, S=[]):
    if S is None:
        S = []
    S.append(s)
    if s == g:
        print S
    else:
        for u in G[s]:
            if u in S:
                continue
            allpaths(G, u, g, S)
    S.remove(s)


def lz_ecj(G, s, sccs=set()):           # elementary circuit by Johnson
    blocked = {v: False for v in G}    # unblock all
    B = {v: set() for v in G}
    stack = []

    def unblock(u):
        blocked[u] = False
        for w in B[u]:
            if blocked[w]:
                unblock(w)
        B[u].clear()

    def circuit(v, stack):
        stack.append(v)
        blocked[v] = True
        for u in G[v]:
            if u == s:
                print 'bottom'
                unblock(v)
                yield len(stack)
            elif not blocked[u]:
                print 'recurse'
                for x in circuit(u, stack):
                    unblock(v)
                    yield x
            else:
                print 'unmet'
                for w in G[v]:
                    B[w].add(v)
        stack.pop()
#    circuit(s,stack)
    for v in circuit(s, stack):
        yield v


def exist_equal_paths(h, G, a, b):
    Sa, Sb, Dff = set(), set(), set()
    ag = iterate_allpaths(G, h, a, 0, [], True)
    bg = iterate_allpaths(G, h, b, 0, [], True)
    for v in izip_longest(ag, bg):
        print v
        Sa.add(v[0])
        Sb.add(v[1])
        if v[0] in Sb or v[1] in Sa:
            return True
    return False


def iexist_equal_paths(h, G, a, b):
    """ checks if there  exist exact length paths from the  head node to there
        nodes at question  by iterative deepining to avoid  oing through all
        paths """
    Sa, Sb = set(), set()
    Pa, Pb = [], []
    ag = iddfs(G, h, a)
    bg = iddfs(G, h, b)
    for v in izip(ag, bg):
        print v
        Sa.add(len(v[0]))
        Pa.append(v[0])
        Sb.add(len(v[1]))
        Pb.append(v[1])
        if len(v[0]) in Sb or len(v[1]) in Sa:
            return True
    return False


def iterate_allpaths(G, s, g, d=0, S=[], c=True):
    if S is None:
        S = []
    S.append(s)
    d += 1
    if s == g:
        if c:
            yield d - 1
        else:
            yield list(S)
    else:
        for u in G[s]:
            if u in S:
                continue
            for v in iterate_allpaths(G, u, g, d, S, c):
                yield v
    S.remove(s)


def iddfs(G, s, g):  # iterative depening DFS paths
    yielded = set()

    def recurse(G, s, g, d, S=None):
        if s not in yielded:
            yielded.add(s)
        if d == 0:
            return
        if S is None:
            S = []
        S.append(s)
        if s == g:
            yield list(S)
        else:
            for u in G[s]:
                if u in S:
                    continue
                for v in recurse(G, u, g, d - 1, S):
                    yield v
        S.remove(s)
    n = len(G)
    for d in range(n):
        # if len(yielded) == n: break
        for u in recurse(G, s, g, d):
            yield u


def ecj_compat(G, p1, p2):
    n = len(p1)
    m = len(p2)
    p2, p1 = [[p1, p2][i] for i in scipy.argsort([n, m])]
    m, n = scipy.sort([n, m])
    delta = n - m
    if not delta:
        return True  # redundant check for 0
    if has_unit_cycle(G, p2):
        return True  # equivalent
    # if the shorter path does not have cycles they are not compatible
    # if the  longer path does not  have cycles: check  if the shorter
    #                                            path    has    cycles
    #                                            divisible by delta

    # otherwise start checking
    print p1, p2, n, m, delta


def wc(n):
    n = n * 3
    a = {str(v): set([str(v + 1), str(v + 2), str(v + 3)]) for v in range(1, n, 3)}
    b = {str(v): set([str(v + 2)]) for v in range(2, n, 3)}
    c = {str(v): set([str(v + 1)]) for v in range(3, n, 3)}
    a.update(b)
    a.update(c)
    a.update({str(n): set()})
    return a


def residue_table(a):
    """ Frobenius number from here: http://cgi.gladman.plus.com/wp/?page_id=563 """
    n = [0] + [None] * (a[0] - 1)
    for i in range(1, len(a)):
        d = gcd(a[0], a[i])
        for r in range(d):
            try:
                nn = min(n[q] for q in range(r, a[0], d) if n[q] != None)
            except:
                continue
            if nn != None:
                for c in range(a[0] // d):
                    nn += a[i]
                    p = nn % a[0]
                    nn = min(nn, n[p]) if n[p] != None else nn
                    n[p] = nn
    return n


def frobenius_number(a):
    return max(residue_table(sorted(a))) - min(a)


def isJclique(G):
    """ Jianyu does not use bidirected edges """
    return (sum([len(G[w].keys()) for w in G]) == len(G) ** 2)


def sample_graph(graph_g, steps=5):
    graph_g_list = [graph_g]
    for i in range(0, steps):
        g = increment_u(graph_g, graph_g_list[-1])
        graph_g_list.append(g)
    return graph_g_list


def reached_at_step(G, s, d):
    """
    Iterate over nodes in G reachable from s in exactly d steps
    """
    yielded = set()

    def recurse(G, s, d, B=None):
        if d == 0:
            if s not in yielded:  # this avoids yielding duplicates
                yielded.add(s)
                yield s
            return
        if B is None:
            B = []  # black - backed out of this path
        for u in G[s]:
            # if u in B: continue
            if G[s][u] == (edge_type['bidirected'], 0):
                continue
            for v in recurse(G, u, d - 1, B):
                yield v
        B.append(s)
    for u in recurse(G, s, d):
        yield u


def d_trek(h, G, a, b, d):
    """
    Does there exist a trek with head h connecting a and b in d steps.
    """
    return set([a, b]).issubset(reached_at_step(G, h, d))


def d_biegde(G, a, b, d):
    """
    Do  a and  b  become connected  by  a bidirectional  edge after  d
    undersamples
    """
    for i in range(1, d + 1):
        for u in G:
            if d_trek(u, G, a, b, i):
                return True
    return False


def undersample(G, d, bid=True):
    """

    """
    N = {}
    for u in G:
        N.update({u: {v: set([(0, 1)]) for v in reached_at_step(G, u, d + 1)}})
    if bid:
        items = G.keys()
        for i in range(len(items)):
            for j in range(i + 1, len(items)):
                u, v = items[i], items[j]
                if d_biegde(G, u, v, d):
                    try:
                        N[u][v].add((edge_type['bidirected'], 0))
                    except KeyError:
                        N[u].update({v: set([(edge_type['bidirected'], 0)])})
                    try:
                        N[v][u].add((edge_type['bidirected'], 0))
                    except KeyError:
                        N[v].update({u: set([(edge_type['bidirected'], 0)])})
    return N







""" from graphkit.py """



def ibedgelist(g):  # directed iterator
    '''
    iterate over the list of tuples for edges of g
    '''
    for v in g:
        for w in g[v]:
            if (2, 0) in g[v][w]:
                yield (v, w)


def edgenumber(g):
    return sum([sum([len(g[y][x]) for x in g[y]]) for y in g])


def iedgelist(g):  # directed iterator
    '''
    iterate over the list of tuples for edges of g
    '''
    for v in g:
        for w in g[v]:
            if (0, 1) in g[v][w]:
                yield (v, w)


def list2dbn(l):
    """ convert list of edge presences/absences (0,1) to a DBN graph
    """
    n = scipy.sqrt(len(l))
    l = scipy.reshape(map(int, l), [n, n])
    G = ecj.adj2DBN(l)
    return G


def rnd_dbn(n):
    return list2dbn(rnd_edges(n))


def sp_rnd_dbn(n, maxindegree=3):
    '''
    a sparse random DBN graph
    '''
    l = sp_rnd_edges(n)
    return list2dbn(l)


def rnd_edges(n):
    """ generate a random uniformly distributed mask
    """
    rnum = std_random.getrandbits(n ** 2)
    l = list(bin(rnum)[2:])
    l = ['0' for i in range(0, n ** 2 - len(l))] + l
    return l


def rnd_adj(n, maxindegree=5):
    l = scipy.zeros([n, n])
    for u in range(0, n):
        cap = scipy.random.randint(min([n, maxindegree + 1]))
        idx = scipy.random.randint(n, size=cap)
        l[u, idx] = 1
    return l


def sp_rnd_edges(n, maxindegree=5):
    '''
    a sparse set of random edges
    '''
    l = rnd_adj(n, maxindegree=maxdegree)
    return scipy.reshape(l, n ** 2)


def list2CG(l):
    """ convert list of edge presences/absences (0,1) to a compressed
    graph (CG) representation
    """
    n = scipy.sqrt(len(l))
    l = scipy.reshape(map(int, l), [n, n])
    G = adj2graph(l)
    return G


def rnd_cg(n):
    return list2CG(rnd_edges(n))


def sp_rnd_CG(n, maxindegree=3, force_connected=False):
    l = sp_rnd_edges(n, maxindegree=maxindegree)
    cg = list2CG(l)
    if force_connected:
        while not connected(cg):
            cg = list2CG(sp_rnd_edges(n, maxindegree=maxindegree))
    return cg


def adj2graph(A):
    G = {str(i): {} for i in range(1, A.shape[0] + 1)}
    idx = np.where(A == 1)
    for i in range(len(idx[0])):
        G['%i' % (idx[0][i] + 1)]['%i' % (idx[1][i]+1)] = set([(0, 1)])
    return G


def emptyG(n):
    A = [[0 for j in range(n)] for i in range(n)]
    return adj2graph(np.asarray(A))


def fullG(n):
    A = [[1 for j in range(n)] for i in range(n)]
    return adj2graph(np.asarray(A))


def CG2uCG(cg):
    """
    convert to an undirected graph
    """
    G = {}
    for u in cg:
        G[u] = cg[u].copy()
    for u in cg:
        for v in cg[u]:
            G[v][u] = cg[u][v]
    return G


def connected(cg):
    n = len(cg)
    return sum(1 for _ in ecj.traverse(CG2uCG(cg), '1')) == n


def fork_mismatch(g):
    be = bedgelist(g)
    benum = len(be) / 2
    forknum = 0
    for v in g:
        fn = len([w for w in g[v] if (0, 1) in g[v][w]])
        forknum += fn * (fn - 1) / 2.
    if benum < len(g) * (len(g) - 1) / 2.:
        return (forknum - benum) > benum
    else:
        return False







""" from pc.py """

def addallb(g):
    n = len(g)
    for i in range(n):
        for j in range(n):
            if str(j+1) in g[str(i+1)]:
                g[str(i+1)][str(j+1)].add((2, 0))
            else:
                g[str(i+1)][str(j+1)] = set([(2, 0)])
    return g




""" from simpleloops.py """


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





""" from traversal.py """

def next_or_none(it):
    try:
        n = it.next()
    except StopIteration:
        return None
    return n
