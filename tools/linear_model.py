import sys
sys.path.append('/home/splis/soft/src/dev/craft/gunfolds/tools/')

import ecj, bfutils
import graphkit as gk
from sympy.matrices import SparseMatrix
from progressbar import ProgressBar, Percentage, \
    Bar, RotatingMarker, ETA, FileTransferSpeed
from scipy import linalg, optimize
import numpy as np
import scipy
from statsmodels.tsa.api import VAR, SVAR
from matplotlib.mlab import amap

def symchol(M): # symbolic Cholesky
    B = SparseMatrix(M)
    t = B.row_structure_symbolic_cholesky()
    B = np.asarray(B)*0
    for i in range(B.shape[0]): B[i,t[i]] = 1
    return B

def G2SVAR(G):
    n = len(G)
    A,B = npG2SVAR(G)
    P,L,U = linalg.lu(B)
    A = linalg.inv(L).tolist()
    B = B.tolist()
    A = listplace(A, 0.0, 0.0)
    for i in range(0,n): A[i][i] = 1
    B = listplace(B, 0.0, 'e')
    for i in range(0,n): B[i][i] = 'e'
    return A,B,P

def G2AH(G):
    n = len(G)
    A,B = npG2SVAR(G)
    P,L,U = linalg.lu(B)
    A = linalg.inv(L).tolist()
    B = B.tolist()
    A = listplace(A, 0.0, 0.0)
    for i in range(0,n): A[i][i] = 1
    B = listplace(B, 0.0, 'e')
    for i in range(0,n): B[i][i] = 'e'
    return A,B,P

def bnf2CG(fname):
    d = eval(open(fname).read())
    G = {}
    for v in d:
        G[v] = {u: set((0,1)) for u in d[v]['pars']}
    G = ecj.tr(G)
    for v in G:
        ld = {u: set([(0,1)]) for u in G[v]}
        G[v] = ld
    return G

def npG2SVAR(G):
    n = len(G)
    A = [[0]*n]*n
    B = [[0]*n]*n
    for i in range(n): B[i][i] = 1
    
    for v in G:
        if G[v]:
            directed = [w for w in G[v] if (0,1) in G[v][w]]
            bidirected = [w for w in G[v] if (2,0) in G[v][w]]
            for w in   directed: A[int(w)-1][int(v)-1] = 1
            for w in bidirected: B[int(w)-1][int(v)-1] = 1

    A = np.asarray(A)
    B = symchol(B)
    return A,B

def x2M(x, A, B, aidx, bidx):
    A[aidx] = x[:len(aidx[0])]
    B[bidx] = x[len(aidx[0]):]
    #B[(bidx[1],bidx[0])] = x[len(aidx[0]):]    
    return A, B

def nllf(x, A, B, Y, aidx, bidx): # negative log likelihood
    A,B = x2M(x, A, B, aidx, bidx)
    T = Y.shape[1]
    X = Y[:,1:] - np.dot(A, Y[:,:-1])
    ldB = T*np.log(abs(1./linalg.det(B)))
    return ldB + 0.5*np.trace( np.dot(np.dot(B.T, B), np.dot(X,X.T)))

def nllf2(x, A, B, YY, XX, YX, T, aidx, bidx): # negative log likelihood
    A,B = x2M(x, A, B, aidx, bidx)
    AYX = np.dot(A, YX.T)
    S = YY - AYX - AYX.T + np.dot(np.dot(A,XX), A.T)
    ldB = T*np.log(abs(1./linalg.det(B)))
    return ldB + 0.5*np.dot(np.dot(B.T, B).T.flat,S.flat)
    #return ldB + 0.5*np.trace( np.dot(np.dot(B.T, B), S))

def VARbic(nllf, K, T):
    return 2*nllf + K*np.log(T)

def listplace(l, a, b):
    return [listplace(x,a,b) if not np.isscalar(x) else b if x != a  else x for x in l]

# -------------------------------------------------------------------
# data generation
# -------------------------------------------------------------------

def randweights(n, c=0.1, factor=9):
    rw = scipy.random.randn(n)
    idx = scipy.where(abs(rw) < factor*c)
    if idx:
        rw[idx] = rw[idx]+scipy.sign(rw[idx])*c*factor
    return rw

def transitionMarix(cg, minstrength=0.1):
    A = gk.CG2adj(cg)
    edges = scipy.where(A==1)
    A[edges] = randweights(edges[0].shape[0], c=minstrength)
    l = linalg.eig(A)[0]
    c = 0
    pbar = ProgressBar(widgets=['Something: ', Percentage(), ' '], maxval=10000).start()    
    while max(l*scipy.conj(l)) > 1:
        A[edges] = randweights(edges[0].shape[0], c=c)
        c += 1
        l = linalg.eig(A)[0]
        pbar.update(c)
    pbar.finish()        
    return A

def sampleWeights(n, minstrength=0.1):
    r = scipy.rand(n)
    s = minstrength/np.min(np.abs(r))
    r = s*r
    return r

def transitionMarix2(cg, minstrength=0.1):
    A = gk.CG2adj(cg)
    edges = scipy.where(A==1)
    A[edges] = sampleWeights(edges[0].shape[0], minstrength=minstrength)
    l = linalg.eig(A)[0]
    c = 0
    pbar = ProgressBar(widgets=['Something: ', Percentage(), ' '], maxval=10000).start()
    while max(l*scipy.conj(l)) > 1:
        A[edges] = sampleWeights(edges[0].shape[0], minstrength=minstrength)
        c += 1
        l = linalg.eig(A)[0]
        pbar.update(c)
    pbar.finish()
    return A

def drawsamplesLG(A, nstd=0.1, samples=100):
    n = A.shape[0]
    data = scipy.zeros([n, samples])
    data[:,0] = nstd*scipy.random.randn(A.shape[0])
    for i in range(1,samples):
        data[:,i] = scipy.dot(A,data[:,i-1]) \
                    + nstd*scipy.random.randn(A.shape[0])
    return data


def getAgraph(n, mp=2, st=0.5):
    keeptrying = True
    while keeptrying:
        G = gk.rnd_CG(n, maxindegree=mp, force_connected=True)
        try:
            A = transitionMarix2(G, minstrength=st)
            keeptrying = False
        except ValueError:
            print "!!! Unable to find strong links for a stable matrix !!!"
            print "*** trying a different graph"
    return {'graph':      G,
            'transition': A,
            'converges':  len(bfutils.call_undersamples(G))}

def getAring(n, density=0.1, st=0.5):
    keeptrying = True
    plusedges = bfutils.dens2edgenum(density,n)
    while keeptrying:
        G = gk.ringmore(n, plusedges)
        try:
            A = transitionMarix2(G, minstrength=st)
            keeptrying = False
        except ValueError:
            print "!!! Unable to find strong links for a stable matrix !!!"
            print "*** trying a different graph"
    return {'graph':      G,
            'transition': A,
            'converges':  len(bfutils.call_undersamples(G))}



# -------------------------------------------------------------------
# estimation
# -------------------------------------------------------------------

def scoreAGraph(G, data, x0 = None):
    A,B = npG2SVAR(G)
    K = scipy.sum(abs(A)+abs(B))
    a_idx = np.where(A != 0)
    b_idx = np.where(B != 0)
    if x0:
        o = optimize.fmin_bfgs(nllf, x0, args=(A, B, data, a_idx, b_idx),
                               disp=False, full_output=True)
    else:
        o = optimize.fmin_bfgs(nllf, scipy.randn(K),
                               args=(np.double(A), np.double(B),
                                     data, a_idx, b_idx),
                               disp=False, full_output=True)
    return 2*o(1) + K*np.log(T) #VARbic(o[1],K,data.shape[1])

def estimateG(G,YY,XX,YX,T,x0=None):
    A,B = npG2SVAR(G)
    K = scipy.sum(abs(A)+abs(B))
    a_idx = np.where(A != 0)
    b_idx = np.where(B != 0)
    try:
        s = x0.shape
        x = x0
    except AttributeError:
        x = scipy.randn(K)
    o = optimize.fmin_bfgs(nllf2, x,
                           args=(np.double(A), np.double(B),
                                 YY,XX,YX,T,a_idx, b_idx),
                           disp=False, full_output=True)
    A,B = x2M(o[0], np.double(A), np.double(B), a_idx, b_idx)
    return  A,B


def data2AB(data,x0=None):
    n = data.shape[0]
    T = data.shape[1]
    YY = np.dot(data[:,1:],data[:,1:].T)
    XX = np.dot(data[:,:-1],data[:,:-1].T)
    YX = np.dot(data[:,1:],data[:,:-1].T)    
    A = np.ones((n,n))
    B = np.ones((n,n))
    np.fill_diagonal(B,0)
    B[np.triu_indices(n)] = 0
    K = scipy.sum(abs(A)+abs(B))
    a_idx = np.where(A != 0)
    b_idx = np.where(B != 0)
    np.fill_diagonal(B,1)    
    try:
        s = x0.shape
        x = x0
    except AttributeError:
        x = scipy.randn(K)
    o = optimize.fmin_bfgs(nllf2, x,
                           args=(np.double(A), np.double(B),
                                 YY,XX,YX,T,a_idx, b_idx),
                           gtol=1e-12,
                           disp=False, full_output=True)
    A,B = x2M(o[0], np.double(A), np.double(B), a_idx, b_idx)
    B = B+B.T
    return  A,B

def AB2intAB(A,B, th=0.07):
    A[amap(lambda x: abs(x) > th, A)] = 1
    A[amap(lambda x: abs(x) < 1, A)] = 0
    B[amap(lambda x: abs(x) > th, B)] = 1
    B[amap(lambda x: np.abs(x) < 1, B)] = 0
    return A,B

def intAB2graph(A,B):
    n = A.shape[0]
    g = {str(i):{} for i in range(1,n+1)}

    for i in range(n):
        for j in range(n):
            if A[j,i]: g[str(i+1)][str(j+1)] = set([(0,1)])

    for i in range(n):
        for j in range(n):
            if B[j,i] and j!=i:
                if str(j+1) in g[str(i+1)]:
                    g[str(i+1)][str(j+1)].add((2,0))                    
                else:
                    g[str(i+1)][str(j+1)] = set([(2,0)])
    return g

def data2graph(data,x0=None):
    A,B = data2AB(data,x0=x0)
    Ab,Bb = AB2intAB(A,B)
    return intAB2graph(Ab,Bb)