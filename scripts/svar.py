import signal
import time, socket
import numpy as np
import scipy
import functools, itertools
import progressbar as pb
import sys,os
TOOLSPATH='~/soft/src/dev/craft/gunfolds/tools/'
sys.path.append(os.path.expanduser(TOOLSPATH))

from multiprocessing import Pool,Process, Queue, cpu_count, current_process
import linear_model as lm
import traversal as trv
import bfutils as bfu
import graphkit as gk
import zickle as zkl

NOISE_STD = '1.'
DEPTH=2
PARALLEL=True
INPNUM = 1 # number of randomized starts per graph
CAPSIZE= 1000 # stop traversing after growing equivalence class tothis size
REPEATS = 100
if socket.gethostname().split('.')[0] == 'leibnitz':
    PNUM=60
    PNUM=max((1,PNUM/INPNUM))
elif socket.gethostname().split('.')[0] == 'mars':
    PNUM=12
    PNUM=max((1,PNUM/INPNUM))
elif socket.gethostname().split('.')[0] == 'hooke':
    PNUM=21
    PNUM=max((1,PNUM/INPNUM))
else:
    # Setting the number  of parallel running processes  to the number
    # of cores minus 7% for breathing room
    PNUM=cpu_count()-int(0.07*cpu_count())
    PNUM=max((1,PNUM/INPNUM))
print 'processes: ',PNUM, INPNUM

def timeout(func, args=(), kwargs={}, timeout_duration=1, default=None):
    import signal

    class TimeoutError(Exception):
        pass

    def handler(signum, frame):
        raise TimeoutError()

    # set the timeout handler
    signal.signal(signal.SIGALRM, handler)
    signal.alarm(timeout_duration)
    try:
        result = func(*args, **kwargs)
    except TimeoutError as exc:
        result = default
    finally:
        signal.alarm(0)

    return result

def examine_bidirected_flips(g2, depth=0):
    s = trv.v2g22g1(g2, capsize=CAPSIZE, verbose=False)
    if s: return s
    bedges = gk.bedgelist(g2)
    all_bedges = [x for x in itertools.combinations(g2,2)]
    for be in all_bedges:
        if be in bedges:
            g2[be[0]][be[1]].remove((2,0))
            g2[be[1]][be[0]].remove((2,0))
            s = trv.v2g22g1(g2, capsize=CAPSIZE, verbose=False)
            if s:
                return s
            else:
                if depth: examine_bidirected_flips(g2, depth=depth-1)
            if s: return s                
            g2[be[0]][be[1]].add((2,0))
            g2[be[1]][be[0]].add((2,0))
        else:
            mask = [be[1] in g2[be[0]], be[0] in g2[be[1]]]
            if mask[0]:
                g2[be[0]][be[1]].add((2,0))
            else:
                g2[be[0]][be[1]] = set([(2,0)])
            if mask[1]:
                g2[be[1]][be[0]].add((2,0))
            else:
                g2[be[1]][be[0]] = set([(2,0)])
            s = trv.v2g22g1(g2, capsize=CAPSIZE, verbose=False)
            if s:
                return s
            else:
                if depth:
                    s = examine_bidirected_flips(g2, depth=depth-1)
            if s: return s
            g2[be[0]][be[1]].remove((2,0))
            g2[be[1]][be[0]].remove((2,0))
    return s

def wrapper(fold,n=10,dens=0.1):
    scipy.random.seed()
    curr_proc=current_process()
    curr_proc.daemon=False
    output = Queue()
    rate = 2

    r = None
    s = set()
    counter = 0
    while not s:
        scipy.random.seed()        
        print 'o',
        sys.stdout.flush()
        sst = 0.5
        r = None        
        while not r:
            r = timeout(lm.getAring, args=(n, dens, sst, False),
                        timeout_duration=60)
            sst -= 0.05
            if sst < 0: break
        g = r['graph']
        true_g2 = bfu.undersample(g, rate-1)        
        data = lm.drawsamplesLG(r['transition'], samples=20000,
                                nstd=np.double(NOISE_STD))
        data = data[:,10000:]
        startTime = int(round(time.time() * 1000))
        g2 = lm.data2graph(data[:,::2])
        print gk.OCE(g2,true_g2)
        #s = examine_bidirected_flips(g2, depth=DEPTH)
        #s = trv.v2g22g1(g2, capsize=CAPSIZE, verbose=False)
        s = trv.edge_backtrack2g1_directed(g2, capsize=CAPSIZE)
        if -1 in s: s=set()
        endTime = int(round(time.time() * 1000))
        counter += 1
    print ''
    oce = [gk.OCE(bfu.num2CG(x,n),g) for x in s]
    cum_oce = [sum(x['directed'])+sum(x['bidirected']) for x in oce]
    idx = np.argmin(cum_oce)
    print "{:2}: {:8} : {:4}  {:10} seconds".\
          format(fold, round(dens,3), cum_oce[idx],
                 round((endTime-startTime)/1000.,3))
    return {'gt':r,
            'eq':s,
            'OCE':oce[idx],
            'estimate': g2,
            'graphs_tried': counter,
            'strength':sst+0.01,
            'ms':endTime-startTime}


def wrapgen(fold,n=10,dens=0.1):
    scipy.random.seed()
    curr_proc=current_process()
    curr_proc.daemon=False
    output = Queue()
    rate = 2

    r = None
    s = set()
    counter = 0
    sst = 0.5
    r = None        
    while not r:
        r = timeout(lm.getAring, args=(n, dens, sst, False),
                    timeout_duration=60)
        sst -= 0.01
        if sst < 0: break
    return r

densities = {6: [0.25, 0.3, 0.35],
             8: [0.15, 0.2, 0.25, 0.3],
             10:[0.15, 0.25, 0.3],
             15:[0.1, 0.15, 0.2],
             20:[0.1],
             25:[0.1],
             30:[0.1],
             35:[0.1]}

wrp = wrapgen

for nodes in [10]:
    z = {}
    pool=Pool(processes=PNUM)
    for dens in densities[nodes]:
        print "{:2}: {:8} : {:10}  {:10}".format('id', 'density', 'OCE', 'time')

        if PARALLEL:
            errors = pool.map(functools.partial(wrp, n=nodes,
                                                dens=dens),
                              range(REPEATS))
        else:
            errors = []
            for i in range(REPEATS):
                errors.append(wrp(i,n=nodes,dens=dens))
        print 'computed'
        z[dens] = errors
        zkl.save(z[dens],
                 socket.gethostname().split('.')[0]+\
                     '_nodes_'+str(nodes)+'_density_'+str(dens)+'_OCE_model_.zkl')                 
#                     '_nodes_'+str(nodes)+'_density_'+str(dens)+'_noise_'+NOISE_STD+'_OCE_model.zkl')
        print ''
        print '----'
        print ''
    pool.close()
    pool.join()
#    zkl.save(z,socket.gethostname().split('.')[0]+'_nodes_'+str(nodes)+'_noise_'+NOISE_STD+'_OCE_model.zkl'
    zkl.save(z,socket.gethostname().split('.')[0]+'_nodes_'+str(nodes)+'_OCE_model_.zkl')
