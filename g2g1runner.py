import sys, os
sys.path.append('./tools/')
import traversal, bfutils
from multiprocessing import Pool,Process, Queue, cpu_count, current_process
import functools
import zickle as zkl
import time, socket
import scipy

INPNUM = 1 # number of randomized starts per graph
CAPSIZE= 1000 # stop traversing after growing equivalence class tothis size
REPEATS = 100
if socket.gethostname().split('.')[0] == 'leibnitz':
    PNUM=20
    PNUM=max((1,PNUM/INPNUM))
elif socket.gethostname().split('.')[0] == 'mars':
    PNUM=12
    PNUM=max((1,PNUM/INPNUM))
else:
    # Setting the number  of parallel running processes  to the number
    # of cores minus 7% for breathing room
    PNUM=cpu_count()-int(0.07*cpu_count())
    PNUM=max((1,PNUM/INPNUM))
print 'processes: ',PNUM, INPNUM

def wrapper(fold, n=10, k=10):
    scipy.random.seed()
    while True:
        try:
            g = bfutils.ringmore(n,k)
            g2 = traversal.increment_u(g,g)
            print fold,': ',traversal.density(g),':',
            startTime = int(round(time.time() * 1000))
            s = traversal.g22g1(g2, capsize=CAPSIZE)
            endTime = int(round(time.time() * 1000))
            print len(s)
        except MemoryError:
            print 'memory error... retrying'
            continue
        break
    return {'gt':g,'eq':s,'ms':endTime-startTime}

def killall(l):
    for e in l: 
        e.join(timeout=0.001)
        if not e.is_alive():
            #print 'first result'            
            for p in l:
                if p != e:
                    #print 'terminating ', p.name
                    p.terminate()
                    p.join()
                else:
                    p.join()
            return True
    return False

def fan_wrapper(fold,n=10,k=10):
    scipy.random.seed()
    curr_proc=current_process()
    curr_proc.daemon=False
    output = Queue()
    while True:
        try:
            g = bfutils.ringmore(n,k)
            g2 = traversal.increment_u(g,g)
            #g2 = bfutils.undersample(g,2)
            def inside_wrapper():
                scipy.random.seed()
                print fold,': ',traversal.density(g),':',
                startTime = int(round(time.time() * 1000))
                s = traversal.v2g22g1(g2, capsize=CAPSIZE)
                endTime = int(round(time.time() * 1000))
                print len(s),'   ',
                print (endTime-startTime)/1000.,' seconds'
                output.put({'gt':g,'eq':s,'ms':endTime-startTime})
            pl = [Process(target=inside_wrapper) for x in range(INPNUM)]
            for e in pl: e.start()
            while True:
                if killall(pl): break
            r = output.get()           
        except MemoryError:
            print 'memory error... retrying'
            for p in pl:
		p.terminate()
		p.join()
            continue
        break
    for p in pl: p.join()
    return r

densities = {6: [0.2, 0.25, ],
             8: [0.15, 0.2, 0.25,],
             10:[0.3],#0.15, 0.2, 0.25, 0.3],
             15:[0.25, 0.30],#0.1, 0.15, 0.2, 0.25, 0.3],
             20:[0.15, 0.2, 0.25, 0.30]}#0.1, 0.15, 0.2, 0.25, 0.3]}

for nodes in [15]:
#for nodes in [8]:
    z = {}
    pool=Pool(processes=PNUM)
    for dens in densities[nodes]:
    #for dens in [0.4, 0.5, 0.6, 0.7, 0.8]:
    #for dens in [0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]:
        e = bfutils.dens2edgenum(dens, n=nodes)
        eqclasses = pool.map(functools.partial(fan_wrapper, n=nodes, k=e), 
                             range(REPEATS))
        z[dens] = eqclasses
        zkl.save(z[dens],
                 socket.gethostname().split('.')[0]+\
                     '_nodes_'+str(nodes)+'_density_'+str(dens)+'_U2.zkl')
        print ''
        print '----'
        print ''
    pool.close()
    pool.join()
    zkl.save(z,socket.gethostname().split('.')[0]+'_nodes_'+str(nodes)+'_U2.zkl')
