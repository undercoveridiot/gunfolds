import sys, os
sys.path.append('./tools/')
import traversal, bfutils
from multiprocessing import Pool,Process, Queue, cpu_count
import functools
import zickle as zkl
import time, socket
import scipy

if socket.gethostname().split('.')[0] == 'leibnitz':
    PNUM=45
else:
    # Setting the number  of parallel running processes  to the number
    # of cores minus 7% for breathing room
    PNUM=cpu_count()-int(0.07*cpu_count()) 
print 'processes: ',PNUM
repeats = 100

def wrapper(fold, n=10, k=10):
    scipy.random.seed()
    while True:
        try:
            g = bfutils.ringmore(n,k)
            g2 = traversal.increment_u(g,g)
            print fold,': ',traversal.density(g),':',
            startTime = int(round(time.time() * 1000))
            s = traversal.eqc(g2, capsize=10000)
            endTime = int(round(time.time() * 1000))
            print len(s)
        except MemoryError:
            print 'memory error... retrying'
            continue
        break
    return {'gt':g,'eq':s,'ms':endTime-startTime}

def wrapper2(g):
    scipy.random.seed()
    g2 = traversal.increment_u(g,g)
    print ': ',traversal.density(g),':',
    startTime = int(round(time.time() * 1000))
    s = traversal.eqc(g2, capsize=10000)
    endTime = int(round(time.time() * 1000))
    print len(s)
    output.put({'gt':g,'eq':s,'ms':endTime-startTime})

def killall(l):
    for e in l: 
        e.join(timeout=0.001)
        if not e.is_alive():
            print 'first result'            
            for p in l:
                if p != e:
                    print 'terminating ', p.name
                    p.terminate()
                    p.join()
            return True
    return False
    
output = Queue()
g = bfutils.ringmore(6,12)
pl = [Process(target=wrapper2, args=(g,)) for x in range(5)]
for e in pl: e.start()
while True:
    if killall(pl): break
print output.get()['eq']
stop

#for nodes in [10, 15, 20, 30, 60]:
for nodes in [6]:
    z = {}
    pool=Pool(processes=PNUM)
    for dens in [0.2, 0.23, 0.28, 0.32, 0.35]:
        e = bfutils.dens2edgenum(dens, n=nodes)
        eqclasses = pool.map(functools.partial(wrapper, n=nodes, k=e), 
                             range(repeats))
        z[dens] = eqclasses
    pool.close()
    pool.join()
    zkl.save(z,socket.gethostname().split('.')[0]+'_nodes_'+str(nodes)+'.zkl')
