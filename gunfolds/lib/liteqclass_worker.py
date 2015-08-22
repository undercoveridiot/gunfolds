import gunfolds.tools.bfutils as bfu
from gunfolds.tools.conversions import num2CG
import gunfolds.tools.unknownrate
from multiprocessing import Process
import signal
import sys

class LitEqClassWorker(Process):

    def __init__(self, H, cp, ccf, in_queue, out_queue, *args, **kwargs):
        super(LitEqClassWorker, self).__init__(*args, **kwargs)
        self.H = H
        self.n = len(H)
        self.cp = cp
        self.ccf = ccf
        self.in_queue = in_queue
        self.out_queue = out_queue

        # Catch shutdown signals to shutdown gracefully
        signal.signal(signal.SIGINT, self.shutdown)
        signal.signal(signal.SIGTERM, self.shutdown)

    def run(self):
        """ Run forever while data is sent over the queue """
        seen_graphs = set() # Cache pre process
        while True:
            dsr = {}
            solutions = set()
            gset = set()
            eset = set()
            try:
                msg = self.in_queue.get() # blocks here until data in q
            except Exception, e:
                sys.exit('Error reading off of in_queue: {}'.format(e))
            else:
                gnum, sloops = msg
                for sloop in sloops:
                    if sloop & gnum == sloop:
                        continue
                    num = sloop | gnum
                    if sloop in self.ccf and gunfolds.tools.unknownrate.skip_conflictors(num, self.ccf[sloop]):
                        continue
                    if not num in seen_graphs:
                        g = num2CG(num, self.n)
                        if not bfu.call_u_conflicts(g, self.H):
                            seen_graphs.add(num)
                            gset.add((num, sloop))
                            eset.add(sloop)
                            if bfu.call_u_equals(g, self.H):
                                solutions.add(num)

                for gn, e in gset:
                    if e in self.cp:
                        dsr[gn] = eset - self.cp[e] - {e}
                    else:
                        dsr[gn] = eset - {e}

                self.out_queue.put((dsr, solutions))


    def shutdown(self, signalNum=0, frame=0):
        sys.exit(0)
