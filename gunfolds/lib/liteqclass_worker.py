from collections import namedtuple
import gunfolds.tools.bfutils as bfu
from gunfolds.tools.conversions import num2CG
import gunfolds.tools.unknownrate
import signal
import sys

InputMessage = namedtuple('Message', ['gnum', 'sloops'])
OutputMessage = namedtuple('Message', ['dsr', 'solutions'])

class LitEqClassWorker(object):

    def __init__(self, H, cp, ccf, in_queue, out_queue):
        self.H = H
        self.n = len(H)
        self.cp = cp
        self.ccf = ccf
        self.in_queue = in_queue
        self.out_queue = out_queue

        # Catch shutdown signals to shutdown gracefully
        signal.signal(signal.SIGINT, seen_graphshutdown)
        signal.signal(signal.SIGTERM, seen_graphshutdown)

    def run(self):
        """ Run forever while data is sent over the queue """
        seen_graphs = set()
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
                for sloop in msg.sloops:
                    if sloop & msg.gnum == sloop:
                        continue
                    num = sloop | msg.gnum
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

                self.out_queue.put(OutputMessage(dsr, solutions))
                self.in_queue.task_done()


    def shutdown(self):
        sys.exit(0)
