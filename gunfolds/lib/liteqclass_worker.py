from collections import namedtuple
import gunfolds.tools.bfutils as bfu
from gunfolds.tools.conversions import num2CG
import gunfolds.tools.unknownrate
import sys

Message = namedtuple('Message', ['gnum', 'sloops'])

class LitEqClassWorker(object):

    def __init__(self, H, cp, ccf, capsize, in_queue, out_queue):
        self.H = H
        self.n = len(H)
        self.cp = cp
        self.ccf = ccf
        self.capsize = capsize
        self.in_queue = in_queue
        self.out_queue = out_queue
        self.dsr = {}
        self.s = set()
        self.ss = set()
        self.gset = set()
        self.eset = set()
        self.currsize = 0

    def run(self):

        while True:
            try:
                msg = self.in_queue.get() # blocks here until data in q
            except Exception, e:
                sys.exit('Error reading off of in_queue: {}'.format(e))

                for sloop in msg.sloops:
                    if sloop & msg.gnum == sloop:
                        continue
                    num = sloop | msg.gnum
                    if sloop in self.ccf and gunfolds.tools.unknownrate.skip_conflictors(num, self.ccf[sloop]):
                        continue
                    if not num in self.s:
                        g = num2CG(num, self.n)
                        if not bfu.call_u_conflicts(g, self.H):
                            self.s.add(num)
                            self.gset.add((num, sloop))
                            self.eset.add(sloop)
                            if bfu.call_u_equals(g, self.H):
                                self.ss.add(num)
                                if self.capsize <= len(self.ss) + currsize:
                                    self.out_queue.put((self.dsr, self.ss))
                                    self.in_queue.task_done()

                for gn, e in self.gset:
                    if e in self.cp:
                        self.dsr[gn] = self.eset - self.cp[e] - {e}
                    else:
                        self.dsr[gn] = self.eset - {e}

                self.out_queue.put((self.dsr, self.ss))
                self.in_queue.task_done()



    def reset(self, cursize):
        # Reset for next round
        self.dsr = {}
        self.s = set()
        self.ss = set()
        self.gset = set()
        self.eset = set()
        self.cursize = cursize
