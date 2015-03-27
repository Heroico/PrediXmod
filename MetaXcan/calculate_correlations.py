#!/usr/bin/env python

import datetime
import gc
import gzip
import numpy as np
import os
import sys

DIR = '1000GP/1000GP_Phase3/'
SNP_LIST = set(x.strip() for x in open('SNP_LIST'))

class SNP:
    def __init__(self, id, position, hap):
        self.id = id
        self.position = position
        self.hap = hap

class DoubleIterator:
    def __init__(self, chrlegend):
        self.chrlegend = chrlegend
        self.chrhap = chrlegend.replace('legend', 'hap')
        self.chrnumber = chrlegend.split('.')[0].split('chr')[1]

        self.legend_file = gzip.open(os.path.join(DIR, self.chrlegend))
        self.legend_file.readline() # ignore the header
        self.hap_file = gzip.open(os.path.join(DIR, self.chrhap))

        self.buffer = []
        self.ret = None

    def next_SNP(self):
        if len(self.buffer) == 0:
            while True:
                legend_line = self.legend_file.readline()
                hap_line = self.hap_file.readline()
                if legend_line == '':
                    return None
                elif legend_line.split(':')[0] in SNP_LIST:
                    break
            id = legend_line.split(':')[0]
            position = int(legend_line.split(':')[1])
            hap = list(map(int, hap_line.strip().split()))
            self.buffer.append(SNP(id, position, hap))
        self.ret = self.buffer.pop()
        return self.ret

    def rewind(self):
        assert(self.ret is not None)
        self.buffer.append(self.ret)
        self.ret = None
            

def calculate_correlations(outfile, chrnumber, a, b):
    corr = np.corrcoef(a.hap, b.hap)[0,1]
    outfile.write('chr%s\t%s\t%s\t%s\n'%(chrnumber, a.id, b.id, corr))

    

correlations_file = gzip.open('correlations.txt.gz', 'w+')
for chrlegend in sorted(x for x in os.listdir(DIR) if '.legend.' in x):
    print "%s Processing %s..."%(datetime.datetime.now(), chrlegend)
    itr = DoubleIterator(chrlegend)
    window = []
    
    # Main cycle
    while True:

        # Fill the window
        while True:
            nxt = itr.next_SNP()
            if nxt is None:
                itr.rewind() ; break
            elif len(window) == 0 or (nxt.position - window[0].position) < 1000000:
                if (len(window) == 0) or (nxt.id != window[-1].id):
                    window.append(nxt)
            else:
                itr.rewind() ; break

        #print "%s window size = %d"%(datetime.datetime.now(), len(window))
        # Calculate all correlations between the leftmost SNP and the rest of the window
        if len(window) > 1:
            for idx in range(1, len(window)):
                calculate_correlations(correlations_file, itr.chrnumber, window[0], window[idx])

        # Remove the leftmost SNP from the window
        del window[0]
        gc.collect() # This could be done much less often, but it doesn't seem to add much to runtime.
        if len(window) == 0:
            break

correlations_file.close()
    
