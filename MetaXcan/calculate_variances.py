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

class Iterator:
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
        return SNP(id, position, hap)


def calculate_variance(outfile, chrnumber, a):
    var = np.var(a.hap)
    outfile.write('chr%s\t%s\t%s\n'%(chrnumber, a.id, corr))

    

variances_file = open('variances.txt', 'w+')
for chrlegend in sorted(x for x in os.listdir(DIR) if '.legend.' in x):
    print "%s Processing %s..."%(datetime.datetime.now(), chrlegend)
    itr = Iterator(chrlegend)
    
    while True:
        nxt = itr.next_SNP()
        if nxt is None:
            break
        calculate_variance(variances_file, itr.chrnumber, nxt)

variances_file.close()
    
