#!/usr/bin/env python

import datetime
import gzip
import math
import numpy as np
import os
import sqlite3
import sys

DATA_DIR = '1000GP/'
GP_DIR = os.path.join(DATA_DIR, '1000GP_Phase3')
PREDICDB_FILE = 'weights.db'
BETASDB_FILE = 'betas.db'

class DB:
    def __init__(self, file):
        self.conn = sqlite3.connect(file)
    
    def __call__(self, query, args=None):
        c = self.conn.cursor()
        if args:
            for row in c.execute(query, args):
                yield row
        else:
            for row in c.execute(query):
                yield row

class PredicDB:
    def __init__(self):
        self.db = DB(PREDICDB_FILE)

    def SNP_iterator(self):
        for snp in self.db('SELECT DISTINCT rsid FROM weights ORDER BY rsid'):
            yield snp

    def gene_iterator(self):
        for snp in self.db('SELECT DISTINCT gene FROM weights ORDER BY rsid'):
            yield snp
predic_db = PredicDB()


class BetasDB:
    def __init__(self):
        self.db = DB(BETASDB_FILE)
beta_db = BetasDB()

def w(l, g):    
    "Weight of SNP l for gene g"
    ret = list(predic_db.db('SELECT weight FROM weights WHERE rsid=? AND gene=?', (l, g)))
    if len(ret) > 0:
        return ret[0][0]
    else:
        return 0

def beta(l):
    "GWAS regression coefficient for SNP l"
    ret = list(beta_db.db('SELECT beta FROM betas WHERE rsid=?', (l, )))
    if len(ret) > 0:
        return ret[0][0]
    else:
        return 0


if __name__ == '__main__':

    Zg = {} # gene -> z-score
    for g in db.gene_iterator():
        numerator = sum(w(l,g) * beta(l) * pow(sigma_l(l), 2) for l in db.SNP_iterator())
        denominator = sqrt(n * sum(w(l,g) * G(l, k) * w(k, g) for l in db.SNP_iterator() for k in db.SNP_iterator()))
        Zg[g] = numerator/denominator
        

    
