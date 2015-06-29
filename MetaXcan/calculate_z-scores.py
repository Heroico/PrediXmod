#!/usr/bin/env python

import datetime
import json
from math import sqrt
import numpy as np
import os
import sqlite3
import sys

DATA_DIR = '1000GP/'
GP_DIR = os.path.join(DATA_DIR, '1000GP_Phase3')
PREDICDB_FILE = 'weights.db'
BETASDB_FILE = 'betas.db'
VARIANCESDB_FILE = 'variances.db'
CORRELATIONSDB_FILE = 'correlations.db'

def get_n():
    for line in open(os.path.join(GP_DIR, [x for x in os.listdir(GP_DIR) if '.hap.' in x][0])):
        return len(line.strip().split())
    

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

    def __call__(self, query, args=None):
        return self.db(query, args)

    def SNP_iterator(self, gene):
        for snp in self('SELECT DISTINCT rsid FROM weights WHERE gene=? ORDER BY rsid', (gene,)):
            yield snp[0]

    def gene_iterator(self):
        for snp in self('SELECT DISTINCT gene FROM weights ORDER BY rsid'):
            yield snp[0]
predic_db = PredicDB()


beta_db = DB(BETASDB_FILE)
variances_db = DB(VARIANCESDB_FILE)
correlations_db = DB(CORRELATIONSDB_FILE)

def w(l, g):    
    "Weight of SNP l for gene g"
    ret = list(predic_db.db('SELECT weight FROM weights WHERE rsid=? AND gene=?', (l, g)))
    if len(ret) > 0:
        return ret[0][0]
    else:
        return 0

def weights_vector(g):
    "Vector of weights and SNPs for geneg"
    res = list(predic_db.db('SELECT rsid, weight FROM weights WHERE gene=?', (g,)))
    SNPs = [x[0] for x in res]
    weights = [x[1] for x in res]
    return SNPs, weights


def beta(l):
    "GWAS regression coefficient for SNP l"
    ret = list(beta_db('SELECT beta FROM betas WHERE rsid=?', (l, )))
    if len(ret) > 0:
        return ret[0][0]
    else:
        return 0

def sigma_l(l):
    "Variance for SNP l"
    ret = list(variances_db('SELECT variance FROM variances WHERE rsid=?', (l, )))
    if len(ret) > 0:
        return ret[0][0]
    else:
        return 0

def G(l, k):
    "Correlation between SNPs l and k"
    ret = list(correlations_db('SELECT correlation FROM correlations WHERE rsid_a=? AND rsid_B=?', (l, k)))
    if len(ret) > 0:
        return ret[0][0]
    else:
        return 0

def correlated_with(l):
    "SNPs correlated with l, and their correlation value"
    yield (l, 1)
    for (k, c) in correlations_db("SELECT rsid_B, correlation FROM correlations WHERE rsid_a=?", (l,)):
        yield (k, c)

def numerator_term(l, g):
    x = w(l, g)
    if x == 0:
        return 0
    y = beta(l)
    if y == 0:
        return 0
    return pow(sigma_l(l), 2)

def calculate_numerator(SNPs, g):
    return sum(numerator_term(l, g) for l in SNPs)

def calculate_denominator_old(n,g):
    acc = 0.0
    valid_SNPs = set()
    for l in predic_db.SNP_iterator(gene=g):
        for k, c in correlated_with(l):
            acc += w(l,g) * c *  w(k, g)
            valid_SNPs.update((l, k))
    return sqrt(n * acc), valid_SNPs

def calculate_denominator(n,g):
    valid_SNPs, weights = weights_vector(g)
    c_snp = dict(zip(valid_SNPs, range(len(valid_SNPs))))
    weights = np.array(weights)
    G = np.zeros((len(weights), len(weights)))
    for x, l in enumerate(valid_SNPs):
        G[x,x] = 1
        for k, c in correlated_with(l):
            if k in c_snp:
                y = c_snp[k]
                G[x,y] = c
    return np.dot(np.dot(np.transpose(weights), G), weights), valid_SNPs

def calculate_score(n, g):
    denominator, valid_SNPs = calculate_denominator(n, g)
    numerator = calculate_numerator(valid_SNPs, g)
    return numerator/denominator if denominator > 0 else 'denominator_0'


if __name__ == '__main__':
    Zg = {} # gene -> z-score * sigma_y
    n = get_n()

    for k, g in enumerate(predic_db.gene_iterator()):
        if k%100 == 0:
            print "%d out of 11239"%k, datetime.datetime.now()
        Zg[g] = calculate_score(n, g)

    with open('zg.json', 'w+') as outfile:
        outfile.write(json.dumps(Zg))

    
