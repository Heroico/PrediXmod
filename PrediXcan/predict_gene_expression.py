#!/usr/bin/env python

import argparse
from collections import defaultdict
import datetime
import gzip
import numpy as np
import os
import sqlite3
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--genelist', action="store", dest="genelist", default=None, help="Text file with chromosome, gene pairs.")
parser.add_argument('--dosages', action="store", dest="dosages", default="data/dosages", help="Path to a directory of gzipped dosage files.")
parser.add_argument('--dosages_prefix', dest="dosages_prefix", default="chr", action="store", help="Prefix of filenames of gzipped dosage files.")
parser.add_argument('--weights', action="store", dest="weights",default="data/weights.db",  help="SQLite database with rsid weights.")
parser.add_argument('--output', action="store", dest="output", default="output.txt", help="Path to the output file.")
args = parser.parse_args()

 
GENE_LIST = args.genelist
DOSAGE_DIR = args.dosages
DOSAGE_PREFIX = args.dosages_prefix
BETA_FILE = args.weights
OUTPUT_FILE = args.output


def get_all_dosages():
    for chrfile in [x for x in sorted(os.listdir(DOSAGE_DIR)) if x.startswith(DOSAGE_PREFIX)]:
        print datetime.datetime.now(), "Processing %s"%chrfile
        for line in gzip.open(os.path.join(DOSAGE_DIR, chrfile)):
            arr = line.strip().split()
            rsid = arr[1]
            refallele = arr[4]
            dosage_row = np.array(map(float, arr[6:]))
            yield rsid, refallele, dosage_row

class WeightsDB:
    def __init__(self):
        self.conn = sqlite3.connect(BETA_FILE)

    def query(self, sql, args=None):
        c = self.conn.cursor()
        if args:
            for ret in c.execute(sql, args):
                yield ret
        else:
            for ret in c.execute(sql):
                yield ret

class GetApplicationsOf:
    def __init__(self):
        self.db = WeightsDB()

    def __call__(self, rsid):
        for tup in self.db.query("SELECT gene, weight, eff_allele FROM weights WHERE rsid=?", (rsid,)):
            yield tup
get_applications_of = GetApplicationsOf()

class TranscriptionMatrix:
    def __init__(self):
        self.D = None

    def get_gene_list(self):
        if GENE_LIST:
            return list(sorted([line.strip().split()[-1] for line in open(GENE_LIST)]))
        else:
            return [tup[0] for tup in WeightsDB().query("SELECT DISTINCT gene FROM weights ORDER BY gene")]

    def update(self, gene, weight, ref_allele, allele, dosage_row):
        if self.D is None:
            self.gene_list = self.get_gene_list()  
            self.gene_index = { gene:k for (k, gene) in enumerate(self.gene_list) }
            self.D = np.zeros((len(self.gene_list), len(dosage_row))) # Genes x Cases
        if gene in self.gene_index:            
            multiplier = 1 if ref_allele == allele else -1
            self.D[self.gene_index[gene],] += dosage_row * weight * multiplier # Update all cases for that gene

    def save(self):
        with open(OUTPUT_FILE, 'w+') as outfile:
            outfile.write('\t'.join(self.gene_list) + '\n') # Nb. this lists the names of rows, not of columns
            for col in range(0, self.D.shape[1]):
                outfile.write('\t'.join(map(str, self.D[:,col]))+'\n')

transcription_matrix = TranscriptionMatrix()
for rsid, allele, dosage_row in get_all_dosages():
    for gene, weight, ref_allele in get_applications_of(rsid):
        transcription_matrix.update(gene, weight, ref_allele, allele, dosage_row)
transcription_matrix.save()        
