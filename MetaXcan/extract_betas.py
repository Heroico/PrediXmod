#!/usr/bin/env python

import gzip
import os
import sqlite3
import sys

DATA_DIR = 'im-lab'
TARGET_FILE = 'betas.db'

if os.path.exists(TARGET_FILE):
    os.unlink(TARGET_FILE)

class DB:
    def __init__(self, file):
        self.conn = sqlite3.connect(file)
    
    def __call__(self, query, args=None):
        c = self.conn.cursor()
        if args:
            ret = list(c.execute(query, args))
        else:
            ret = list(c.execute(query))
        self.conn.commit()
        return ret        

    def close(self):
        self.conn.close()

if __name__ == '__main__':

    os.system("""zcat `find %s -iname "*dosage.gz"` | cut -f2,4 > betas.txt"""%DATA_DIR)
    with open('load.sql', 'w+') as outfile:
        outfile.write('CREATE TABLE betas (rsid TEXT, beta real);\n')
        outfile.write('CREATE INDEX beta_index ON betas (rsid);\n')
        outfile.write('.separator "\\t"\n')
        outfile.write('.import betas.txt betas\n')
    os.system('sqlite3 betas.db < load.sql')
            
    
