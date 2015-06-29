#!/usr/bin/env python

import gzip
import os
import sqlite3
import sys

DATA_DIR = 'im-lab'

if __name__ == '__main__':
    os.system("""zcat `find %s -iname "*dosage.gz"` | cut -f2,4 > betas.txt"""%DATA_DIR)
    with open('load.sql', 'w+') as outfile:
        outfile.write('CREATE TABLE betas (rsid TEXT, beta real);\n')
        outfile.write('CREATE INDEX beta_index ON betas (rsid);\n')
        outfile.write('.separator "\\t"\n')
        outfile.write('.import betas.txt betas\n')
    os.system('sqlite3 betas.db < load.sql')
            
    
