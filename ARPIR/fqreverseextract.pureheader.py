#!/usr/bin/python
# -*- coding: utf-8 -*-
# da usare quando vuoi estrarre FASTQ avendo i dati dal BED.

import sys 
# acquire ID from file as list, adding the prefix @
ids, lineno, flag = set('@' + x.strip() for x in open(sys.argv[1])), 0, False
#print ids
for l in sys.stdin:
    lineno += 1
    if lineno%4 == 1: flag = (l.split(' ')[0].strip() not in ids)
    if flag: print l,

