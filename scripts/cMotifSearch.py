#!/usr/bin/python
#0

""" Module/Script Description

Copyright (c) 2008 Xiaowei Chen <bighanchen2008@gmail.com>
This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).
@status:  experimental
@version: $Revision$
@author:  Xiaowei Chen
@contact: bighanchen2008@gmail.com
"""

# ------------------------------------
# python modules
# ------------------------------------

import sys
import string
import os
#from wbed import ColumnReader

# ------------------------------------
# constants
# ------------------------------------
MOTIF=["GGACA","GGACC","GGACT","GAACA","GAACC","GAACT","AGACA","AGACC","AGACT","AAACA","AAACC","AAACT"]
# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main
# ------------------------------------

if __name__=="__main__":
    if len(sys.argv)==1:
        print "Usage: "+sys.argv[0]+" .bed .fa > motif.bed"
    else:
        chrom=[]
        start=[]
        end=[]
        ID=[]
        score=[]
        strand=[]
        fh=open(sys.argv[1])
        for line in fh:
            line=line.rstrip(os.linesep)
            line=line.split("\t")
            chrom.append(line[0])
            start.append(line[1])
            end.append(line[2])
            ID.append(line[3])
            score.append(line[4])
            strand.append(line[5])
        ft=open(sys.argv[2])
        count=1
        for line in ft:
            line=line.rstrip(os.linesep)
            if line[0]==">":
                seqID=line.split(">")[1]
                idx=ID.index(seqID)
            else:
                for i in range(len(line)-5+1):
                    if line[i:i+5] in MOTIF:
                        print chrom[idx]+"\t"+str(string.atoi(start[idx])+i)+"\t"+str(string.atoi(start[idx])+i+5)+"\tM6A"+str(count).zfill(6)+"\t"+score[idx]+"\t"+strand[idx]+"\t"+seqID
                        count=count+1
