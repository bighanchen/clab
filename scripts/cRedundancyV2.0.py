#!/usr/bin/python
#Last-modified: 14 Dec 2011 10:08:37 PM

""" Module/Scripts Description

Copyright (c) 2008 Xiaowei Chen <bighanchen2008@gmail.com>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  experimental
@version: $Revision$
@author:  Xiaowei Chen
@contact: bighanchen2008@gmail.com

Thank Liqing Tian for the improved algorithm in this version. 

"""

# ------------------------------------
# python modules
# ------------------------------------

from __future__ import division
import sys
import string

# ------------------------------------
# constants
# ------------------------------------

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
        print "Usage: "+sys.argv[0]+" sequence_ID.list sorted_blatoutput.psl seqLenRatio similarity > noRedundancy.list"
        print "sorted_blatoutput.psl is the sorted output of blat in psl format. This programme can find sequences with similarity >=95%(default)"
    else:
        try:
            seqLenRatio=sys.argv[3]
            if (seqLenRatio>1) or (seqLenRatio<=0)
                seqLenRatio=0.95
        except:
            seqLenRatio=0.95
        try:
            similarity=sys.argv[4]
            if (similarity>1) or (similarity<=0)
                similarity=0.95
        except:
            similarity=0.95
        fSeqID=open(sys.argv[1])
        seq_ID=fSeqID.readline().rstrip("\n\r")
        seq_ID_list=[]
        while seq_ID:
            seq_ID_list.append(seq_ID)
            seq_ID=fSeqID.readline().rstrip("\n\r")
        fBlatOutput=open(sys.argv[2])
        while 
        




        identicalNc={}
        queue=[]
        keyqueue=[]
        valuequeue=[]
        fh=open(sys.argv[1])
        line=fh.readline().rstrip("\n\r")
        while line:
            line=line.split("\t")
            n0=string.atoi(line[0])
            n10=string.atoi(line[10])
            n14=string.atoi(line[14])
            if (n10>=n14):
                longer=n10
                shorter=n14
            else:
                longer=n14
                shorter=n10
            if(shorter/longer>=0.95) and (n0/shorter>=0.95):
                if (line[9] not in queue) and (line[13] not in queue):
                    identicalNc[line[9]]=line[13]
                    keyqueue.append(line[9])
                    valuequeue.append(line[13])
                    queue.append(line[9])
                    queue.append(line[13])
                elif (line[9] in queue) and (line[9] in identicalNc) and (line[13] not in queue):
                    identicalNc[line[9]]=identicalNc[line[9]]+":"+line[13]
                    keyqueue.append(line[9])
                    valuequeue.append(line[13])
                    queue.append(line[13])
                elif (line[9] in queue) and (line[9] in identicalNc) and (line[13] in queue):
                    pass
                elif (line[9] in queue) and (line[9] not in identicalNc) and (line[13] in queue):
                    pass
                elif (line[9] in queue) and (line[9] not in identicalNc) and (line[13] not in queue):
                    idx=valuequeue.index(line[9])
                    identicalNc[keyqueue[idx]]=identicalNc[keyqueue[idx]]+":"+line[13]
                    keyqueue.append(keyqueue[idx])
                    valuequeue.append(line[13])                    
                    queue.append(line[13])
                elif (line[9] not in queue) and (line[13] in identicalNc):
                    identicalNc[line[13]]=identicalNc[line[13]]+":"+line[9]
                    keyqueue.append(line[13])
                    valuequeue.append(line[9])
                    queue.append(line[9])
                elif (line[9] not in queue) and (line[13] not in identicalNc) and (line[13] in queue):
                    idx=valuequeue.index(line[13])
                    identicalNc[keyqueue[idx]]=identicalNc[keyqueue[idx]]+":"+line[9]
                    keyqueue.append(keyqueue[idx])
                    valuequeue.append(line[9])                    
                    queue.append(line[9])
            line=fh.readline().rstrip("\n\r")
        for key in identicalNc.keys():
            print key+"\t"+identicalNc[key]

            #line=line.rstrip("\n\r").split("\t")
