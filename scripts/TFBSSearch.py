#!/usr/bin/python
#0

""" Script Description

This Script can help to search TFBS in lncRNA gene promoters.

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
import re
import math
import numpy
import time
import random
import itertools
#from wbed import ColumnReader

# ------------------------------------
# constants
# ------------------------------------

# The frequency of A,C,G,T in hg19.
HG19BK=[0.2950886,0.2046611,0.2047848,0.2954655]
DEFAULT=math.log(1.0/101)
BASE=['A','C','G','U']


# ------------------------------------
# Misc functions
# ------------------------------------
def GetMaxScore(Matrix,MotifLen):
    #Initialize
    #scorestr="" #String of scores for print
    #score=[0]*(len(PromoterSeq)-MotifLen+1) #List to store normalized scores
    #Len=len(PromoterSeq)-MotifLen+1
    #curscore=[0]*Len#List to store original scores
    maxscore=0#store max score
    for j in range(MotifLen):
        maxscore+=max(Matrix[j])
    return maxscore

def GetScore(Matrix,PromoterSeq,MotifLen):
    #Initialize
    #scorestr="" #String of scores for print
    #score=[0]*(len(PromoterSeq)-MotifLen+1) #List to store normalized scores
    Len=len(PromoterSeq)-MotifLen+1
    curscore=[0]*Len#List to store original scores
    for i in range(Len):
        flag=1
        for k in range(i,i+MotifLen):
            if PromoterSeq[k] not in BASE:
                flag=0
                break
        if flag:
            for j in range(MotifLen):
                curscore[i]+=Matrix[j][BASE.index(RNASeq[i+j])]
    return curscore
    #for i in range(len(score)):
    #    score[i]=(curscore[i]-minscore)/(maxscore-minscore)
    #    scorestr=scorestr+"\t"+str(score[i])
    #print RNAID+"\t"+MotifID+scorestr

def SelectTFBS(Lodscore,Maxscore,GeneID,MotifID,MotifLen,Quantile):
    Len=len(Lodscore)
    threshold=Maxscore*Quantile
    if max(Lodscore)>=threshold:
        matchstart=""
        matchend=""
        matchscore=""
        for s in range(Len):
            if Lodscore[s]>=threshold:
                matchstart=matchstart+str(s+1)+","
                matchend=matchend+str(s+MotifLen)+","
                matchscore=matchscore+str(Lodscore[s])+","
        print GeneID+"\t"+str(Len+MotifLen-1)+"\t"+MotifID+"\t"+str(MotifLen)+"\t"+matchstart+"\t"+matchend+"\t"+matchscore+"\t"+str(threshold)
    else:
        print GeneID+"\t"+str(Len+MotifLen-1)+"\t"+MotifID+"\t"+str(MotifLen)+"\tNa\tNa\t"+str(max(Lodscore))+"\t"+str(threshold)

# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main
# ------------------------------------

if __name__=="__main__":
    if len(sys.argv)==1:
        print "Usage: "+sys.argv[0]+" MotifPWM.txt PromoterSequence.fa Quantile> PredictedTFBS.txt"
    else:
        #PermuSize=string.atoi(sys.argv[3])
        Quantile=string.atof(sys.argv[3])
        #timestart=time.clock()
        GeneIDStart=re.compile(r'^>')
        MotifIDStart=re.compile(r'^#')
        MotifPos=re.compile(r'^\d')
        #Initialize Motif 
        MotifID=""
        MotifLen=0
        Matrix=[]
        Cmplodscore=[]
        # Read PFM 
        pfm=open(sys.argv[1])
        for line in pfm:
            line=line.rstrip(os.linesep)
            if MotifIDStart.match(line):
                #After read a new PFM, calculate the score of previous PFM
                if Matrix:
                    #Initialize RNAID and RNASeq
                    GeneID=""
                    PromoterSeq=""
                    #Read fasta
                    fasta=open(sys.argv[2])
                    for row in fasta:
                        row=row.rstrip(os.linesep)
                        if GeneIDStart.match(row):
                            #After read a new sequence, calculate the score of previous Seq
                            if PromoterSeq:
                                Lodscore=GetScore(Matrix,PromoterSeq,MotifLen) if len(PromoterSeq)==MotifLen else [0]
                                Maxscore=GetMaxScore(Matrix,MotifLen)
                                SelectTFBS(Lodscore,Maxscore,GeneID,MotifID,MotifLen,Quantile)
                            row=row.split(">")
                            GeneID=row[1]
                            PromoterSeq=""
                        #Read sequence
                        else:
                            PromoterSeq=PromoterSeq+row
                    #Calculate the score of last motif
                    Lodscore=GetScore(Matrix,PromoterSeq,MotifLen) if len(PromoterSeq)==MotifLen else [0]
                    Maxscore=GetMaxScore(Matrix,MotifLen)
                    SelectTFBS(Lodscore,Maxscore,GeneID,MotifID,MotifLen,Quantile)
                line_tmp=line.split("\t")
                line=line_tmp[0].split("#")
                MotifID=line[1]
                Matrix=[]
                Cmplodscore=[]
            elif MotifPos.match(line):
                line=line.split("\t")
                Matrix.append([])
                Matrix[-1].append(math.log(string.atof(line[1])/HG19BK[0]) if string.atof(line[1])>0 else DEFAULT)
                Matrix[-1].append(math.log(string.atof(line[2])/HG19BK[1]) if string.atof(line[2])>0 else DEFAULT)
                Matrix[-1].append(math.log(string.atof(line[3])/HG19BK[2]) if string.atof(line[3])>0 else DEFAULT)
                Matrix[-1].append(math.log(string.atof(line[4])/HG19BK[3]) if string.atof(line[4])>0 else DEFAULT)
                MotifLen=string.atoi(line[0])
        GeneID=""
        PromoterSeq=""
        #Read fasta
        fasta=open(sys.argv[2])
        for row in fasta:
            row=row.rstrip(os.linesep)
            if GeneIDStart.match(row):
                #After read a new sequence, calculate the score of previous Seq
                if PromoterSeq:
                    Lodscore=GetScore(Matrix,PromoterSeq,MotifLen) if len(PromoterSeq)==MotifLen else [0]
                    Maxscore=GetMaxScore(Matrix,MotifLen)
                    SelectTFBS(Lodscore,Maxscore,GeneID,MotifID,MotifLen,Quantile)
                row=row.split(">")
                GeneID=row[1]
                PromoterSeq=""
                #Read the motif matrix into a dictionary
            else:
                PromoterSeq=PromoterSeq+row
        #Calculate the score of last motif
        Lodscore=GetScore(Matrix,PromoterSeq,MotifLen) if len(PromoterSeq)==MotifLen else [0]
        Maxscore=GetMaxScore(Matrix,MotifLen)
        SelectTFBS(Lodscore,Maxscore,GeneID,MotifID,MotifLen,Quantile)
        #elspsed=(time.clock()-timestart)
        #print "Time used:"+str(elspsed)
