#!/usr/bin/python
#0

""" Script Description

This Script can help to search RNA binding sites of RBPs in given RNA sequences.

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
#from wbed import ColumnReader

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------
def GetScore(Matrix,RNAID,RNASeq,MotifID,MotifLen):
    #Initialize score
    scorestr=""
    score=[0]*(len(RNASeq)-MotifLen+1)
    curscore=[0]*(len(RNASeq)-MotifLen+1)
    InforVec=[0]*MotifLen
    maxscore=0
    minscore=0
    for j in range(MotifLen):
        InforVec[j]=Matrix["A"+str(j+1)]*math.log(4*Matrix["A"+str(j+1)])+Matrix["C"+str(j+1)]*math.log(4*Matrix["C"+str(j+1)])+Matrix["G"+str(j+1)]*math.log(4*Matrix["G"+str(j+1)])+Matrix["U"+str(j+1)]*math.log(4*Matrix["U"+str(j+1)])
        maxscore+=InforVec[j]*max(Matrix["A"+str(j+1)],Matrix["C"+str(j+1)],Matrix["G"+str(j+1)],Matrix["U"+str(j+1)])
        minscore+=InforVec[j]*min(Matrix["A"+str(j+1)],Matrix["C"+str(j+1)],Matrix["G"+str(j+1)],Matrix["U"+str(j+1)])
    for i in range(len(RNASeq)-MotifLen+1):
        for j in range(MotifLen):
            curscore[i]+=InforVec[j]*Matrix[RNASeq[(i+j):(i+j+1)]+str(j+1)]
    for i in range(len(score)):
        score[i]=(curscore[i]-minscore)/(maxscore-minscore)
        scorestr=scorestr+"\t"+str(score[i])
    print RNAID+"\t"+MotifID+scorestr
    
# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main
# ------------------------------------

if __name__=="__main__":
    if len(sys.argv)==1:
        print "Usage: "+sys.argv[0]+" RNASequence.fa MotifPWM.txt > MotifSearchScore.txt"
    else:
        #RNA ID start with ">", Motif ID start with "#"
        RNAIDStart=re.compile(r'^>')
        MotifIDStart=re.compile(r'^#')
        #Initialize RNAID and RNASeq
        RNAID=""
        RNASeq=""
        # Read RNA Sequences 
        fasta=open(sys.argv[1])
        for line in fasta:
            line=line.rstrip(os.linesep)
            if RNAIDStart.match(line):
                if RNASeq:
                    #Initialize Motif 
                    MotifID=""
                    MotifLen=0
                    Matrix={}
                    #Read PWM
                    pwm=open(sys.argv[2])
                    for row in pwm:
                        row=row.rstrip(os.linesep)
                        if MotifIDStart.match(row):
                            if Matrix:
                                GetScore(Matrix,RNAID,RNASeq,MotifID,MotifLen)
                            row=row.split("#")
                            MotifID=row[1]
                            Matrix={}
                        elif not cmp(row[0],"A"):
                            row=row.split("\t")
                            MotifLen=len(row)-1
                            for i in range(len(row)-1):
                                key="A"+str(i+1)
                                Matrix[key]=string.atof(row[i+1])
                        elif not cmp(row[0],"C"):
                            row=row.split("\t")
                            for i in range(len(row)-1):
                                key="C"+str(i+1)
                                Matrix[key]=string.atof(row[i+1])
                        elif not cmp(row[0],"G"):
                            row=row.split("\t")
                            for i in range(len(row)-1):
                                key="G"+str(i+1)
                                Matrix[key]=string.atof(row[i+1])
                        elif not cmp(row[0],"U"):
                            row=row.split("\t")
                            for i in range(len(row)-1):
                                key="U"+str(i+1)
                                Matrix[key]=string.atof(row[i+1])
                    GetScore(Matrix,RNAID,RNASeq,MotifID,MotifLen)
                line=line.split(">")
                RNAID=line[1]
                RNASeq=""
            else:
                line=line.rstrip(os.linesep)
                RNASeq=RNASeq+line
        pwm=open(sys.argv[2])
        MotifLen=0
        Matrix={}
        for row in pwm:
            row=row.rstrip(os.linesep)
            if MotifIDStart.match(row):
                if Matrix:
                    GetScore(Matrix,RNAID,RNASeq,MotifID,MotifLen)
                row=row.split("#")
                MotifID=row[1]
            elif not cmp(row[0],"A"):
                row=row.split("\t")
                MotifLen=len(row)-1
                for i in range(len(row)-1):
                    key="A"+str(i+1)
                    Matrix[key]=string.atof(row[i+1])
            elif not cmp(row[0],"C"):
                row=row.split("\t")
                for i in range(len(row)-1):
                    key="C"+str(i+1)
                    Matrix[key]=string.atof(row[i+1])
            elif not cmp(row[0],"G"):
                row=row.split("\t")
                for i in range(len(row)-1):
                    key="G"+str(i+1)
                    Matrix[key]=string.atof(row[i+1])
            elif not cmp(row[0],"U"):
                row=row.split("\t")
                for i in range(len(row)-1):
                    key="U"+str(i+1)
                    Matrix[key]=string.atof(row[i+1])
        GetScore(Matrix,RNAID,RNASeq,MotifID,MotifLen)
