#!/usr/bin/python
#0

""" Module/Script Description

This Module/Script can help to calculate partial correlation coefficient combined with information theory.

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
import numpy
#from wbed import ColumnReader

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------
def CalPC(a,b,c):
    return (a-b*c)/numpy.sqrt((1-b**2)*(1-c**2))

def Eps(a,b,c):
    return (CalPC(a,b,c)/a+CalPC(b,a,c)/b+CalPC(c,a,b)/c)/3.0
# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main
# ------------------------------------

if __name__=="__main__":
    if len(sys.argv)==1:
        print "Usage: "+sys.argv[0]+" expression > corrcoef.txt"
    else:
        i=0
        geneid=[]
        fh=open(sys.argv[1])        
        for line in fh:
            line=line.rstrip(os.linesep)
            line=line.split("\t")
            if cmp(line[0],"ProbeID"):
                i+=1
                geneid.append(line[0])
                exp=numpy.array(map(float,line[1:]))
                if i==1:
                    Exp=exp
                else:
                    Exp=numpy.vstack((Exp,exp))
        CorrCoefMat=numpy.corrcoef(Exp)
        #for i in range(10):
        #    for j in range(10):
        #        print CorrCoefMat[i,j]
        nCorr=CorrCoefMat.shape[0]-1
        CorrFlag=numpy.ones(nCorr,int)
        #CorrCoef=numpy.zeros(nCorr,int)
        for i in range(1,nCorr+1):
            for j in numpy.hstack((range(1,i),range(i+1,nCorr+1))):
                epsilon=Eps(CorrCoefMat[0,i],CorrCoefMat[0,j],CorrCoefMat[i,j])
                #print epsilon
                if numpy.abs(CorrCoefMat[0,i])<=numpy.abs(epsilon*CorrCoefMat[0,j]) and numpy.abs(CorrCoefMat[0,i])<=numpy.abs(epsilon*CorrCoefMat[i,j]):
                    CorrFlag[i-1]=0
                    break
            #CorrCoef[i]=CorrCoefMat[0,i]
            #print CorrCoefMat[0,i]
            #print CorrFlag[i-1]
        for i in range(nCorr):
            print "TARDBP\t"+geneid[i]+"\t"+str(CorrFlag[i])+"\t"+str(CorrCoefMat[0,i+1])
