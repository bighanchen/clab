#!/usr/bin/python
#Last-modified: 07 Mar 2012 04:55:18 PM

""" Module/Script Description

This Module/Script can help to

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
#from wbed import ColumnReader

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
        print "Usage: "+sys.argv[0]+" file1 file2..."
    else:
        fh=open(sys.argv[1])
        start={}
        end={}
        chrom={}
        strand={}
        line=fh.readline().rstrip("\n\r")
        while line:
            line=line.split("\t")
            start[line[3]]=line[1]
            end[line[3]]=line[2]
            chrom[line[3]]=line[0]
            strand[line[3]]=line[5]
            line=fh.readline().rstrip("\n\r")
        ft=open(sys.argv[2])
        row=ft.readline().rstrip("\n\r")
        while row:
            row=row.split("\t")
            if not cmp(strand[row[1]],"+"):
                print chrom[row[1]]+"\t"+str(string.atoi(start[row[1]])+string.atoi(row[2])-1)+"\t"+str(string.atoi(start[row[1]])+string.atoi(row[3]))+"\t"+row[4]+"\t"+str(0)+"\t"+strand[row[1]]
            else:
                print chrom[row[1]]+"\t"+str(string.atoi(end[row[1]])-string.atoi(row[3]))+"\t"+str(string.atoi(end[row[1]])-string.atoi(row[2])+1)+"\t"+row[4]+"\t"+str(0)+"\t"+strand[row[1]]
            if row[6]!="":
                if not cmp(strand[row[1]],"+"):
                    print chrom[row[1]]+"\t"+str(string.atoi(start[row[1]])+string.atoi(row[6])-1)+"\t"+str(string.atoi(start[row[1]])+string.atoi(row[7]))+"\t"+row[8]+"\t"+str(0)+"\t"+strand[row[1]]
                else:
                    print chrom[row[1]]+"\t"+str(string.atoi(end[row[1]])-string.atoi(row[7]))+"\t"+str(string.atoi(end[row[1]])-string.atoi(row[6])+1)+"\t"+row[8]+"\t"+str(0)+"\t"+strand[row[1]]
            if len(row)>10:
                if row[10]!="":
                    if not cmp(strand[row[1]],"+"):
                        print chrom[row[1]]+"\t"+str(string.atoi(start[row[1]])+string.atoi(row[10])-1)+"\t"+str(string.atoi(start[row[1]])+string.atoi(row[11]))+"\t"+row[12]+"\t"+str(0)+"\t"+strand[row[1]]
                    else:
                        print chrom[row[1]]+"\t"+str(string.atoi(end[row[1]])-string.atoi(row[11]))+"\t"+str(string.atoi(end[row[1]])-string.atoi(row[10])+1)+"\t"+row[12]+"\t"+str(0)+"\t"+strand[row[1]]
            row=ft.readline().rstrip("\n\r")
