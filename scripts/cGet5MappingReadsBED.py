#!/usr/bin/python
#Last-modified: 07 Mar 2012 09:41:08 PM

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
import re
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
        positionNo={}
        fh=open(sys.argv[1])
        line=fh.readline().rstrip("\n\r")
        while line:
            line=line.split("\t")
            if positionNo.has_key(line[0]):
                positionNo[line[0]]=positionNo[line[0]]+1
            else:
                positionNo[line[0]]=1
            line=fh.readline().rstrip("\n\r")
        ft=open(sys.argv[1])
        row=ft.readline().rstrip("\n\r")
        while row:
            row=row.split("\t")
            if positionNo[row[0]]<6:
                if len(row[4])>16 and len(row[4])<27:
                    end=string.atoi(row[3])+len(row[4])
                    print row[2]+"\t"+row[3]+"\t"+str(end)+"\t"+row[0]+"\t"+str(0)+"\t"+row[1]
            row=ft.readline().rstrip("\n\r")
        #       for line in ColumnReader(open(sys.argv[1])):
            #line=line.rstrip("\n\r").split("\t")
