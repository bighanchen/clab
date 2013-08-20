#!/usr/bin/python
#Last-modified: 06 Mar 2012 10:56:09 AM

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
        line=fh.readline().rstrip("\n\r")
        while line:
            line=line.split("\t")
            miRNA=line[0]+line[1]+line[2]+line[3]
            readslist=""
            readscounts=0
            length=string.atoi(line[2])-string.atoi(line[1])
            ft=open(sys.argv[2])
            row=ft.readline().rstrip("\n\r")
            while row:
                row=row.split("\t")
                match=row[0]+row[1]+row[2]+row[3]
                if not cmp(miRNA,match):
                    if string.atoi(row[12])>17:
                        readslist=readslist+row[9]+";"
                        reads=row[9].split(":")
                        readscounts=readscounts+string.atoi(reads[1])
                row=ft.readline().rstrip("\n\r")
            print line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+line[5]+"\t"+readslist+"\t"+str(readscounts)
            line=fh.readline().rstrip("\n\r")
