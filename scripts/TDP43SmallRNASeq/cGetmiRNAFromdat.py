#!/usr/bin/python
#Last-modified: 22 Feb 2012 11:43:36 PM

""" Module/Script Description

This Module/Script can help to extract miRNA information from .dat file

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
        fh=open(sys.argv[1])
        line=fh.readline().rstrip("\n\r")
        hsa=re.compile("^ID   mmu")
        spaces=re.compile("\s+")
        dots=re.compile("\.+")
        miRNA=re.compile("^FT   miRNA")
        stop=re.compile("^//")
        while line:
            if hsa.match(line):
                line=spaces.sub("\t",line)
                line=line.split("\t")
                print line[1]+"\t",
                next=fh.readline().rstrip("\n\r")
                next=fh.readline().rstrip("\n\r")
                next=re.sub(";","",next)
                next=spaces.sub("\t",next)
                next=next.split("\t")
                print next[1]+"\t",
                next=fh.readline().rstrip("\n\r")
                while next:
                    if miRNA.match(next):
                        next=spaces.sub("\t",next)
                        next=dots.sub("\t",next)
                        next=next.split("\t")
                        print next[2]+"\t"+next[3]+"\t",
                        next=fh.readline().rstrip("\n\r")
                        next=next.split("\"")
                        print next[1]+"\t",
                        next=fh.readline().rstrip("\n\r")
                        next=next.split("\"")
                        print next[1]+"\t",
                        next=fh.readline().rstrip("\n\r")
                    if stop.match(next):
                        print
                        break
                    next=fh.readline().rstrip("\n\r")
            line=fh.readline().rstrip("\n\r")
        #for line in ColumnReader(open(sys.argv[1])):
            #line=line.rstrip("\n\r").split("\t")

