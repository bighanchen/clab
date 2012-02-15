#!/usr/bin/python
#Last-modified: 14 Feb 2012 01:17:01 PM

""" Module/Scripts Description

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
from wbed import ColumnReader
import BeautifulSoup

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
        j=0
        line=fh.readline()
        while line:
            line=line.split("\t")
            if line[8]>line[9]:
                start=int(line[9])-1
                strand="-"
                print line[1]+"\t"+str(start)+"\t"+line[8]+"\t"+str(j)+"\t"+str(0)+"\t"+strand
            else:
                start=int(line[8])-1
                strand="+"
                print line[1]+"\t"+str(start)+"\t"+line[9]+"\t"+str(j)+"\t"+str(0)+"\t"+strand
            j=j+1
            line=fh.readline().rstrip("\n\r")
