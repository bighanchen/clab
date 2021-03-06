#!/usr/bin/python
#Last-modified: 13 Nov 2012 10:00:56 AM

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
#from wbed import ColumnReader
#import BeautifulSoup

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
            item=line[0].split(":")
            if string.atof(line[2])==100.00 and string.atoi(line[3])==string.atoi(item[1]) and string.atoi(line[6])==1 and string.atoi(line[7])==string.atoi(item[1]):
                if string.atof(line[8])>string.atof(line[9]):
                    start=int(line[9])-1
                    strand="-"
                    print line[1]+"\t"+str(start)+"\t"+line[8]+"\t"+line[0]+"\t"+str(0)+"\t"+strand
                else:
                    start=int(line[8])-1
                    strand="+"
                    print line[1]+"\t"+str(start)+"\t"+line[9]+"\t"+line[0]+"\t"+str(0)+"\t"+strand
            line=fh.readline().rstrip("\n\r")
