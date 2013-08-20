#!/usr/bin/python
#Last-modified: 20 Feb 2012 05:38:04 PM

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
#import wbed

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
        line=fh.readline()
        #for line in ColumnReader(open(sys.argv[1])):
        while line:
            line=line.split("\t")
            if len(line[4])>17:
                end=string.atoi(line[3])+len(line[4])
                print line[2]+"\t"+line[3]+"\t"+str(end)+"\t"+line[0]+"\t"+str(0)+"\t"+line[1]
            line=fh.readline()
