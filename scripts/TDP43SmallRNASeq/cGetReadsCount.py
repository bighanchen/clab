#!/usr/bin/python
#Last-modified: 02 Mar 2012 05:19:54 PM

""" Module/Script Description

This Module/Script can help to get Reads Count file from Trim-adapter file

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
        j=0
        while line:
            robj=re.compile("^\s+")
            line=robj.sub("",line)
            robj=re.compile("\s")
            line=robj.split(line)
            if string.atoi(line[0])>2:
                print ">Read_"+str(j)+":"+str(len(line[1]))+":"+line[0]
                print line[1]
                j=j+1
            line=fh.readline().rstrip("\n\r")
#        for line in ColumnReader(open(sys.argv[1])):
