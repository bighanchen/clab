#!/usr/bin/python
#Last-modified: 21 Mar 2012 11:06:22 AM

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
        flag=0
        start=re.compile("^>")
        space=re.compile("\s")
        rRNA=re.compile(".*miRNA.*")

        fh=open(sys.argv[1])
        line=fh.readline().rstrip("\n\r")
        while line:
            if start.match(line) is not None:
                if rRNA.match(line) is not None:
                    item=space.split(line)
                    print item[0]
                    flag=1
                else:
                    flag=0
            else:
                if flag==1:
                    print line
            line=fh.readline().rstrip("\n\r")
        #for line in ColumnReader(open(sys.argv[1])):
            #line=line.rstrip("\n\r").split("\t")
