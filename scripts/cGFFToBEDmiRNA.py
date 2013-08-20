#!/usr/bin/python
#0

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
from wbed import ColumnReader
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
        for line in ColumnReader(open(sys.argv[1])):
            if line[0]!='#':
                line=line.rstrip("\n\r").split("\t")
                group=line[8].split("\"")
                start=string.atoi(line[3])-1
                print "chr"+line[0]+"\t"+str(start)+"\t"+line[4]+"\t"+group[1]+"\t"+str(0)+"\t"+line[6]
