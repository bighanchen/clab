#!/usr/bin/python
#Last-modified: 22 Aug 2011 01:56:58 PM

""" Module/Scripts Description

Copyright (c) 2008 Yunfei Wang <tszn1984@gmail.com>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  experimental
@version: $Revision$
@author:  Yunfei Wang
@contact: tszn1984@gmail.com
"""

# ------------------------------------
# python modules
# ------------------------------------

import sys
import string
from wbed import ColumnReader

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
            #line=line.rstrip("\n\r").split("\t")

