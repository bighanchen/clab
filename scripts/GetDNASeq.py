#!/usr/bin/python
#0

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
from wbed import ColumnReader,Utils

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
        print "Usage: "+sys.argv[0]+" *.2bit *.bed/*.tab > DNA.fa"
    else:
        i=157763
        source=sys.argv[1]
        for genebed in ColumnReader(open(sys.argv[2]),'bed'):
            seq=genebed.getSeq(source)
            rcseq=Utils.rc(seq)
            print ">"+genebed.id
            print rcseq

