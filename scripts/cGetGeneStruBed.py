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
        print "Usage: "+sys.argv[0]+" *.tab exon/intron/gene/UTR3/UTR5> *.bed"
    else:
        if not cmp(sys.argv[2],"gene"):
            for genebed in ColumnReader(open(sys.argv[1]),'tab'):
                print genebed.toBed()
        elif not cmp(sys.argv[2],"exon"):
            for genebed in ColumnReader(open(sys.argv[1]),'tab'):
                for exon in genebed.exons():
                    print exon
        elif not cmp(sys.argv[2],"intron"):
            for genebed in ColumnReader(open(sys.argv[1]),'tab'):
                for intron in genebed.introns():
                    print intron
        elif not cmp(sys.argv[2],"UTR3"):
            for genebed in ColumnReader(open(sys.argv[1]),'tab'):
                print genebed.getUTR3()
        elif not cmp(sys.argv[2],"UTR5"):
            for genebed in ColumnReader(open(sys.argv[1]),'tab'):
                print genebed.getUTR5()
