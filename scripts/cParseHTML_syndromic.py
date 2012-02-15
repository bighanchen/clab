#!/usr/bin/python
#Last-modified: 15 Feb 2012 02:26:59 PM

""" Module/Script Description

This Module/Script can help to get autism-associated syndromic genes from http://autismkb.cbi.pku.edu.cn 

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
import re
import string
import urllib2
#from wbed import ColumnReader
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
    if len(sys.argv)==2:
        print "Usage: "+sys.argv[0]
    else:
        #for line in ColumnReader(open(sys.argv[1])):
            #line=line.rstrip("\n\r").split("\t")
        for i in range(1,6):
            URL="http://autismkb.cbi.pku.edu.cn/syndromic_gene.php?page="+str(i)
            page=urllib2.urlopen(URL)
            soup=BeautifulSoup.BeautifulSoup(page)
            for j in range(len(soup('td'))):
                if (j+1)%9==1:
                    print soup('td')[j].a.string,
                elif (j+1)%9==2:
                    print soup('td')[j].a.string,
                elif (j+1)%9==3:
                    print soup('td')[j].a.string,
                elif (j+1)%9==4:
                    print soup('td')[j].string,
                elif (j+1)%9==5:
                    print soup('td')[j].string,
                elif (j+1)%9==6:
                    print soup('td')[j].string,
                elif (j+1)%9==7:
                    print soup('td')[j].a.string,
                elif (j+1)%9==0:
                    print soup('td')[j].span.string
