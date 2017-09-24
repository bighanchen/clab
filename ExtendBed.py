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
import os
#from wbed import ColumnReader

# ------------------------------------
# constants
# ------------------------------------
#hg19
CHR_LENGTH={"chr1":249250621,"1":249250621,"chr2":243199373,"2":243199373,"chr3":198022430,"3":198022430,"chr4":191154276,"4":191154276,"chr5":180915260,"5":180915260,"chr6":171115067,"6":171115067,"chr7":159138663,"7":159138663,"chr8":146364022,"8":146364022,"chr9":141213431,"9":141213431,"chr10":135534747,"10":135534747,"chr11":135006516,"11":135006516,"chr12":133851895,"12":133851895,"chr13":115169878,"13":115169878,"chr14":107349540,"14":107349540,"chr15":102531392,"15":102531392,"chr16":90354753,"16":90354753,"chr17":81195210,"17":81195210,"chr18":78077248,"18":78077248,"chr19":59128983,"19":59128983,"chr20":63025520,"20":63025520,"chr21":48129895,"21":48129895,"chr22":51304566,"22":51304566,"chrX":155270560,"X":155270560,"chrY":59373566,"Y":59373566}
#hg38
#CHR_LENGTH={"chr1":248956422,"chr2":242193529,"chr3":198295559,"chr4":190214555,"chr5":181538259,"chr6":170805979,"chr7":159345973,"chr8":145138636,"chr9":138394717,"chr10":133797422,"chr11":135086622,"chr12":133275309,"chr13":114364328,"chr14":107043718,"chr15":101991189,"chr16":90338345,"chr17":83257441,"chr18":80373285,"chr19":58617616,"chr20":64444167,"chr21":46709983,"chr22":50818468,"chrX":156040895,"chrY":57227415}
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
        print "Usage: "+sys.argv[0]+"size .bed > _extended.bed"
    else:
        size=string.atoi(sys.argv[1])
        fh=open(sys.argv[2])
        for line in fh:
            line=line.rstrip(os.linesep)
            line=line.split("\t")
            if line[5]=="+":
                if line[0] in CHR_LENGTH.keys():
                    start=string.atoi(line[1])-size if string.atoi(line[1])-size>0 else 0
                    #end=string.atoi(line[2])+10000 if string.atoi(line[2])+10000<CHR_LENGTH[line[0]] else CHR_LENGTH[line[0]]
                    print line[0]+"\t"+str(start)+"\t"+line[1]+"\t"+line[3]+"\t"+line[4]+"\t"+line[5]
            elif line[5]=="-":
                if line[0] in CHR_LENGTH.keys():
                    #start=string.atoi(line[1])-10000 if string.atoi(line[1])-10000>0 else 0
                    end=string.atoi(line[2])+size if string.atoi(line[2])+size<CHR_LENGTH[line[0]] else CHR_LENGTH[line[0]]
                    print line[0]+"\t"+line[2]+"\t"+str(end)+"\t"+line[3]+"\t"+line[4]+"\t"+line[5]
