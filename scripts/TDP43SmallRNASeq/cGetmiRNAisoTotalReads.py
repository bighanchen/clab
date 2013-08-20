#!/usr/bin/python
#Last-modified: 06 Mar 2012 01:57:56 PM

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
        while line:
            line=line.split("\t")
            miRNA=line[0]+line[1]+line[2]+line[3]
            miRNAreadslist=""
            a5readslist=""
            a3readslist=""
            a53readslist=""
            t5readslist=""
            t3readslist=""
            t53readslist=""
            a3t5readslist=""
            a5t3readslist=""
            miRNAreadscounts=0
            a5readscounts=0
            a3readscounts=0
            a53readscounts=0
            t5readscounts=0
            t3readscounts=0
            t53readscounts=0
            a3t5readscounts=0
            a5t3readscounts=0
            #length=string.atoi(line[2])-string.atoi(line[1])
            ft=open(sys.argv[2])
            row=ft.readline().rstrip("\n\r")
            while row:
                row=row.split("\t")
                match=row[0]+row[1]+row[2]+row[3]
                if not cmp(miRNA,match):
                    if string.atoi(row[12])>16:
                        reads=row[9].split(":")
                        if string.atoi(row[1])==string.atoi(row[7]) and string.atoi(row[2])==string.atoi(row[8]): 
                            miRNAreadslist=miRNAreadslist+row[9]+";"
                            miRNAreadscounts=miRNAreadscounts+string.atoi(reads[2])
                        elif string.atoi(row[1])>string.atoi(row[7]) and string.atoi(row[2])==string.atoi(row[8]):
                            if not cmp(line[5],"+"):
                                a5readslist=a5readslist+row[9]+";"
                                a5readscounts=a5readscounts+string.atoi(reads[2])
                            else:
                                a3readslist=a3readslist+row[9]+";"
                                a3readscounts=a3readscounts+string.atoi(reads[2])
                        elif string.atoi(row[1])==string.atoi(row[7]) and string.atoi(row[2])<string.atoi(row[8]):
                            if not cmp(line[5],"+"):
                                a3readslist=a3readslist+row[9]+";"
                                a3readscounts=a3readscounts+string.atoi(reads[2])
                            else:
                                a5readslist=a5readslist+row[9]+";"
                                a5readscounts=a5readscounts+string.atoi(reads[2])
                        elif string.atoi(row[1])>string.atoi(row[7]) and string.atoi(row[2])<string.atoi(row[8]):
                            a53readslist=a53readslist+row[9]+";"
                            a53readscounts=a53readscounts+string.atoi(reads[2])
                        elif string.atoi(row[1])==string.atoi(row[7]) and string.atoi(row[2])>string.atoi(row[8]):
                            if not cmp(line[5],"+"):
                                t3readslist=t3readslist+row[9]+";"
                                t3readscounts=t3readscounts+string.atoi(reads[2])
                            else:
                                t5readslist=t5readslist+row[9]+";"
                                t5readscounts=t5readscounts+string.atoi(reads[2])
                        elif string.atoi(row[1])<string.atoi(row[7]) and string.atoi(row[2])==string.atoi(row[8]):
                            if not cmp(line[5],"+"):
                                t5readslist=t5readslist+row[9]+";"
                                t5readscounts=t5readscounts+string.atoi(reads[2])
                            else:
                                t3readslist=t3readslist+row[9]+";"
                                t3readscounts=t3readscounts+string.atoi(reads[2])
                        elif string.atoi(row[1])<string.atoi(row[7]) and string.atoi(row[2])>string.atoi(row[8]):
                            t53readslist=t53readslist+row[9]+";"
                            t53readscounts=t53readscounts+string.atoi(reads[2])
                        elif string.atoi(row[1])>string.atoi(row[7]) and string.atoi(row[2])>string.atoi(row[8]):
                            if not cmp(line[5],"+"):
                                a5t3readslist=a5t3readslist+row[9]+";"
                                a5t3readscounts=a5t3readscounts+string.atoi(reads[2])
                            else:
                                a3t5readslist=a3t5readslist+row[9]+";"
                                a3t5readscounts=a3t5readscounts+string.atoi(reads[2])
                        elif string.atoi(row[1])<string.atoi(row[7]) and string.atoi(row[2])<string.atoi(row[8]):
                            if not cmp(line[5],"+"):
                                a3t5readslist=a3t5readslist+row[9]+";"
                                a3t5readscounts=a3t5readscounts+string.atoi(reads[2])
                            else:
                                a5t3readslist=a5t3readslist+row[9]+";"
                                a5t3readscounts=a5t3readscounts+string.atoi(reads[2])
                row=ft.readline().rstrip("\n\r")
            print line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+line[5]+"\t"+miRNAreadslist+a5readslist+a3readslist+a53readslist+t5readslist+t3readslist+t53readslist+a3t5readslist+a5t3readslist+"\t"+str(miRNAreadscounts+a5readscounts+a3readscounts+a53readscounts+t5readscounts+t3readscounts+t53readscounts+a3t5readscounts+a5t3readscounts)+"\t"+miRNAreadslist+"\t"+str(miRNAreadscounts)+"\t"+a5readslist+"\t"+str(a5readscounts)+"\t"+a3readslist+"\t"+str(a3readscounts)+"\t"+a53readslist+"\t"+str(a53readscounts)+"\t"+t5readslist+"\t"+str(t5readscounts)+"\t"+t3readslist+"\t"+str(t3readscounts)+"\t"+t53readslist+"\t"+str(t53readscounts)+"\t"+a3t5readslist+"\t"+str(a3t5readscounts)+"\t"+a5t3readslist+"\t"+str(a5t3readscounts)
            line=fh.readline().rstrip("\n\r")
