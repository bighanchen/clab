#!/usr/bin/python
# programmer :Chen Xiaowei 
# usage:

import sys
import string
import re

class snoGPSRead:
    def __init__(self):
        self.id=""
        self.score=""
        self.length=0
        self.seq=""
        self.structure=""
    
def readsnoGPSresult(fp):
    result=snoGPSRead()
    while(True):
        line=fp.readline()
        if not line:
            return result
        if("".join(line[0:4])=="####"):
            lineset=line
            line=fp.readline()
            if("(W)" in line):
                lineset=lineset+line
                for i in range(16):
                    line=fp.readline()
                    lineset=lineset+line
                linegroup=re.compile(r"####(.*)\n>(.*)\n([ATCGU].*)\n#(.*)\n#(.*)\n#(.*)\n#(.*)\n#(.*)\n"\
                r"#(.*)\n#(.*)\n#(.*)\n#upstream(.*)\n([ATCGU]*)\n#downstream(.*)\n([ATCGU]*)\n#stem1(.*)"\
                r"\n([ATCGU]*)\n(.*)\n",re.M).match(lineset).groups()

                q=re.compile("CeN\d*-*\d*\.\d*").match(linegroup[1]).group()
                result.id=q
#                p=re.compile("\s\d{2}\.\d{2}\s").match(linegroup[1]).group()
#                result.score=p
                p=linegroup[1].split("\t")
                result.score=p[1]


                upspace=[" " for i in range(len(linegroup[12]))]
                upspace="".join(upspace)
                downspace=[" " for i in range(len(linegroup[14]))]
                downspace="".join(downspace)
                result.length=len(linegroup[16])
                result.seq=linegroup[12]+linegroup[16]+linegroup[14]
                result.structure=upspace+linegroup[17]+downspace
                return result


if __name__=="__main__":
    if len(sys.argv)==1:
        print "Usage: "+sys.argv[0]+" snoGPS.rst >snoGPS.tab"
    else:
        fp=open(sys.argv[1])
        while(True):
            result=readsnoGPSresult(fp)
            if(result.id == ""):
                break
            print ">"+result.id+"|"+result.score+"|"+str(result.length)
            print result.seq
            print result.structure
