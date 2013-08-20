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
@modified from: wbed by tszn1984@gmail.com and zbed by nimezhu@163.com
"""

# ------------------------------------
# python modules
# ------------------------------------

import os,sys,string,random,cPickle,gzip,copy,wRNA,numpy
from math import log,sqrt
from bisect import bisect
from zSeqIO import *
from wWigIO import *

# ------------------------------------
# constants
# ------------------------------------

_numBins   = 37450;
_binLevels = 7;

# bins range in size from 16kb to 512Mb 
# Bin  0          spans 512Mbp,   # Level 1
# Bins 1-8        span 64Mbp,     # Level 2
# Bins 9-72       span 8Mbp,      # Level 3
# Bins 73-584     span 1Mbp       # Level 4
# Bins 585-4680   span 128Kbp     # Level 5
# Bins 4681-37449 span 16Kbp      # Level 6
_binOffsetsExtended = [32678+4096+512+64+8+1, 4096+512+64+8+1, 512+64+8+1, 64+8+1, 8+1, 1, 0]
_binFirstShift = 14;      # How much to shift to get to finest bin. 
_binNextShift  = 3;       # How much to shift to get to next larger bin.

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------


class Bed:
    '''This class is used to process BED format lines. The default BED format is at least three fields.'''
    def __init__(self,x,description=None):
        '''Initiate the bed from either line or list.'''
        try:
            x=x.rstrip("\n\r").split("\t")
        except:
            if isinstance(x[-1],basestring):
                x[-1].rstrip("\n\r")
        self.chr=x[0].strip()
        self.start=int(x[1])
        if self.start<0:
            self.start=0
        self.stop=int(x[2])
        try:
            self.id=x[3]
        except:
            self.id="NONAME"
        try:
            self.score=float(x[4])
        except:
            self.score=1
        try:
            self.strand=x[5]
        except:
            self.strand="."
        try:
            self.otherfields=x[6:]
        except:
            self.otherfields=[]
        self.description=description
    def __str__(self):
        '''Return the bed in basestring format.'''
        string=self.chr+"\t"+str(self.start)+"\t"+str(self.stop)+"\t"+str(self.id)+"\t"+("%-5.2f\t"% self.score)+self.strand
        return string
    def __add__(A,B):
        '''Add Bed A and B together. Please test isOverlap for bed merge!'''
        if not A and not B: return None
        if not A: return B   #B+0=B
        if not B: return A   #A+0=A
        if A.chr==B.chr: #A+B
            start=min(A.start,B.start)
            stop=max(A.stop,B.stop)
            return Bed([A.chr,start,stop,A.id+","+B.id,(A.score*A.length()+B.score*B.length())/(stop-start),A.strand if A.strand==B.strand else "."]) #score= (lenA*scoreA+lenB*scoreB)/len; if strands are different, result strand is ".".
        return A # A and B are not in the same chromosome.
    def __cmp__(self,other):
        return cmp(self.chr,other.chr) or cmp(self.start,other.start) or cmp(self.stop,other.stop) or cmp(other.strand,self.strand) or cmp(self.score,other.score)#'.'<'-'<'+'
    def length(self):
        '''Return the bed range.'''
        return self.stop-self.start
    def isOverlap(self,B):
        '''Return a bool value of the status of if two bed are overlapped.'''
        if not B: return False
        if(self.chr != B.chr) : return False
        if (self.stop <= B.start) : return False
        if (B.stop <= self.start) : return False
        return True
    def overlapLength(self,B):
        '''Return the overlapped length. >0: overlap; =0: neighbour; <0: -distance; -1000000000: not in the same chromosome.'''
        if not B: return -1000000000
        if self.chr==B.chr:
            return self.length()+B.length()-max(self.stop,B.stop)+min(self.start,B.start)
        return -1000000000
    def testBoundary(self):
        '''Swap the start and stop position and reverse the strand if start is bigger than stop.'''
        if self.start>self.stop:
            print >>sys.stderr, "start is bigger than stop, swap them and reverse the strand."
            self.stop,self.start=self.start,self.stop
            if self.strand=="+": 
                self.strand="-"
            elif self.strand=="-":
                self.strand="+"        
    def distance(self,B):
        '''Get the distance between two Beds.'''
        return -self.overlapLenght(B)
        #if(self.chr != bed.chr): return 1000000000 #In different chromosome
        #if(self.isOverlap(bed)): return -1 #Overlapped
        #return min(abs(self.start-bed.stop),abs(self.stop-bed.start))
    def getSeq(self,fn="/disk/Genome/hg18/hg18.2bit"):
        '''Get fasta sequence from 2bit file'''
        if(fn is None): 
            print >>sys.stderr,"2Bit file not specified."
            return ""
        seq=getSeq(fn,self.chr,self.start,self.stop) #Default is lowercase ???
        seq=seq.upper()
        if "-" in self.strand:
            return Utils.rc(seq)
        return seq
    def getWig(self,fn):
        '''get base value from bigWig file.'''
        return getDepthList(fn,self.chr,self.start,self.stop)
    def updownExtend(self,up=0,down=0):
        '''Return the a new Bed with upstream up and downstream down'''
        if self.strand=="-":
            start=self.start-down
            stop=self.stop+up
        else:
            start=self.start-up
            stop=self.stop+down
        tbed=Bed([self.chr,start,stop,self.id+("_up"+str(up) if up!=0 else "")+("_down"+str(down) if down!=0 else ""),self.score,self.strand])
        tbed.testBoundary()
        return tbed
    def stringGraph(self,scale=1):
        '''Illustrate the bed in graph mode.'''
        n=int(self.length()*scale)
        if n==0: n=1
        if self.strand == "+": return ">"*n
        if self.strand == "-": return "<"*n
        return "|"*n
    def setDepth(self,cover):
        '''Add depth for each base.'''
        self.depth=cover
    def strandCmp(self,bed):
        '''Test if the same strand.'''
        if bed.strand == "." or self.strand==".": return "."
        return "+" if self.strand == bed.strand else "-"
    def getBIN(self):
        '''Get the genome BIN.'''
        start=self.start>>_binFirstShift
        stop=(self.stop-1)>>_binFirstShift
        for i in range(_binLevels):
            if start==stop:
                return _binOffsetsExtended[i] + start
            else:
                start >>= _binNextShift
                stop  >>= _binNextShift
        assert "Error! Bed range is out of 512M."
    def toTSS(self):
        '''Change to TSS.'''
        return TSS([self.chr,self.stop-1 if self.strand=='-' else self.start,self.id,self.score,self.strand])
    def toSite(self,stype='TSS'):
        '''Change to Site.'''
        if stype=='TSS':
            return Site([self.chr,self.stop-1 if self.strand=='-' else self.start,self.id,self.score,self.strand],stype)
        if stype=='EnzymeDigest':
            return Site([self.chr,self.stop if self.strand=='-' else self.start,self.id,self.score],stype)
        if stype=='TIS':
            pass
        if stype=='TTS':
            pass
        return None
        
class BedCoverage(Bed):
    '''Bed6 format with coverage infomation.'''
    def __init__(self,x,description=None):
        '''Initiation'''
        Bed.__init__(self,x,description)
        if isinstance(self.otherfields[0],str):
            self.depth=[float(i) for i in self.otherfields[0].split(',')]
            self.otherfields=self.otherfields[1:]
        elif isinstance(self.otherfields[0],list):
            self.depth=self.otherfields[0]
            self.otherfields=self.otherfields[1:]
        else:
            self.depth=[self.score for i in xrange(self.length())]

class BedList(list):
    '''List class for hold line objects.'''
    def __init__(self,data=[],description=None):
        '''Initiate data by list class.'''
        list.__init__(self,data)
        self.sorted=0
        self.description=description
    def readfile(self,infile,format='bed'):
        '''Read data fromfile by  ColumnReader.'''
        for item in ColumnReader(infile,format):
            self.append(item)
    def sort(self):
        '''sort BedList.'''
        list.sort(self)
        self.sorted=1
    def bisect(self,item):
        '''Find the nearest item for comparation.'''
        if not self.sorted:
            self.sort()
        return bisect(self,item)
    def clear(self):
        '''Clear List.'''
        del self[:]
        self.sorted=0
    def mergeSort(bedfiles): #generator
        '''Merge Beds from multiple files. The Bed in each file should be sorted.'''
        print >>sys.stderr, "Make sure each file input is sorted!"
        if isinstance(bedfiles,str):
            bedfiles=[bedfiles]
        if len(bedfiles)==0:
            print >>sys.stderr, "No bed file names provided...."
        elif len(bedfiles)==1: # for single file
            for tbed in ColumnReader(bedfiles[0],'bed'):
                yield tbed
        else: #for multiple file
            cr=[] #List for iteraters
            bedfiledict={} #record bedfile index for iteration
            beds=BedList()
            for index,bedfile in enumerate(bedfiles):
                bedfiledict[bedfile]=index
                cr.append(ColumnReader(bedfile,'bed'))
                try:
                    beds.append(cr[index].next())
                    beds[index].description=bedfile #record the source file number for iteration
                except Exception,e:
                    print >>sys.stderr,"Read Bed error!",e
                    raise
            beds.sort()
            while True:
                if len(beds)>0:
                    yield beds[0] # yield the minimum bed
                    try:
                        tbed=cr[bedfiledict[beds[0].description]].next()
                        tbed.description=beds[0].description
                        beds.insert(beds.bisect(tbed),tbed) # insert(pos,item) insert tbed in  the right position
                    except StopIteration:
                        print >>sys.stderr,bedfiles[bedfiledict[beds[0].description]]+" is finished..."
                    del beds[0]                
                else:
                    break
    mergeSort=staticmethod(mergeSort)
    def mergeBeds(bedfiles,forcestrand=False): #generator
        '''Merge Beds from multiple files. The overlapped beds are combined. The Bed in each file should be sorted.'''
        beds=[None,None]
        bedcount=[0,0]
        bedindex={".":0,"+":0,"-":(1 if forcestrand else 0)}
        for tbed in BedList.mergeSort(bedfiles):
            if tbed.isOverlap(beds[bedindex[tbed.strand]]):
                beds[bedindex[tbed.strand]]+=tbed
            else:
                if beds[bedindex[tbed.strand]]:
                    bedcount[bedindex[tbed.strand]]+=1
                    beds[bedindex[tbed.strand]].description=beds[bedindex[tbed.strand]].id
                    beds[bedindex[tbed.strand]].id="Region_"+str(bedcount[bedindex[tbed.strand]])
                    yield beds[bedindex[tbed.strand]]
                beds[bedindex[tbed.strand]]=tbed
        for tbed in beds:
            if tbed:
                bedcount[bedindex[tbed.strand]]+=1
                tbed.description="Region_"+str(bedcount[bedindex[tbed.strand]])
                tbed.id,tbed.description=tbed.description,tbed.id
                yield tbed
        assert "Reach this line."
    mergeBeds=staticmethod(mergeBeds)
            
class GeneBed(Bed):
    '''UCSC GenePred format.'''
    def __init__(self,x,description=None):
        '''Initiate from GeneBed lines. GeneBed names column are not allowed to be numbers.'''
        try:
            self.bin=int(x[0])
            if self.bin<10000:
                x=x[1:]
        except:
            pass
        self.id=x[0]
        self.chr=x[1]
        self.strand=x[2]
        self.start=int(x[3])
        self.stop=int(x[4])
        self.txstart=int(x[5])
        self.txstop=int(x[6])
        self.exoncount=int(x[7])
        if isinstance(x[8],basestring):
            self.exonstarts=[int(p) for p in x[8].split(",")[0:-1]]
            self.exonstops=[int(p) for p in x[9].split(",")[0:-1]]
        else:
            self.exonstarts=[int(p) for p in x[8]]
            self.exonstops=[int(p) for p in x[9]]                    
        self.score=0
        try:
            self.protein_id=x[10]
        except:
            self.protein_id=None
        self.description=description
    def __str__(self):
        '''Return GeneBed line.'''
        return "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s,\t%s," % (self.id,self.chr,self.strand,self.start,self.stop,self.txstart,self.txstop,self.exoncount,",".join([str(p) for p in self.exonstarts]),",".join([str(p) for p in self.exonstops]))
    def toBed(self):
        '''Transform to Bed format.'''
        return Bed([self.chr,self.start,self.stop,self.id,self.score,self.strand])
    def getExon(self,i):
        '''Get the Bed format of ith exon.'''
        if i>self.exoncount or i<1:
            return None
        p= (self.exoncount-i if self.strand=='-' else i-1)
        return Bed([self.chr,self.exonstarts[p],self.exonstops[p],self.id+":exon_"+str(i),0,self.strand])
    def getIntron(self,i):
        '''Get the Bed format of ith intron.'''
        if i>0 and i<self.exoncount:
            p=(self.exoncount-i if self.strand=='-' else i)
            #Notice: GeneBed is 1 based and Bed is 0 based.
            return Bed([self.chr,self.exonstops[p-1],self.exonstarts[p],self.id+":intron_"+str(i),0,self.strand])
        else:
            return None
    def _getUTRs(self,end='left'):
        '''Get the UTR.'''
        if end=='left': #Get the UTR located in the left end of chromosome
            if self.start==self.txstart: #No UTR
                return None
            utr=copy.deepcopy(self)
            for i in range(self.exoncount):
                if self.txstart<=self.exonstops[i]: #the txStart site locates in (i+1)th exon.
                    break
            utr.exonstarts=utr.exonstarts[0:i+1]
            utr.exonstops=utr.exonstops[0:i+1]
            utr.stop=utr.txstop=utr.exonstops[i]=self.txstart
            utr.txstart=self.start
            utr.exoncount=i+1
            return utr
        else: #get the UTR located in the right end of chromosome
            if self.stop==self.txstop: #No UTR
                return None
            utr=copy.deepcopy(self)
            for i in range(self.exoncount-1,-1,-1):
                if self.txstop>=self.exonstarts[i]:
                    break
            utr.exonstarts=utr.exonstarts[i:self.exoncount]
            utr.exonstops=utr.exonstops[i:self.exoncount]
            utr.start=utr.txstart=utr.exonstarts[0]=self.txstop
            utr.txstop=self.stop
            utr.exoncount=self.exoncount-i
            return utr
    def getUTR5(self):
        '''Get the 5UTR.'''
        if self.strand=='-':
            utr5=self._getUTRs('right')
        else:
            utr5=self._getUTRs('left')
        if utr5:
            utr5.id+=':UTR5'
        return utr5
    def getUTR3(self):
        '''Get the 3'UTR.'''
        if self.strand=='-':
            utr3=self._getUTRs('left')
        else:
            utr3=self._getUTRs('right')
        if utr3:
            utr3.id+=':UTR3'
        return utr3
    def overlapLength(self,B):
        '''Return the overlap length between genebed and bed. Only exons are considered.'''
        l=0
        for exon in self.exons():
            l+=exon.overlapLength(B)
        return l
    def exons(self):
        '''Iterate all exons.'''
        for i in range(1,self.exoncount+1):
            yield self.getExon(i)
    def introns(self):
        '''Iterate all introns.'''
        for i in range(1,self.exoncount):
            yield self.getIntron(i)
    def getcDNALength(self):
        '''Return cdna length.'''
        l=0
        cdsbed=Bed([self.chr,self.txstart,self.txstop])
        for exon in self.exons():
            l+=exon.length()
        return l
    def getcDNASeq(self,fn="/disk/Genome/hg18/hg18.2bit"):
        '''Get cDNA Sequence.'''
        seq=""
        for i in range(self.exoncount):
            seq+=self.getExon(i+1).getSeq(fn)
        return seq
    def getCDSLength(self):
        '''Return CDS length.'''
        l=0
        for exon in self.exons():
            l+=cdsbed.OverlapLength(exon)
        return l
    def getCDSSeq(self,fn="/disk/Genome/hg18/hg18.2bit"):
        '''get CDS sequence.'''
        seq=""
        cdsbed=Bed([self.chr,self.txstart,self.txstop])
        for exon in self.exons():
            if cdsbed.isOverlap(exon):
                seq+=Bed([self.chr,max(self.txstart,exon.start),min(self.txstop,exon.stop),'',0,self.strand]).getSeq(fn)
        return seq
    def getGeneLength(self):
        '''Return gene length.'''
        return self.stop-self.start+1

class GeneBedList(BedList): #useless
    '''A list for Genes.'''
    def __init__(self,data=[],description=None):
        BedList.__init__(self,data,description)
        
class Site:
    '''site.'''
    def __init__(self,x,stype=None,description=None):
        '''Initiation'''
        self.chr=x[0]
        self.pos=int(x[1])
        try:
            self.id=x[2]
        except:
            self.id="NONAME"
        try:
            self.score=float(x[3])
        except:
            self.score=1.0
        try:
            self.strand=x[4]
        except:
            self.strand="."
        try:
            self.otherfields=x[5:]
        except:
            self.otherfileds=[]
        self.type=stype
        self.description=description
    def __str__(self):
        '''Return string of Site.'''
        return "%s\t%d\t%s\t%-5.3f\t%s\t%s" %(self.chr,self.pos,self.id,self.score,self.strand,self.type if self.type else '')
    def __cmp__(self,site):
        '''Compare self and other Site.'''
        return cmp(self.chr,site.chr) or cmp(self.pos,site.pos) or cmp(site.strand,self.strand) or cmp(self.score,site.score)

class SiteList(list):
    '''A list of Sites.'''
    def __init__(self,data=[],description=None):
        '''Initiate data by list class.'''
        list.__init__(self,data)
        self.sorted=0
        self.description=description
    def readfile(self,infile):
        '''Read data fromfile by  ColumnReader.'''
        for line in ColumnReader(infile):
            self.append(Site(item))
    def mergeSites(self):
        '''Merge sites in the same position.'''
        if not self.sorted:
            self.sort()
        tsite=None
        for site in self:
            if not tsite:
                tsite=site
            elif site==tsite:
                tsite.score+=site.score
            else:
                yield tsite
                tsite=site
        if tsite:
            yield tsite

class Fasta:
    '''Fasta format.'''
    def __init__(self,name=None,seq=None,description=None):
        '''Initiate the fasta record.'''
        self.id=name
        self.seq=seq
        self.description=description
        self.evalue=None
        self.structure=None
    def length(self):
        '''get the length of sequence.'''
        return len(self.seq)
    def __str__(self):
        '''String for output of Fasta.'''
        return ">"+self.id+"\n"+self.seq
    def Fold(self,constraints="",temperature=20):
        '''Fold the struture for the sequence with stem-loop constraints at fixed temperature.'''
        (self.evalue,self.structure)=wRNA.fold(self.seq,constraints,temperature)

class ColumnReader:
    '''Read column files. Support format:Bed,Gene/GeneBed/Tab/GenePred,Bowtie and SOAP'''
    def __init__(self,infile,format=None):
        #Initiate file handle from pwd.
        if isinstance(infile,basestring):
            try:
                self.infile=open(infile)
            except Exception,e:
                print >>sys.stderr,"Can't open file.",e
                raise
        #Initiate file handle from other file handle
        else:
            self.infile=infile
        try:
            self.line=self.infile.readline()
        except Exception,e:
            print >>sys.stderr,"Can't read file:",e
            raise
        #Move to the first valid line
        while self.line and self.line[0]=='#':
            self.line=self.infile.readline()
        self.format=format.lower() if format else None
    def next(self):
        '''Read one valid line from file.'''
        if self.line:
            x=self.line.rstrip(os.linesep).split('\t')
            self.line=self.infile.readline()
            if not self.format: return x
            if self.format=='bed':
                return Bed(x)
            if self.format=='bowtie':
                return Utils.BowtieToBed(x)
            if self.format=='soap':
                return Utils.SOAPToBed(x)
            if self.format in ['gene','genebed','genepred','tab']:
                return GeneBed(x)
        else:
            self.infile.close()
            raise StopIteration
    def __iter__(self)       '''Iterator'''
        return self

class SeqReader:
    '''Read DNA/RNA sequence in Fasta/Fastaq/CSF format.'''
    def __init__(self,infile,format='fasta'):
        '''Initiate file handle.'''
        try:
            self.infile=open(infile)
        except:
            self.infile=infile
        #self.line=self.infile.readline()
        #while self.line and "#" in self.line:
        #    self.line=self.infile.readline()
        self.format=format.lower()
        #self.record=None
    def next(self):
        '''Next Record.'''
        line=self.infile.readline()
        if line:
            if self.format=='fasta' or self.format=='csf':
                description=line.lstrip('>').rstrip("\n\r")
                name=description.split()[0] #split by blanks
                os.linesep='>'
                line=self.infile.readline()
                os.linesep='\n'
                seq=line.rstrip('>').replace('\n','').replace('\r','').upper()
                return Fasta(name,seq,description)
            if self.format=='fastaq':
                name=line.lstrip('@').rstrip("\n\r")
                seq=self.infile.readline().rstrip("\n\r")
                self.infile.readline()
                description=self.infile.readline().rstrip("\n\r")
                return Fasta(name,seq,description)
                return 
            assert False, "File format error!"
        else:
            raise StopIteration
    def __iter__(self):
        return self
                          
class BedMap(dict):
    '''Map Beds with chromsomes and BIN values.'''
    def __init__(self,organism='ce6'):
        '''Default genome is ce6.'''
        dict.__init__(self)
        for chro in Utils.genomeSize(organism).keys():
            self[chro]={}
    def findOverlap(self,bed,fraction=0.5,forcestrand=False):
        '''Find overlaps with bed and put into bedlist.'''
        bedlist=BedList()
        minover=bed.length()*fraction
        maxover=0
        startBin,stopBin= bed.start >> _binFirstShift, (bed.stop-1) >> _binFirstShift
        for i in range(_binLevels):
            offset = _binOffsetsExtended[i]
            for j in xrange(startBin+offset,stopBin+offset+1):
                if not self[bed.chr].has_key(j):
                    continue
                for item in self[bed.chr][j]:
                    overlen=bed.overlapLength(item)
                    if overlen>=minover:
                        if not forcestrand or (forcestrand and bed.strand==item.strand):
                            if maxover<overlen:
                                maxover=overlen
                                bed.description=item
                            bedlist.append(item)
            startBin >>= _binNextShift
            stopBin >>= _binNextShift
        return bedlist
    def intersectBed(self,bed,fraction=0.5,outputoption=1,forcestrand=False):
        '''Intersect bed with items in BedMap.'''
        #outputoption=0 for No overlap,1 for best overlap and 2 for all overlaps}
        bedlist=self.findOverlap(bed,fraction,forcestrand)
        if outputoption==0: #No overlap
            return (True if len(bedlist)==0 else False)
        if outputoption==1: #best overlap
            return bed.description # if description==None, no overlap found.
        #return all overlaps
        return bedlist
    def loadBedToMap(self,bedlist=None,bedtype='bed'):
        '''Load Bed to BedMap from either Bedlist or ColumnReader handle or bedfile. Load one bedlist once.'''
        if bedlist:
            if isinstance(bedlist,str) or isinstance(bedlist,file):
                bedlist=ColumnReader(bedlist,bedtype)
            for bed in bedlist:
                _bin=bed.getBIN()
                self[bed.chr].setdefault(_bin,BedList())
                self[bed.chr][_bin].append(bed)
    def __iter__(self): #useless !!!
        '''Output the bedMap by iteration'''
        for chro in sorted(self.keys()):
            for _bin in sorted(self[chro].keys()):
                for bed in self[chro][_bin]:
                    yield bed

class Utils:
    def rc(seq):
        comps = {'A':"T", 'C':"G", 'G':"C", 'T':"A",
                'B':"V", 'D':"H", 'H':"D", 'K':"M",
                'M':"K", 'R':"Y", 'V':"B", 'Y':"R",
                'W':'W', 'N':'N', 'S':'S'}
        return ''.join([comps[x] for x in seq.upper()[::-1]])
    rc=staticmethod(rc)
    def MW(seq):
        mws={'A':313.21,'C':289.19,'G':329.21,'T':304.2,'I':314.2,'N':308.95,'R':321.21,'Y':296.69,'M':301.2,'K':316.7,'S':309.2,'W':308.71,'H':302.2,'B':307.53,'D':315.54,'V':310.53,'p':79.98,'X':0,'U':290.17}
        return sum([mws[x] for x in seq.upper()])
    MW=staticmethod(MW)
    def TM(seq,Na=100):
        seq=seq.upper()
        if len(seq)<25:
            tm={'A':2,'T':2,'C':4,'G':4}
            return sum([tm[x] for x in seq])
        else:
            N=float(len(seq))
            gc=(seq.count('C')+seq.count('G'))/N
            return 81.5+16.6*log(Na/1000.0,10)+0.41*gc+600.0/N
    TM=staticmethod(TM)
    def toRNA(seq):
        return seq.upper().replace('T','U')
    toRNA=staticmethod(toRNA)
    def toDNA(seq):
        return seq.upper().replace('U','T')
    toDNA=staticmethod(toDNA)
    def toProtein(seq,table):
        '''Translate DNA or RNA to protein according to standard translation table.'''
        seq=seq.upper().rstrip()
        if "U" in seq:
            seq=toDNA(seq)
        if len(seq)%3!=0:
            print >>sys.stderr, "Sequcence length should be 3*N."
            return None
        p=""
        for i in xrange(len(seq)/3):
            p+=table[seq[i*3:(i+1)*3]]
        return p
    toProtein=staticmethod(toProtein)
    def fastaToCSFasta(seq,starter='T'):
        trans= [ ['0','1','2','3'], ['1','0','3','2'], ['2','3','0','1'], ['3','2','1','0']]
        bases= {'A':0,'C':1,'G':2,'T':3}
        csf=starter
        seq=seq.upper()
        tseq=starter+seq
        for i in range(len(seq)):
            csf+=trans[bases[tseq[i]]][bases[tseq[i+1]]]
        return csf
    fastaToCSFasta=staticmethod(fastaToCSFasta)
    def csFastaToFasta(seq):
        trans=[ [0,1,2,3], [1,0,3,2], [2,3,0,1], [3,2,1,0]]
        bases={'A':0,'C':1,'G':2,'T':3}
        csbases={0:'A',1:'C',2:'G',3:'T'}
        f=seq[0]
        for i in range(len(seq)-1):
            f+=csbases[trans[bases[f[i]]][int(seq[i+1])]]
        return f[1:]
    csFastaToFasta=staticmethod(csFastaToFasta)
    def BowtieToBed(x):
        '''Bowtie map result to Bed.Bowtie is 0 based.'''
        x[3]=int(x[3])
        return Bed([x[2],x[3],x[3]+len(x[4]),x[0],1,x[1]])
    BowtieToBed=staticmethod(BowtieToBed)
    def SOAPToBed(x):
        '''SOAP map result to Bed. SOAP is 1 based.'''
        x[8]=int(x[8])-1
        return Bed([x[7],x[8],x[8]+int(x[5]),x[0],1,x[6]])
    SOAPToBed=staticmethod(SOAPToBed)
    def genomeSize(gversion='ce6'):
        '''Genome size dictionary.'''
        genome={}
        genome['ce6']={'chrI':15072421,'chrII':15279323,'chrIII':13783681,'chrIV':17493785,'chrM':13794,'chrV':20919568,'chrX':17718854}
        return genome[gversion]
    genomeSize=staticmethod(genomeSize)
    def translateTables(tabletype="standard"):
        '''Translation tables. @ for stop condons.'''
        tables={}
        tables["standard"]={
            'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
            'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
            'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
            'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
            'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I',
            'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T',
            'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K',
            'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A',
            'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D',
            'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
            'GGG': 'G', 'TAA': '@', 'TAG': '@', 'TGA': '@'}
        return tables[tabletype]
    translateTables=staticmethod(translateTables)
    def dump(handle,filename=None):
        '''Dump object into hard disk.'''
        if not filename:
            filename="tmp"+str(random.randint(0,99))+".dmp.gz"
        fh=gzip.open(filename,'wb')
        cPickle.dump(handle,fh)
        fh.close()
    dump=staticmethod(dump)
    def load(filename):
        '''Load dumpped object into memory.'''
        fh=gzip.open(filename,'rb')
        th=cPickle.load(fh)
        fh.close()
        return th
    load=staticmethod(load)
    def GFFReader_celegans(filename):
        '''Read GFF format file and transform to GeneBed format.'''
        chroTab={'I':'chrI','II':'chrII','III':'chrIII','IV':'chrIV','MtDNA':'chrM','V':'chrV','X':'chrX'}
        gchr=None
        for line in ColumnReader(open(filename)):
            order=line[8][0:2]
            name=line[8].split(';')[0].split(':')[1]
            if order=='ID':
                if gchr:
                    gexons.sort()
                    if gtxstart-1 in gexons:
                        cur=gexons.index(gtxstart)
                        del gexons[cur]
                        del gexons[cur-1]
                    if gtxstop+1 in gexons:
                        cur=gexons.index(gtxstop+1)
                        del gexons[cur]
                        del gexons[cur-1]
                    exonstarts=[]
                    exonstops=[]
                    exoncount=len(gexons)/2
                    for i in range(exoncount):
                        exonstarts.append(gexons[2*i]-1)
                        exonstops.append(gexons[i*2+1])
                    yield GeneBed([gname,gchr,gstrand,gstart-1,gstop,gtxstart-1,gtxstop,exoncount,exonstarts,exonstops])
                #Initiate new record
                gchr=chroTab[line[0]]
                gname=name
                gstrand=line[6]
                gstart=int(line[3])
                gstop=int(line[4])
                gtxstart=gstop
                gtxstop=gstart
                gexons=[]
            elif order=='Pa':
                start=int(line[3])
                stop=int(line[4])
                gexons.append(start)
                gexons.append(stop)
                if line[2]=='coding_exon':
                    gtxstart=min(gtxstart,start)
                    gtxstop=max(gtxstop,stop)
        gexons.sort()
        if gtxstart-1 in gexons:
            cur=gexons.index(gtxstart)
            del gexons[cur]
            del gexons[cur-1]
        if gtxstop+1 in gexons:
            cur=gexons.index(gtxstop+1)
            del gexons[cur]
            del gexons[cur-1]
        exonstarts=[]
        exonstops=[]
        exoncount=len(gexons)/2
        for i in range(exoncount):
            exonstarts.append(gexons[2*i]-1)
            exonstops.append(gexons[i*2+1])
        yield GeneBed([gname,gchr,gstrand,gstart-1,gstop,gtxstart-1,gtxstop,exoncount,exonstarts,exonstops])
    GFFReader=staticmethod(GFFReader)
    def GFF3Reader(filename):
        '''Read GFF3 format file and transform to GenePred format.'''
        #chroTab={'I':'chrI','II':'chrII','III':'chrIII','IV':'chrIV','MtDNA':'chrM','V':'chrV','X':'chrX'}
        tschr=None
        for line in ColumnReader(open(filename)):
            order=line[8][0:6]
            if not cmp(order,"ID=rna"):
                if tschr:
                    exonstarts=[]
                    exonends=[]
                    exoncount=len(tsexons)/2
                    for i in range(exoncount):
                        if not cmp(tsstrand,"+"):
                            exonstarts.append(tsexons[2*i]-1)
                            exonends.append(tsexons[i*2+1])
                        else:
                            exonstarts.append(tsexons[len(tsexons)-2*i-2]-1)
                            exonends.append(tsexons[len(tsexons)-i*2-1])
                    if not cmp(tsstrand, "+"):
                        tsstart=tsexons[0]-1
                        tsend=tsexons[-1]
                        if tscds:
                            cdsstart=tscds[0]-1
                            cdsend=tscds[-1]
                        else:
                            cdsstart=tsend
                            cdsend=tsend
                    else:
                        tsstart=tsexons[-2]-1
                        tsend=tsexons[1]
                        if tscds:
                            cdsstart=tscds[-2]-1
                            cdsend=tscds[1]
                        else:
                            cdsstart=tsend
                            cdsend=tsend
                    description=tsdesc+"|"+genedesc
                    yield GeneBed([tsname,tschr,tsstrand,tsstart,tsend,cdsstart,cdsend,exoncount,exonstarts,exonends,0,genename,None,None,None,description])
                #Initiate new record
                tschr=line[0]
                tsstrand=line[6]
                tsexons=[]
                tscds=[]
                tsname=line[8].split(';')[1].split('=')[1]
                tsdesc=line[8]
                genename=lastgenename
                genedesc=lastgenedesc
            elif not cmp(line[8][0:7],"ID=gene"):
                lastgenename=line[8].split(';')[1].split('=')[1]
                lastgenedesc=line[8]
            else:
                if not cmp(line[2],"exon"):
                    tsexons.append(int(line[3]))
                    tsexons.append(int(line[4]))
                elif not cmp(line[2],"CDS"):
                    tscds.append(int(line[3]))
                    tscds.append(int(line[4]))
        exonstarts=[]
        exonends=[]
        exoncount=len(tsexons)/2
        for i in range(exoncount):
            if not cmp(tsstrand,"+"):
                exonstarts.append(tsexons[2*i]-1)
                exonends.append(tsexons[i*2+1])
            else:
                exonstarts.append(tsexons[len(tsexons)-2*i-2]-1)
                exonends.append(tsexons[len(tsexons)-i*2-1])
        if not cmp(tsstrand,"+"):
            tsstart=tsexons[0]-1
            tsend=tsexons[-1]
            if tscds:
                cdsstart=tscds[0]-1
                cdsend=tscds[-1]
            else:
                cdsstart=tsend
                cdsend=tsend
        else:
            tsstart=tsexons[-2]-1
            tsend=tsexons[1]
            if tscds:
                cdsstart=tscds[-2]-1
                cdsend=tscds[1]
            else:
                cdsstart=tsend
                cdsend=tsend
        description=tsdesc+"|"+genedesc
        yield GeneBed([tsname,tschr,tsstrand,tsstart,tsend,cdsstart,cdsend,exoncount,exonstarts,exonends,0,genename,None,None,None,description])
    GFF3Reader=staticmethod(GFF3Reader)
    def GTFReader(filename):
        '''Read GTF format file and transform to GeneBed format.'''
        #chroTab={'I':'chrI','II':'chrII','III':'chrIII','IV':'chrIV','MtDNA':'chrM','V':'chrV','X':'chrX'}
        tscrpt=[]
        tschr=None
        for line in ColumnReader(open(filename)):
            tsid=line[8].split(';')[1].split('"')[1]
            if tsid not in tscrpt:
                tscrpt.append(tsid)
                if tschr:
                    tsid=tscrpt[-2]
                    exonstarts=[]
                    exonends=[]
                    exoncount=len(tsexons)/2
                    for i in range(exoncount):
                        if not cmp(tsstrand,"+"):
                            exonstarts.append(tsexons[2*i]-1)
                            exonends.append(tsexons[i*2+1])
                        else:
                            exonstarts.append(tsexons[len(tsexons)-2*i-2]-1)
                            exonends.append(tsexons[len(tsexons)-i*2-1])
                    if not cmp(tsstrand, "+"):
                        tsstart=tsexons[0]-1
                        tsend=tsexons[-1]
                        if tscds:
                            cdsstart=tscds[0]-1
                            cdsend=tscds[-1]
                        else:
                            cdsstart=tsend
                            cdsend=tsend
                        if startcodon:
                            cdsstart=startcodon[0]-1
                        if stopcodon:
                            cdsend=stopcodon[1]
                    else:
                        tsstart=tsexons[-2]-1
                        tsend=tsexons[1]
                        if tscds:
                            cdsstart=tscds[-2]-1
                            cdsend=tscds[1]
                        else:
                            cdsstart=tsend
                            cdsend=tsend
                        if startcodon:
                            cdsend=startcodon[1]
                        if stopcodon:
                            cdsstart=stopcodon[0]-1
                    description="Gene_id:"+geneid+";"+"Gene_name:"+genename+";"+"Gene_biotype:"+genebiotype+";"+"Transcript_id:"+tsid+";"+"Transcript_name:"+tsname+";"+"Transcript_biotype:"+tsbiotype+";"
                    yield GeneBed([tsid,tschr,tsstrand,tsstart,tsend,cdsstart,cdsend,exoncount,exonstarts,exonends,0,geneid,None,None,None,description])
                #Initiate new record
                tschr=line[0]
                tsstrand=line[6]
                tsexons=[]
                tscds=[]
                startcodon=[]
                stopcodon=[]
                tsname=line[8].split(';')[5].split('"')[1]
                tsbiotype=line[1]
                geneid=line[8].split(';')[0].split('"')[1]
                genename=line[8].split(';')[3].split('"')[1]
                genebiotype=line[8].split(';')[4].split('"')[1]
                if not cmp(line[2],"exon"):
                    tsexons.append(int(line[3]))
                    tsexons.append(int(line[4]))
                elif not cmp(line[2],"CDS"):
                    tscds.append(int(line[3]))
                    tscds.append(int(line[4]))
                elif not cmp(line[2],"start_codon"):
                    startcodon.append(int(line[3]))
                    startcodon.append(int(line[4]))
                elif not cmp(line[2],"stop_codon"):
                    stopcodon.append(int(line[3]))
                    stopcodon.append(int(line[4]))
            else:
                if not cmp(line[2],"exon"):
                    tsexons.append(int(line[3]))
                    tsexons.append(int(line[4]))
                elif not cmp(line[2],"CDS"):
                    tscds.append(int(line[3]))
                    tscds.append(int(line[4]))
                elif not cmp(line[2],"start_codon"):
                    startcodon.append(int(line[3]))
                    startcodon.append(int(line[4]))
                elif not cmp(line[2],"stop_codon"):
                    stopcodon.append(int(line[3]))
                    stopcodon.append(int(line[4]))
        exonstarts=[]
        exonends=[]
        exoncount=len(tsexons)/2
        for i in range(exoncount):
            if not cmp(tsstrand,"+"):
                exonstarts.append(tsexons[2*i]-1)
                exonends.append(tsexons[i*2+1])
            else:
                exonstarts.append(tsexons[len(tsexons)-2*i-2]-1)
                exonends.append(tsexons[len(tsexons)-i*2-1])
        if not cmp(tsstrand,"+"):
            tsstart=tsexons[0]-1
            tsend=tsexons[-1]
            if tscds:
                cdsstart=tscds[0]-1
                cdsend=tscds[-1]
            else:
                cdsstart=tsend
                cdsend=tsend
        else:
            tsstart=tsexons[-2]-1
            tsend=tsexons[1]
            if tscds:
                cdsstart=tscds[-2]-1
                cdsend=tscds[1]
            else:
                cdsstart=tsend
                cdsend=tsend
        description="Gene_id:"+geneid+";"+"Gene_name:"+genename+";"+"Gene_biotype:"+genebiotype+";"+"Transcript_id:"+tsid+";"+"Transcript_name:"+tsname+";"+"Transcript_biotype:"+tsbiotype+";"
        yield GeneBed([tsid,tschr,tsstrand,tsstart,tsend,cdsstart,cdsend,exoncount,exonstarts,exonends,0,geneid,None,None,None,description])
    GTFReader=staticmethod(GTFReader)
