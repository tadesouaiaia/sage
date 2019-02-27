#!/usr/bin/env python2.7


import os
import sys
import difflib
from collections import defaultdict as dd
import numpy as np
from subprocess import call
import subprocess 
import gzip 
#sys.path.append('/home/rcf-47/souaiaia/analysis/tade/gtfar_source/modules')


##########################################################################################################################################
#####################################################  FASTQ-FILE CLASS START ############################################################
##########################################################################################################################################

def qc_error(eString):
    sys.stderr.write('QCError: '+eString+'\n')
    sys.exit()




class ReadStream:
    def __init__(self,readfile,options): #,strictLevel):
	
	self.options = options 
	self.readlen = options.readlen
	self.readfile = readfile
	self.filename = readfile.name

	self.ext =  readfile.name.split('.')[-1] 


	if self.ext in ['gz']: 
		if readfile.name.split('.')[-2] in ['fastq','fq']: self.ftype = 'fastq' 
		zc = subprocess.Popen(('zcat', self.filename), stdout=subprocess.PIPE)
		output = subprocess.check_output(('wc', '-l'), stdin=zc.stdout)

		self.readfile.close() 
		self.stream = gzip.open(self.filename) 



	else:
		if self.ext in ['fastq','fq']: self.ftype = 'fastq'
		wcl = subprocess.Popen(('wc','-l', self.filename), stdout=subprocess.PIPE)
		output = subprocess.check_output(('head'), stdin = wcl.stdout).split()[0]
		self.stream = readfile 



	if self.ftype == 'fastq': 
		self.total = int(output) / 4 

	else:
		self.total = int(output) 

	self.open = True 


    def read(self):

	name  = self.stream.readline().strip()
	seq   = self.stream.readline().strip()[0:self.readlen] 
	plus  = self.stream.readline().strip()
	qual  = self.stream.readline().strip()[0:self.readlen]

	if len(seq) == 0:    self.open = False 
	elif name[0] != '@': qc_error('Read Format Does not Appear too be valid fastq') 
	else: 		     self.total += 1
        # print name
        # print len(seq)
        # print plus 
        # print qual 


	return name,seq,plus,qual,[]
	

















































class ReadFile:
    def __init__(self,readfile,options,outpath,outprefix,readprefix): #,strictLevel):
        

	self.options = options 

	
	self.prefix  = outpath+'/'+outprefix  
#	self.readprefix  = options.outpath+'/reads/'+options.prefix  
	self.readlen = options.readlen
	self.trimMin = options.readlen
	self.readprefix = readprefix
	self.adaptor_reject_cnt =0 	
	self.adaptor_trim_cnt = 0 
	self.adaptor_found_cnt = 0 

#        self.fileHandle,self.readLen,self.prefix = fileHandle,rlen,prefix

        self.cnt_Entropy,self.cnt_qualMean,self.cnt_qualMax,self.cnt_qualMin=0,0,0,0
	

        #self.BUFFERSIZE=50000

        self.writeTable = {0: open(self.readprefix+'_reject.fastq',"w"), self.readlen: open(self.readprefix+'_1_'+str(self.readlen)+'.fastq',"w")}
	self.writeIdxs = {self.readlen: 1}
	self.writeLines = dd(int) 



 	self.stats = open(self.prefix+".stats","w")
        self.cnt_total,self.cnt_passed,self.cnt_adaptorREJECT = 0,0,0
        self.cnt_interiorN, self.cnt_trimN, self.cnt_lowEntropy, self.dataAdd = 0,0,0,0
        self.cnt_polyA, self.cnt_polyC, self.cnt_polyK = 0,0,0
        self.polyLen=0.66*self.readlen
        
        #self.eval_Entropy,self.eval_qualMean,self.eval_qualMax,self.eval_qualMin = [[]],[[]],[[]],[[]]


        self.polyK_Dict = dd(lambda: dd(int))
        for i in range(1,6): self.polyK_Dict[i]["NA"]=0 
       
	self.readstream = ReadStream(readfile,options)
	self.total = self.readstream.total


 
	self.max_size = self.total 
	self.max_size = self.total / (self.options.cores-1)


	self.cnt_table = dd(int) 
	self.fail_table = dd(int) 
	self.trim_table = dd(int) 
	self.note_table = dd(int) 




 

    def addTrimLengths(self,trims):
        
        MIN_LEN=50
#        try:  trims = sorted([int(s) for s in trimString.split(",") if int(s) < self.readlen and int(s)>=MIN_LEN])
 #       except ValueError:  self.dataFailure("Trim lengths supplied contain non-integer value")
	trims.sort() 
        if len(trims)>0:
            if trims[-1] >= self.readlen: self.dataFailure("Trim lengths supplied contain value equal to or exceeding read length") 
            for t in trims:    
		self.writeTable[t] = open(self.readprefix+"_1_"+str(t)+".fastq","w")
		self.writeIdxs[t] = 1        

	self.trimHash = {} 
        self.trimHash  = dd(bool) 

        trims.append(self.readlen)
        for i in range(len(trims)-1):   self.trimHash.update([(j,trims[i]) for j in range(trims[i],trims[i+1])])
      	for i in range(self.readlen,self.readlen+50):
		self.trimHash[i] = self.readlen  
	self.trimMin = min(trims) 



    def addAdaptors(self,tpAdaptor,fpAdaptor,strictness):

        self.tpAdaptors,self.fpAdaptors = [s for s in tpAdaptor.split(',')],[s for s in fpAdaptor.split(',')]   
        self.tpLens = [len(x) for x in self.tpAdaptors]; self.fpLens = [len(x) for x in self.fpAdaptors]
        if len(self.tpAdaptors) != 2 or len(self.fpAdaptors) != 2:                          self.dataFailure("Adaptors must include two pieces w/ min lengths of 12bp")
        elif min([len(a) for a in self.tpAdaptors]+[len(b) for b in self.fpAdaptors]) < 12: self.dataFailure("Adaptors must include two pieces w/ min lengths of 12bp")
        elif strictness <1 or strictness > 5:                                               self.dataFailure("Supplied strictness must be between 1 and 5")
        else:
            self.minMatch = 3 + strictness
            self.minScr   = [0.8, 0.85, 0.90, 0.925, 0.95, 0.975][strictness]
            self.fpCnt, self.tpCnt, self.doubleCnt, self.fpTrim, self.tpTrim, self.doubleTrim = 0,0,0,0,0,0
            self.fpCnt, self.tpCnt, self.doubleCnt, self.fpTrim, self.tpTrim, self.doubleTrim = 0,0,0,0,0,0
    










    def removeNs(self):



	if 'N' in self.seq: 

            	tmpSeq=sorted(self.seq.split("N"),key=len,reverse=True)[0]
            	tmpQual=self.qual[0:len(tmpSeq)]

		#self.qc_table['N'] = True 
		self.notes.append('N-BASES') 

		self.seq, self.qual = tmpSeq, tmpQual

	return 

    def removeAdaptors(self):


        self.adaptorPass,self.foundAdaptors=False,[]

	fp,tp = 0,len(self.seq) 


        if   self.queryScr(self.fpAdaptors[1],self.fpLens[1]):    fp = self.adaptorLocs[1]
        elif self.queryScr(self.fpAdaptors[0],self.fpLens[0]):    fp = len(self.seq) 

        if self.queryScr(self.tpAdaptors[0],self.tpLens[0]):    tp = self.adaptorLocs[0]
        elif self.queryScr(self.tpAdaptors[1],self.tpLens[1]):  tp = 0 
	

	if fp > 0 or tp < len(self.seq): 

	
		if fp > 0: 		 self.notes.append('FPA')	 
		if tp < len(self.seq):   self.notes.append('TPA') 
		
		self.seq = self.seq[fp:tp] 
		self.qual = self.qual[fp:tp] 

	return 




            


	











    def filter(self):

	#self.name, self.seq, self.plus, self.qual = self.readstream.read() 



	self.cnt_table = dd(int)
	self.qc_table = dd(bool) 
	self.name, fullseq, self.plus, fullqual, self.notes  = self.readstream.read() 

	self.seq,self.qual = fullseq, fullqual 

	

 	while self.readstream.open: 

		self.removeNs() 
		self.removeAdaptors() 


		th = self.trimHash[len(self.seq)]

		notestr = ",".join(self.notes) 
		

		if not th: 
                	self.writeTable[0].write("%s %s\n%s\n%s\n%s\n" % (self.name,notestr,fullseq,self.plus,fullqual))
			self.cnt_table['FAIL'] += 1 
			self.fail_table[notestr]+=1
		else:
			self.seq = self.seq[0:th] 
			self.qual = self.qual[0:th] 
			slen = len(self.seq) 
			self.cnt_table['VALID']+=1
        		self.calculateBasicMetrics()

			if len(self.notes) == 0:
				self.writeTable[slen].write("%s\n%s\n%s\n%s\n" % (self.name, self.seq, self.plus, self.qual)) 
				self.writeLines[slen]   +=1
				self.cnt_table['PASS'] +=1
				if self.writeLines[slen] >= self.max_size and self.writeLines[slen] > 10000:



					self.writeTable[slen].close() 
					self.writeIdxs[slen]+=1
					self.writeTable[slen] = open(self.readprefix+'_'+str(self.writeIdxs[slen])+'_'+str(slen)+'.fastq','w')
					self.writeLines[slen] = 0 


			else:
				self.cnt_table['TRIM']+=1
				self.trim_table[notestr]+=1
				self.writeTable[slen].write("%s %s\n%s\n%s\n%s\n" % (self.name, notestr,self.seq, self.plus, self.qual)) 
		

		for n in self.notes:	self.note_table[n] += 1


		self.name, self.seq, self.plus, self.qual, self.notes = self.readstream.read() 
		



	self.printSummaryStats()












    def queryScr(self,query,qLen,adaptor_name='generic'):
         

        qFind=difflib.SequenceMatcher(None,query,self.seq).find_longest_match(0,qLen,0,len(self.seq))
	 
        if qFind[2] >= self.minMatch:
            seqIdxs   = max(0,qFind[1]-qFind[0]),min(self.readlen,qFind[1]+(qLen-qFind[0]))
            qIdxs = max(0,qFind[0]-qFind[1]),min(qLen,self.readlen-qFind[1]+1)
            myScr = sum(a==b for a,b in zip(self.seq[seqIdxs[0]:seqIdxs[1]],query[qIdxs[0]:qIdxs[1]]))/float(seqIdxs[1]-seqIdxs[0])
            if myScr>self.minScr:
                self.adaptorLocs = seqIdxs
                self.foundAdaptors.append(query)
                return True
        return False        
        







            



    def printSummaryStats(self):
    	


        self.stats.write("TOTAL:RAW-READS %s\n" % self.readstream.total)
        self.stats.write("TOTAL:REJECTED-READS %s\n" %    (self.cnt_table['FAIL']))
        self.stats.write("TOTAL:TRIMMED-READS %s\n" %    (self.cnt_table['TRIM']))
        self.stats.write("TOTAL:UNTRIMMED-READS %s\n" %    (self.cnt_table['PASS']))

	
        self.stats.write("TOTAL:N-READS %s\n" %    (self.note_table['N-BASES']))
        self.stats.write("TOTAL:FPA-READS %s\n" %    (self.note_table['FPA']))
        self.stats.write("TOTAL:TPA-READS %s\n" %    (self.note_table['TPA']))


	VR = float(self.cnt_table['VALID']) 

        self.stats.write("AVG:MEAN-QUALITY %s\n" % (self.cnt_qualMean/VR))
        self.stats.write("AVG:MAX-QUALITY %s\n" % (self.cnt_qualMax/VR))
        self.stats.write("AVG:MIN-QUALITY %s\n" % (self.cnt_qualMin/VR))
      
	return
 
        polyTops = {}
        for i in range(1,6):
            tmpTops=sorted([(self.polyK_Dict[i][x],x) for x in self.polyK_Dict[i].keys()],reverse=True)
            for j in range(5):
                if tmpTops[j][0]>0: self.stats.write("KMER:Read-Dominating-%smer %s %s\n" % (i,tmpTops[j][1],tmpTops[j][0]))
                else:               break


    def calculateBasicMetrics(self):
        self.kMer_Dict=dd(lambda: dd(int))

        self.phredQuals = [ord(x)-35 for x in self.qual]

        for i in range(0,len(self.seq),-4):
            for j in range(1,6):
                self.kMer_Dict[j][self.seq[i:i+j]]+=1
        
        self.cnt_qualMean+=np.mean(self.phredQuals)
        self.cnt_qualMax+=max(self.phredQuals)
        self.cnt_qualMin+=min(self.phredQuals) 
        
        
    def calculateEntropy(self):

        self.kMerObs =   [0]+[len(self.kMer_Dict[x].keys()) for x in range(1,6)]
        self.kMerMaxes = [0]+[max([self.kMer_Dict[x][y] for y in self.kMer_Dict[x].keys()]) for x in range(1,6)]
        

        Pvals = [self.calcObsPvalue(3,self.kMerObs[3]),self.calcObsPvalue(4,self.kMerObs[4]),self.calcObsPvalue(5,self.kMerObs[5])]+[self.calcMaxPvalue(x,self.kMerMaxes[x]) for x in range(1,6)]

        self.seqEntropy = sum(Pvals)/4.0
            



    def calcObsPvalue(self,obs,scr):

        myScrs = self.obs[obs]
        for s in myScrs:
            if scr < s[0]: return s[1]
        return 0.5
    
    def calcMaxPvalue(self,obs,scr):
        myScrs = self.kmax[obs]
        for s in myScrs:
            if scr > s[0]: return s[1]
        return 0.5

    def loadEntropyDist(self):
        
        if self.readlen == 100:
            self.obs,self.kmax={},{}
            self.obs[1] = [(3,0.0),(4,0.00001),(6,0.5)]
            self.obs[2] = [(12,0.0),(13,0.00001),(14,0.0001),(15,0.001),(16,0.01),(17,0.5),(99,0.5)]
            self.obs[3] = [(37,0.0),(38,0.00001),(39,0.0001),(40,0.001),(44,0.01),(46,0.05),(47,0.1),(48,0.2),(49,0.3),(50,0.4),(99,0.5)]
            self.obs[4] = [(62,0.0),(63,0.00001),(64,0.0001),(65,0.001),(69,0.01),(72,0.05),(75,0.1),(76,0.2),(77,0.3),(78,0.4),(99,0.5)]
            self.obs[5] = [(75,0.0),(76,0.00001),(78,0.0001),(80,0.001),(84,0.01),(86,0.05),(87,0.1),(88,0.2),(89,0.3),(90,0.4),(99,0.5)]
            self.kmax[1] = [(44,0.0),(43,0.00001),(42,0.0001),(40,0.001),(37,0.01),(34,0.05),(33,0.1),(32,0.2),(31,0.3),(30,0.4),(0,0.5)]
            self.kmax[2] = [(22,0.0),(21,0.00001),(20,0.0001),(18,0.001),(16,0.01),(14,0.05),(13,0.1),(12,0.2),(11,0.3),(11,0.4),(0,0.5)]
            self.kmax[3] = [(13,0.0),(12,0.00001),(11,0.0001),(10,0.001),(8,0.01),(7,0.05),(6,0.1),(6,0.2),(5,0.3),(5,0.4),(0,0.5)]
            self.kmax[4] = [(11,0.0),(10,0.00001),(9,0.0001),(8,0.001),(7,0.01),(6,0.05),(5,0.1),(4,0.2),(3,0.3),(3,0.4),(0,0.5)]
            self.kmax[5] = [(9,0.0),(8,0.00001),(7,0.0001),(6,0.001),(5,0.01),(4,0.05),(3,0.1),(2,0.2),(2,0.3),(2,0.4),(0,0.5)]
        else:
            self.obs,self.kmax={},{}
            self.obs[1] = [(3,0.0),(4,0.00001),(6,0.5)]
            self.obs[2] = [(7,0.0),(8,0.00001),(10,0.0001),(13,0.001),(15,0.01),(17,0.5),(99,0.5)]
            self.obs[3] = [(25,0.0),(28,0.00001),(29,0.0001),(30,0.001),(34,0.01),(36,0.05),(37,0.1),(38,0.2),(39,0.3),(44,0.4),(99,0.5)]
            self.obs[4] = [(30,0.0),(33,0.00001),(34,0.0001),(35,0.001),(39,0.01),(40,0.05),(40,0.1),(42,0.2),(43,0.3),(44,0.4),(99,0.5)]
            self.obs[5] = [(35,0.0),(36,0.00001),(38,0.0001),(40,0.001),(41,0.01),(42,0.05),(42,0.1),(43,0.2),(43,0.3),(44,0.4),(99,0.5)]
            self.kmax[1] = [(44,0.0),(43,0.00001),(42,0.0001),(40,0.001),(37,0.01),(34,0.05),(33,0.1),(32,0.2),(31,0.3),(30,0.4),(0,0.5)]
            self.kmax[2] = [(22,0.0),(21,0.00001),(20,0.0001),(18,0.001),(16,0.01),(14,0.05),(13,0.1),(12,0.2),(11,0.3),(11,0.4),(0,0.5)]
            self.kmax[3] = [(13,0.0),(12,0.00001),(11,0.0001),(10,0.001),(8,0.01),(7,0.05),(6,0.1),(6,0.2),(5,0.3),(5,0.4),(0,0.5)]
            self.kmax[4] = [(11,0.0),(10,0.00001),(9,0.0001),(8,0.001),(7,0.01),(6,0.05),(5,0.1),(4,0.2),(3,0.3),(3,0.4),(0,0.5)]
            self.kmax[5] = [(9,0.0),(8,0.00001),(7,0.0001),(6,0.001),(5,0.01),(4,0.05),(3,0.1),(2,0.2),(2,0.3),(2,0.4),(0,0.5)]



    def dataFailure(self,msg="FOO"):
        sys.stderr.write(msg+"\n")
        sys.exit(2)
