#!/usr/bin/env python


import os
import sys
import difflib
from collections import defaultdict as dd


#sys.path.append('/home/rcf-47/souaiaia/analysis/tade/gtfar_source/modules')


##########################################################################################################################################
#####################################################  FASTQ-FILE CLASS START ############################################################
##########################################################################################################################################



def parse_error(eString):
    sys.stderr.write('ParseError: '+eString+'\n')
    sys.exit()




class Annotation:
    def __init__(self,chrName,fileHandle,outpath):


        self.isoScores = dd(lambda: [dd(int),dd(int)])
        self.geneScores = dd(lambda: [dd(int),dd(int)])

	self.chr = chrName

	self.outpath = outpath 
	self.samfile = open(outpath+'/sams/full.'+self.chr+'.sam','w')
	self.vertfile = open(outpath+'/verts/full.'+self.chr+'.vert','w')
	self.genefile = open(outpath+'/cnts/'+self.chr+'.genes','w')
	self.isofile = open(outpath+'/cnts/'+self.chr+'.isos','w')


        self.refList,self.jxnKey = [],dd(list)


        for line in open(fileHandle):
            line = line.strip().split()
            start,stop = int(line[0]),int(line[1])
            if line[2] != 'JXN':                  self.refList.append([(start,stop),line[2],line[3]])
            else:                                 self.jxnKey[start] = AnnotatedJxn(line)


        #self.data,self.k,self.fin = [[0,['chrStart']]] + [[int(line.split()[0]),line.split()[1::]] for line in open(fileHandle)]+[[9999999999999,['chrEnd']]]  ,0,False



    def process_mappings(self,map_files): 

 
	for m in map_files:
		self.add_valid_maps(m) 
		while self.open:
			mapping = self.iterate() 
			if mapping.valid:
            			self.samfile.write("%s\n" % ("\t".join(mapping.sam_line)))
            			self.vertfile.write("%s\n" % ("\t".join(mapping.vert_line)))

		
			
		for m,c in mapping.gene_cnt.items():
			self.genefile.write("%s %s\n" % (m,c))
		for m,c in mapping.iso_cnt.items():
			self.isofile.write("%s %s\n" % (m,c))












    def add_valid_maps(self,mapHandle):

        self.mappings = SamMappings(mapHandle)

	

        self.open = True
        self.k=0
        self.refList.sort()
        self.cnt_key = dd(int)


    def iterate(self):

        mapping,m = self.mappings.next(),0

        if not mapping.open:
            self.open = False
            return mapping
        else:



            while mapping.loc > self.refList[self.k][0][1]:  self.k+=1


	    #print mapping.loc, self.k, self.refList[self.k] 
	    
	    
 

            if mapping.loc    < self.refList[self.k][0][0]:

		

                if self.refList[self.k-1][0][-1]+1 == mapping.loc and mapping.loc == self.refList[self.k][0][0]-1 and self.jxnKey[mapping.loc] != []:
		    jk,m = self.jxnKey[mapping.loc], 1

		    

   		    if jk.f_starts == False and jk.f_ends == False:
		    	mapping.add_data(self.refList[k]) 
		    else:
			if jk.f_starts != False:    singleBPexon = "|".join([jG+"="+",".join([x[0]+"@"+x[1] for x in jI]) for jG,jI in jk.f_starts.items()])
			else:  			    singleBPexon = "|".join([jG+"="+",".join([x[0]+"@"+x[1] for x in jI]) for jG,jI in jk.f_ends.items()])
			mapping.add_data([(mapping.loc,mapping.loc),"EXONIC",singleBPexon])
                else:
                    parse_error('smaller?')



            while not mapping.covered:
                mapping.add_data(self.refList[self.k+m])
                m+=1
            mapping.verify()
            mapping.create_lines()
            return mapping












class SamMappings:
    def __init__(self,samFile,verbose=False):

        self.gene_cnt,self.iso_cnt = dd(int),dd(int)
        self.open,self.mapHandle,self.data,self.singular = True,open(samFile),[],False

	#mapLines = self.mapHandle.readlines() 
	#if self.mapHandle.name.split('.')[-2] == 'chr07':

#		for i,x in enumerate(mapLines):
#			print x.split()[0], len(x.split('\t')),i,'wtf'



	#self.mapLines = sorted(mapLines, key = lambda X: int(X.split("\t")[3]), reverse=True)


	#self.mapLines = sorted(self.mapHandle.readlines(), key = lambda x: int(x.split('\t')[3]),reverse=True)
	self.mapLines = sorted(self.mapHandle.readlines(), key = lambda x: int(x.split('\t')[3]),reverse=True)


#	for i in range(3):
#		print self.mapLines[i] 



    def add_data(self,myData):
        self.data.append(myData)
        if self.data[-1][0][1] > self.span[-1]: self.covered = True



    def verify(self):

        if len(self.data) == 1: self.region,self.type,self.subtype,self.features,self.singular = [self.data[0][0]],self.data[0][1],self.data[0][1],self.data[0][2],True
        else:
            self.singular, gCnt,iCnt,self.overlap= False,dd(int),dd(int),False
            self.region,types,features = [d[0] for d in self.data],[d[1] for d in self.data], [d[2] for d in self.data]
            self.type,self.subtype = "/".join([str(x) for x in sorted(set([d[1] for d in self.data]))]),"/".join([str(x) for x in sorted(set([d[1] for d in self.data]))])

            for i,F in enumerate(features):
                for g in list(set([a.split("=")[0] for b in [f.split(">") for f in F.split("|")] for a in b])): gCnt[g] +=1
                for g in list(set([a for b in [f.split(">") for f in F.split("|")] for a in b if a.split("=")[-1]!='NA'])): 
                    for gf in  g.split("=")[-1].split(","): iCnt[g.split("=")[0]+"@"+gf]+=1
                    #iCnt[g]+=1
            gPass,iPass = [gD for gD,gC in gCnt.items() if gC == len(self.data)], [gD for gD,gC in iCnt.items() if gC == len(self.data)]
            if len(iPass)>0:    self.features = '|'.join(iPass)
            elif len(gPass)>0:  self.features = '|'.join(gPass)
            else:               self.features = 'NA'
            if self.type == 'EXONIC' and len(iPass)>0:  self.subtype = 'EXONIC-SPECIFIC'
            elif 'EXONIC' in types:                     self.subtype = 'EXONIC-OVERLAP'
            if 'INTERGENIC' in types:                    self.type ='INTERGENIC'
            elif 'INTRONIC' in types or len(iPass) == 0: self.type = 'INTRONIC'


        self.valid = True
        self.region = ",".join([str(s[0])+'-'+str(s[1]) for s in self.region])




	

    def next(self):
        self.data,self.covered,self.valid = [],False,False


	try: 
		self.line = self.mapLines.pop() 
            	self.sam_start = self.line.split()
            	s = self.line.strip().split()
            	self.id,self.strand,self.chr,self.loc,self.cigar,self.seq,self.places,self.subs,self.rank = s[0],s[1],s[2],int(s[3]),s[5],s[9],s[11],s[12],s[13]
            	self.locations = int(self.places.split(":")[-1])
            	self.span      = (self.loc,self.loc+len(self.seq)-1)


            	return self

	except IndexError:
		self.open = False
		return self



    def create_lines(self):

        gene_dict = dd(list)

       # EXONIC EXONIC RALY2=T1@9,T2@10
       # INTRONIC EXONIC-OVERLAP SNRPB@T1@4|SNRPB@T201@6|SNRPB@T2@4|SNRPB@T3@4

        if self.chr == 'chrF':
            self.gene_cnt[self.chr+":"+self.type+":"+self.features]+=1
            self.valid = False
            return 





        if self.type == 'EXONIC' and self.subtype == 'EXONIC':
            for f in self.features.split("|"):
                gene_dict[f.split("=")[0]].extend(f.split("=")[1].split(","))

        elif self.type == 'EXONIC' and self.subtype == 'EXONIC-SPECIFIC':
            for f in self.features.split("|"):  gene_dict[f.split("@")[0]].append("@".join(f.split("@")[1::]))

        elif self.type == 'INTRONIC' and self.subtype == 'EXONIC-OVERLAP':
            for f in self.features.split("|"):  gene_dict[f.split("@")[0]].append("@".join(f.split("@")[1::]))

        elif self.type == 'INTRONIC' and self.subtype == 'INTRONIC':
            genes = "/".join([f.split("=")[0] for f in self.features.split("|")])
            trans = self.features


        elif self.type == 'INTERGENIC':    
            genes,trans = self.features,self.features

        else:
            print self.type,self.subtype,self.features


        if len(gene_dict)>0:
            genes = "/".join(gene_dict.keys())
            trans = "/".join([gene+"="+"|".join(trans) for gene,trans in gene_dict.items()])

        self.vert_line = [self.region,trans,self.type,self.subtype,"-".join([str(s) for s in self.span]),self.places]
	
        MT='MT:Z:'+self.type
        GE='GE:Z:'+genes
        TR="TR:Z:"+trans
        FT="FT:Z:"+self.subtype
        AD="AD:Z:NA"


        self.sam_line = self.sam_start+[MT,GE,TR,FT,AD]


	if self.type == 'EXONIC':
		multi_trans = trans.split('/') 
		if len(multi_trans) > 1:
			m_valid = [mt for mt in multi_trans if mt.split('=')[-1] != 'NA'] 
			if len(m_valid) > 0 and len(m_valid) < len(multi_trans): 
				genes = "/".join([mv.split('=')[0] for mv in m_valid]) 

        self.gene_cnt[self.chr+":"+self.type+":"+genes]+=1
        self.iso_cnt[self.chr+":"+self.type+":"+trans]+=1
        



















class AnnotatedJxn:
    def __init__(self,annoLine,verbose=False):

        self.gene_starts = annoLine[4].split(",")
        self.gene_ends   = annoLine[5].split(",")
        self.arrival_key,self.departure_key = {},{}
        self.f_starts = self.create_key(annoLine[6])
        self.f_ends = self.create_key(annoLine[7])

        for arrival in annoLine[8].split("&"):
            if arrival != 'NA':
                self.arrival_key[int(arrival.split('>')[0])] = self.create_key(arrival.split(">")[1])
        for departure in annoLine[9].split("&"):
            if departure != 'NA':
                self.departure_key[int(departure.split('>')[0])] = self.create_key(departure.split(">")[1])




    def create_key(self,kStr):

        jKey = {}
        if kStr != 'NA':
            for gene,isoforms in [grp.split("=") for grp in kStr.split("|")]: 
	    	jKey[gene] = [isoform.split("@") for isoform in isoforms.split(",")]
            return jKey
        return False











#!/usr/bin/env python



#sys.path.append('/home/rcf-47/souaiaia/analysis/tade/gtfar_source/modules')


##########################################################################################################################################
#####################################################  FASTQ-FILE CLASS START ############################################################
##########################################################################################################################################



def parse_error(eString):
    sys.stderr.write('ParseError: '+eString+'\n')
    sys.exit()




class NovelAnnotation:
    def __init__(self,chrName,fileHandle,outpath,readlen=100):


        self.readlen = readlen
        self.refList,self.jxnKey = [],dd(list)
        for line in open(fileHandle):
            line = line.strip().split()
            start,stop = int(line[0]),int(line[1])
            if line[2] == 'JXN' and line[9] != 'NA':
                for jxn in line[9].split("&"):
                    dest = int(jxn.split(">")[0])
                    self.jxnKey[start].append(dest);self.jxnKey[dest].append(start) 
            elif line[2] != 'JXN':  self.refList.append([(start,stop),line[2],line[3]])
 

        self.refList.sort()
	self.chr = chrName
	self.candfile = open(outpath+'/cands/'+self.chr+'.cands','w')
	self.isofile = open(outpath+'/cnts/'+self.chr+'.novel.isos','w')


    def add_ref(self,filePath):

        self.seq = []
        c=open(filePath)
        self.chr = c.readline().strip().split()[0].split(">")[1]
        for line in c: self.seq.extend([s for s in line.strip()])
        c.close()
	return self

    def add_novel_maps(self,mapHandle):
        self.clip_mappings = ClipMappings(mapHandle)
        self.open = True



    def process_mappings(self,mappings): 
	 
	for m in mappings:
		self.add_novel_maps(m) 
		while self.open:
			cands = self.iterate()
			for header,seq,cType,xDist,yDist,cCnt in cands:
        			self.candfile.write('%s\n%s\n'%('>'+header,seq))
        			self.isofile.write('%s %s\n'%('>'+header,cCnt))




















    def iterate(self):

	output = [] 
        readlen = self.readlen
        extlen  = self.readlen -3
        cands,cA,cT,siteKey,k=[],set([]),{},{},0
        while self.clip_mappings.open:
            mapping,m = self.clip_mappings.next(),0
            if not mapping:
                self.open = False
                break
            seq1 = "".join(self.seq[mapping.loc-1:mapping.loc+self.readlen-1])
            seq2 = "".join(self.seq[mapping.loc+mapping.cigar[1]-1:mapping.loc+mapping.cigar[1]+mapping.readlen-1])

            for cand,candType in mapping.compare_seqs(seq1,seq2):

                if cand.jxn[1] not in self.jxnKey[cand.jxn[0]]:
                    cands.append(cand); cA.add(cand.jxn[0]); cA.add(cand.jxn[1]); cT[cand] = candType
        for site in sorted(list(cA)):
            while self.refList[k][0][1] < site: k+=1
            siteKey[site] = self.refList[k]



        cnt_key = dd(int)
        seq_key = {}
        left_set,right_set = set([]),set([])
        for cand in cands:

            ct = cT[cand]
            x,y = cand.jxn
            sX,sY,tX,tY=siteKey[x],siteKey[y],siteKey[x][1],siteKey[y][1]
            if tX == 'INTERGENIC':      gX = [sX[-1].split(">")[-1]]
            else:                       gX = list(set([s.split("=")[0] for s in sX[-1].split("|")]))
            if tY == 'INTERGENIC':      gY = [sY[-1].split(">")[0]]
            else:                       gY = list(set([s.split("=")[0] for s in sY[-1].split("|")]))

            if sX[0][0] > x-extlen: jX = [sX[0][0],x]
            else:                    jX = [x-extlen,x]
            if sY[0][1] < y+extlen: jY = [y,sY[0][1]]
            else:                    jY = [y,y+extlen]

            sharedGenes = [gene for gene in gX if gene in gY]
            if len(sharedGenes)>0: sGenes = ":".join(sharedGenes)
            else:                  sGenes = 'NA'
            tGenes = "|".join(gX+gY)
            if (jX[1]-jX[0]) + (jY[1]-jY[0]) < readlen: continue



            header =  self.chr+":NOVEL:"+str(jX[0])+"-"+str(jX[1])+'-'+str(jY[0])+'-'+str(jY[1])+":"+tX+":"+tY+":"+sGenes+"|"+tGenes
            seqX = "".join(self.seq[jX[0]-1:jX[1]])
            seqY = "".join(self.seq[jY[0]-1:jY[1]])
            total_seq = seqX+seqY

            runScr,dubScr,tripScr = self.calc_entropy(total_seq)

	    if runScr>100 or dubScr > 0.20 or tripScr > 0.10: continue


            cnt_key[(x,y)]+=1
            seq_key[(x,y)] = [header,total_seq,ct]
            left_set.add(x)
            right_set.add(y) 

        dist_key=dd(int)
        left_list = sorted(list(set(left_set)))
        right_list = sorted(list(set(right_set)))


        for i,s in enumerate(left_list):
            if i == len(left_list)-1:dist_key[s] = s-left_list[i-1]
            elif i == 0:    dist_key[s] = left_list[i+1]-s
            else: dist_key[s] = min(s-left_list[i-1],left_list[i+1]-s)

        for i,s in enumerate(right_list):
            if i == len(right_list)-1:dist_key[s] = s-right_list[i-1]
            elif i == 0:    dist_key[s] = right_list[i+1]-s
            else: dist_key[s] = min(s-right_list[i-1],right_list[i+1]-s)


        for (x,y),(header,seq,sType) in seq_key.items():
#            if sType == "CANON" or (dist_key[x]>10 and dist_key[y]>10 and cnt_key[(x,y)]>1):	
		output.append((header,seq,sType,dist_key[x],dist_key[y],cnt_key[(x,y)]))	
        return output  
		

        #sys.stdout.write('%s\n%s\n'%(header,seq))


    def calc_entropy(self,seq):

        s,hRun,hiRuns,x=seq[0],0,[],1
        kMers = dd(int)
        kTrips = dd(int)
        for i in range(1,len(seq)):
            kMers[seq[i:i+2]]+=1
            kTrips[seq[i:i+3]]+=1
            if seq[i] == s: hRun+=1
            else:
                if hRun > 3: hiRuns.append(hRun) 
                s,hRun = seq[i],0
    
        kDub,kTrip = max(kMers.values()),max(kTrips.values())
        for h in hiRuns: x*=h
        sl = float(len(seq))
        return x/sl,kDub/sl,kTrip/sl



class ClipMappings:
    def __init__(self,samFile,verbose=False):

        self.open,self.mapHandle,self.data,self.singular = True,open(samFile),[],False

	'''
    def add_data(self,myData):
        if self.map_type == 'continuous':
            self.data.append(myData)
            if self.data[-1][0][1] > self.span[-1]: self.covered = True
        else:
            print self.map_type
            print 'hmmm'
            sys.exit()
	'''


    def next(self):
        self.readlen=100
        self.data,self.covered,self.line = [],False,self.mapHandle.readline().strip()
        if len(self.line) == 0:
            self.open = False
            return False
        else:
            sp = self.line.split()
            self.id,self.loc,self.seq,self.chr,cigar= sp[0],int(sp[3]),sp[9],sp[3],sp[5].split("M")
            self.cigar =[int(m) for m in [cigar[0]]+cigar[1].split("N")]
            return self

    def compare_seqs(self,seq1,seq2):
            topLocs,nextLocs = [],[]
            fq1,fq2 = self.sub_cnt(self.seq[0:25],seq1[0:25]),self.sub_cnt(self.seq[0:25],seq2[0:25])
            tq1,tq2 = self.sub_cnt(self.seq[75::],seq1[75::]),self.sub_cnt(self.seq[75::],seq2[75::])
#            if fq1 + tq2 == 0 and fq2 > 10 and tq1 > 10:

	    if fq1 + tq2 > -1000: 
                pivotScores,firstSubs,secondSubs = [],[],[]
                for i in range(25,76):
                    if self.seq[i] != seq1[i]:  firstSubs.append(i)
                    if self.seq[i] != seq2[i]: secondSubs.append(i) 
                for i in range(25,76):  pivotScores.append([len([f for f in firstSubs if f<=i])+len([f for f in secondSubs if f>=i]),i])
        #        pivotPos = sorted([p[1] for p in pivotScores if p[0] == 0])
                pivotPos = sorted([p[1] for p in pivotScores if p[0] < 10])
	
		
                if len(pivotPos)>0:
	
                    for p in pivotPos:
                        jSeq = seq1[p]+seq1[p+1]+seq2[p-2]+seq2[p-1]
                        endBreak,startBreak = self.loc+p-1,self.loc+p+self.cigar[1]
                        startTuple,endTuple = (self.loc,self.loc+p-1),(startBreak,startBreak+self.readlen-p-1)
                        myseq1  = seq1[0:p]
                        myseq2  = seq2[p::]


                        if len(set([jx for jx in  jSeq]))>1:


                            if jSeq == 'CTAC' or jSeq == 'GTAG':                    topLocs.append((ClipCand(startTuple,endTuple,myseq1,myseq2,jSeq),"CANON"))
                            else:
                                mix = myseq1[-2::]+myseq2[0:2]
                                if mix[1] != mix[2] and len(set([x for x in mix]))>2 :
                                    nextLocs.append((ClipCand(startTuple,endTuple,myseq1,myseq2,jSeq),"NONCANON"))

            if len(topLocs)>0: return topLocs
            else:              return nextLocs



    def sub_cnt(self,mySeq,xSeq):

        if mySeq == xSeq: return 0
        return len([1 for i in range(len(mySeq)) if mySeq[i]!=xSeq[i]])


class ClipCand:
    def __init__(self,startTuple,endTuple,myseq1,myseq2,jSeq):

        self.splice_acceptor = jSeq
        self.siteA, self.siteB = startTuple,endTuple
        self.jxn = (startTuple[-1],endTuple[0])
        self.seq1,self.seq2=myseq1,myseq2















