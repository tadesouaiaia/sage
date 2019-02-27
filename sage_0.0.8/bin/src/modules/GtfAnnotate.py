#!/usr/bin/env python

import os
import sys
from collections import defaultdict as dd
import string

##########################################################################################################################################
##################################################  GTF FILE - CLASS START  ##############################################################
##########################################################################################################################################



def geneJoin(myList,sep=','):
    if len(myList) == 0: return 'NA'
    else:                return sep.join([str(s) for s in myList]) 



def gtf_error(eString):
    sys.stderr.write('GtfError: '+eString+'\n')
    sys.exit()




class GtfFile:
    def __init__(self, gtf, chrPath, options, readlen=100, flanklen = 1000 ):


	

        self.finished,self.j_files = {},{}

        self.chrPath,self.readLen, self.flanklen =  chrPath,readlen, flanklen
        self.open, self.minLen, self.maxLen,self.line = True, 35, 400,'#'
        self.genes = GtfGenes(gtf,self.readLen)

	if self.chrPath[-1] == '/': self.chrPath = self.chrPath[0:-1]
	

	self.options = options

	if self.options.prefix == None: self.options.prefix = 'sagePre'
	self.listPrefix = self.options.outpath+'/lists'
        self.refPrefix =  self.options.outpath+'/refs/'+self.options.prefix
        self.annoPrefix = self.options.outpath+'/anno/'+self.options.prefix
	self.statPrefix = self.options.outpath+'/stats/'+self.options.prefix

############################################################################################################################################
######################################## Annotation Methods ################################################################################
############################################################################################################################################




    def printGeneNicks(self):

        self.name_key = open(self.annoPrefix+'.'+self.chr+'.feature.key','w')
        for gene,nick in self.genes.gene_nicks.items():
            self.genes.reverse_nicks[nick] = gene 
            self.name_key.write("%s %s %s %s\n" % (self.chr,'GENE',nick,gene))
            for t,tv in self.genes.tran_nicks[nick].items():    self.name_key.write("%s %s %s %s %s\n" % (self.chr,'ISOFORM',nick,tv,t))
        self.name_key.close()


    def printGeneSeqs(self):
	u,v,w = open(self.statPrefix+"_"+self.chr+'_gene_stats.txt','w'), open(self.statPrefix+"_"+self.chr+'_transcript_reads.fa','w'), open(self.statPrefix+"_"+self.chr+'_transcript_seqs.fa','w')
	u.write('%s %s %s %s %s\n' % ('---','transcripts','uniq-exon-bases','uniq-mapping-bases','sim-reads'))
	trapLen = self.readLen + 5
	for n,gene in enumerate(self.genes.members):	
		gFull = "~".join(self.genes.reverse_nicks[gene.name].split('@'))
		r_total, baselens, exlens, tran_seqs, x_obs, j_obs = 0,0, 0, dd(list) , dd(bool), dd(bool)




		'''

		t_refs =  {T: "".join(["".join(self.seq[a:b]) for (a,b) in E]) for T,E in gene.exons.items()}

		for t,sR in t_refs.items():
			if len(sR) > 2*self.readLen: 
				w.write('>~%s~\n%s\n' % ("~".join([self.chr,str(n+1),gFull,t]),sR))
				exlens += len(sR) 
				if len(sR) < 500: stepsize = 10 
				elif len(sR) < 2000: stepsize = 20 
				else:		     stepsize = 50 
				

				for k,i in enumerate(range(10,len(sR)-self.readLen,stepsize)):
					s = sR[i:i+self.readLen]
					baselens += len(s) 
					r_total += 1
					v.write('>~%s~\n%s\n' % ( "~".join([self.chr,str(n+1),gFull,t,str(k+1)]),s))
		u.write('~%s~ %s %s %s %s\n' % ("~".join([self.chr,str(n+1),gFull]),len(tran_seqs.keys()),exlens,baselens,r_total))
		continue 
		'''


 
		for T,E in gene.exons.items():
			if len(E) == 1: 
				eS,eE = E[0] 
				if x_obs[(eS,eE)]: continue 
				x_obs[(eS,eE)] = True 
				tran_seqs[T].append("".join(self.seq[eS:eE]))
			else:
				if not x_obs[(E[0][0],E[0][1])]: 		myseq = [E[0]]
				elif not j_obs[(E[0][1],E[1][0])]:
					if E[0][1]-E[0][0]  <= self.readLen: 	myseq = [E[0]]
					else:			      		myseq = [[E[0][1]-self.readLen,E[0][1]]]
				else:						myseq = [] 
				for i,e in [Xt for Xt in enumerate(E[1:-1])][1::]:
					pB,nA,eA,eB = E[i-1][1],E[i+1][0],e[0],e[1]
					if not x_obs[(e[0],e[1])]: 
						myseq.append(e) 
						x_obs[(eA,eB)] = True 
						continue	
					elif len(myseq) == 0:
						if j_obs[(eB,nA)]: continue 
						elif eB-eA <= self.readLen: myseq = [e] 
						else:			    myseq = [[eB-self.readLen,eB]]
					else:
						if eB-eA <= self.readLen * 2: 
							if (not j_obs[(pB,eA)]) or (not j_obs[(eB,nA)]): 
								myseq.append(e) 
							else:
								tran_seqs[T].append(["".join(self.seq[start:stop]) for start,stop in myseq])
								myseq = [] 
						else:
							if (not j_obs[(pB,eA)]): myseq.append([eA,eA+self.readLen])
							tran_seqs[T].append(["".join(self.seq[start:stop]) for start,stop in myseq])
							if (not j_obs[(eB,nA)]): myseq = [[eB-self.readLen,eB]]
							else:			 myseq = [] 
					j_obs[(pB,eA)] = True 
				if (not x_obs[(E[-1][0],E[-1][1])]): myseq.append(E[-1])				
				elif len(myseq) != 0: 
					if (not j_obs[(E[-2][1],E[-1][0])]):
						if E[-1][1]-E[-1][0] < self.readLen: myseq.append(E[-1]) 
						else: 					  myseq.append([E[-1][0],E[-1][0]+self.readLen])
						j_obs[(E[-2][1],E[-1][0])] = True 
				x_obs[(E[-1][0],E[-1][1])] = True
				if len(myseq) > 0: tran_seqs[T].append("".join(["".join(self.seq[start:stop]) for start,stop in myseq]))
				


		t_refs =  {T: "".join(["".join(self.seq[a:b]) for (a,b) in E]) for T,E in gene.exons.items()}




		for t,sR in t_refs.items():
			if len(sR) > 2*self.readLen: 
				w.write('>~%s~\n%s\n' % ("~".join([self.chr,str(n+1),gFull,t]),sR))

				exlens += len(sR) 
				SEQS = tran_seqs[t] 	
				for z,S in enumerate(SEQS):
					baselens += len(S) 
					if len(S) <= self.readLen+10: continue 
					stepsize = int(len(S)/100.0)
					if stepsize   > 25: stepsize = 25
					elif stepsize < 10: stepsize = 10
					for k,i in enumerate(range(5,len(S)-self.readLen,stepsize)):
						r_total += 1
						s = S[i:i+self.readLen]
						v.write('>~%s~\n%s\n' % ( "~".join([self.chr,str(n+1),gFull,t,str(z+1),str(k+1)]),s))
			#w.write('>%s\n%s\n' % ("~".join([self.chr,str(n+1),gFull,t,str(z+1)]),S))

#		uniq_bases = sum([E[1]-E[0] for E in list(set([a for b in [[tuple(x) for x in X] for X in gene.exons.values()] for a in b]))])
		u.write('~%s~ %s %s %s %s\n' % ("~".join([self.chr,str(n+1),gFull]),len(tran_seqs.keys()),exlens,baselens,r_total))


		'''
		for t,SEQS in tran_seqs.items():
			for z,S in enumerate(SEQS):
				baselens += len(S) 
				if len(S) <= self.readLen+1: continue 
				stepsize = int(len(S)/100.0)
				if stepsize   > 25: stepsize = 25
				elif stepsize < 10: stepsize = 10
				for k,i in enumerate(range(0,len(S)-self.readLen,stepsize)):
					r_total += 1
					s = S[i:i+self.readLen]
					v.write('>%s\n%s\n' % ( "~".join([self.chr,str(n+1),gFull,t,str(z+1),str(k+1)]),s))
				#w.write('>%s\n%s\n' % ("~".join([self.chr,str(n+1),gFull,t,str(z+1)]),S))

		uniq_bases = sum([E[1]-E[0] for E in list(set([a for b in [[tuple(x) for x in X] for X in gene.exons.values()] for a in b]))])
		u.write('%s %s %s %s %s\n' % ("~".join([self.chr,str(n+1),gFull]),len(tran_seqs.keys()),uniq_bases,baselens,r_total))
		'''

			




    def getGeneArrivalsAndDeparatures(self,jumps):

	arrivals, departures = dd(list), dd(list)
    	for gene in jumps[0]:
        	for loc,tran in jumps[0][gene]:
                	trans = [t[0]+'@'+str(t[1])+'-'+str(t[2]) for t in tran]
                    	arrivals[loc].append(gene+'='+geneJoin(trans))
        arriveStr = geneJoin([str(a)+">"+geneJoin(arrivals[a],"|") for a in arrivals],'&')

        for gene in jumps[1]:
        	for loc,tran in jumps[1][gene]:
                	trans = [t[0]+'@'+str(t[1])+'-'+str(t[2]) for t in tran]
                	departures[loc].append(gene+'='+geneJoin(trans))
        depStr = geneJoin([str(a)+">"+geneJoin(departures[a],"|") for a in departures],'&')
	return arriveStr,depStr




    def printData(self):

	if not self.genes.open:
		self.open = False 
		return 
        self.chr,jumps,starts,ends,jxns,gene_info = self.genes.load()

        self.addFasta(self.chrPath+'/'+self.chr+'.fa')
	self.printGeneNicks()
	if self.options.print_gene_seqs:	self.printGeneSeqs() 
        self.finished[self.chr] = self.chrPath+"/"+self.chr+'.fa'
        self.anno_file = open(self.annoPrefix+ "."+self.chr+".annotation",'w')
        spots = sorted(set(starts.keys()+ends.keys()))
        lastSpot,lastGenes,iMap,lastIsos =  -1,['NA'],dd(list),dd(list)



        for i in spots:


#            print 'hi', i, iMap.keys() 


            geneStarts,geneEnds,lastGenes,sMap,eMap = [gene for gene in starts[i].keys() if gene not in iMap.keys()],[gene for gene in ends[i].keys() if gene_info[gene][1] == i],[],dd(list),dd(list)
            for gene in starts[i].keys():
                for isoList in starts[i][gene]:
                    for iso,n in isoList:  sMap[gene].append(iso+'@'+str(n))
            for gene in ends[i].keys():
                for isoList in ends[i][gene]:
                    for iso,n in isoList:  eMap[gene].append(iso+'@'+str(n))

            isoStarts,isoEnds = geneJoin([g+'='+geneJoin(isos) for g,isos in sMap.items()],"|"),geneJoin([g+'='+geneJoin(isos) for g,isos in eMap.items()],"|"),
    	    arriveStr, depStr = self.getGeneArrivalsAndDeparatures(jumps[i])


            if len(iMap.keys()) == 0:
                spData = [lastSpot+1,i-1,'INTERGENIC',geneJoin(sorted(lastGenes),"|")+">"+geneJoin(sorted(sMap.keys()),"|")]
            elif len([a for b in iMap.values() for a in b]) == 0: 
                spData = [lastSpot+1,i-1,'INTRONIC',geneJoin(sorted([g+'='+geneJoin(sorted(lastIsos[g])) for g in sorted(lastIsos)]),"|")+">"+geneJoin(sorted([g+'='+geneJoin(sorted(sMap[g])) for g in sorted(sMap.keys())]),"|")]
            else:                                                 spData = [lastSpot,i,'EXONIC',geneJoin(sorted([g+'='+geneJoin(sorted(iMap[g])) for g in sorted(set(iMap.keys()))]),"|")]




            self.anno_file.write('%s\n' % "\t".join([str(X) for X in spData]))
            spotData = [i,i,'JXN','GS,GE,IS,IE,Arr,Dep',geneJoin(geneStarts),geneJoin(geneEnds),isoStarts,isoEnds,arriveStr,depStr]
            self.anno_file.write('%s\n' % "\t".join([str(X) for X in spotData]))

#	    print 'sp'," ".join([str(X) for X in spData])
#	    print 'spot'," ".join([str(X) for X in spotData]) 
#	    print ""


            for gene in set(sMap.keys()+eMap.keys()): iMap[gene] = sorted(list(set(sMap[gene]+eMap[gene]+iMap[gene])))
            for gene,isos in iMap.items():





                if i == gene_info[gene][1]:
                    lastGenes.append(gene)
                    iMap.pop(gene)
                else:
                    if gene in lastIsos.keys(): keyL = {I.split("@")[0]:I.split("@")[1] for I in lastIsos[gene]}
                    else:                       keyL = {} 
                    keyI,keyS,keyE = {I.split("@")[0]:I.split("@")[1] for I in isos},{I.split("@")[0]:I.split("@")[1] for I in sMap[gene]},{I.split("@")[0]:I.split("@")[1] for I in eMap[gene]}

                    for tName,tCnt in keyI.items():
                        if tName in keyE.keys():
                            if keyI[tName] != keyE[tName]: 


				gtf_error('Invalid Naming') 
                            else: 			   
				keyL[tName] = keyI.pop(tName)


                        elif tName in keyS.keys():
                            if tName in keyL.keys(): keyL.pop(tName)
                    iMap[gene] = [tName+'@'+tCnt for tName,tCnt in keyI.items()]
                    lastIsos[gene] = [tName+'@'+tCnt for tName,tCnt in keyL.items()]
            lastIsos = {gene: lastIsos[gene] for gene in lastIsos.keys() if len(lastIsos[gene])>0 and gene not in lastGenes}
            sMap,eMap,lastSpot = dd(list),dd(list),i
        self.anno_file.write('%s\n' % ("\t".join([str(X) for X in [lastSpot+1,lastSpot*100,'INTERGENIC',geneJoin(lastGenes,"|")+">"+geneJoin([],"|")]])))

        if len(jxns.keys())>0:
            self.jxn_file = open(self.refPrefix+"_"+self.chr+"_known_jxns.fa",'w') 
            for jxn,genes in jxns.items():
                seq = "".join(["".join(self.seq[s[0]-1:s[1]]) for s in jxn])
                if len(seq) < self.readLen: continue
                ref,ref2 = '>'+self.chr+':KJXN='+"-".join([str(s[0])+'-'+str(s[1]) for s in jxn]),[]
                for gene,transcripts in genes.items():  ref2.append(gene+'='+"|".join([t[0]+':'+t[1] for t in [a for b in transcripts for a in b]]))
                sName = ref+':'+"/".join(ref2)
                self.jxn_file.write('%s\n%s\n' % (sName,seq))
            self.j_files[self.chr] =  self.refPrefix+"_"+self.chr+"_known_jxns.fa"
            self.jxn_file.close()
        self.anno_file.close()










    def print_unannotated_chrs(self,filterDir):

        from shutil import copyfile
        if filterDir:
            #fDir = os.path.dirname(__file__)+'/'+filterDir
            for F in os.listdir(filterDir):
                f,fSuffix = F.split(".")[0],F.split(".")[-1]
                if fSuffix in ['fa','fasta'] and f not in self.finished.keys():
                    if f+'.annotation' in os.listdir(fDir):
                        copyfile(fDir+'/'+f+'.annotation',self.annoPrefix+'.'+f+'.annotation')
                        self.finished[f] = fDir+'/'+F
       # print os.path(self.chrPath)
        for F in os.listdir(self.chrPath):
            f,fSuffix = F.split(".")[0],F.split(".")[-1]
            if fSuffix in ['fa','fasta'] and f not in self.finished.keys():
                    self.finished[f] = self.chrPath+"/"+F
                    self.anno_file = open(self.annoPrefix+ "."+f+".annotation",'w')
                    self.anno_file.write('%s %s INTERGENIC NA>NA\n' % (0,9999999999))
                    self.anno_file.close()

        c_file = open(self.listPrefix+"/CHRS.txt",'w')
        j_file = open(self.listPrefix+"/JXNS.txt",'w')
        c_file.write("\n".join(self.finished.values())+"\n")
        j_file.write("\n".join(self.j_files.values())+"\n")





    def addFasta(self,filePath):
        self.seq = []
        try:
            c=open(filePath)
        except IOError:
            if len(filePath.split("/")) > 1:
                failChr = filePath.split("/")[-1]
            else: 
		failChr = filePath
            gtf_error("Error: "+failChr+" is referenced by the suppled gtf-file but not found in the supplied chromosome path: "+filePath)


        fChr = c.readline().strip().split()[0].split(">")[1]
        if self.chr != fChr:
            print self.chr,fChr
            print "Wrong Chromosome File Error"; sys.exit()
        else:
            for line in c:
                self.seq.extend([s for s in line.strip()])
        c.close()



class GtfGenes:

    def __init__(self, gtfFile,readlen, flanklen = 1000 ):


        self.t_dict = dd(int)

        self.GTF_STYLE = "NA"
	self.members = [] 
        self.filehandle,self.readLen,self.chr,self.line = gtfFile,readlen,None,["#"]

	self.check_gtf_file()

        self.gene_nicks, self.tran_nicks = {},dd(lambda: {}) 
	self.reverse_nicks = {} 



    def check_gtf_file(self):


	line = self.filehandle.readline().strip() 
	while line[0] == '#': line = self.filehandle.readline().strip() 
	self.chr,self.open = line.split()[0], True 
	if len(line.split('gene_description'))> 1:	
		self.getGeneFromFile = self.getGeneFromRhesusFile					
		self.GTF_STYLE = 'RHESUS'
	else:	
		self.getGeneFromFile = self.getGeneFromGencodeFile					
		self.GTF_STYLE = 'HUMAN'
	self.line = line.split() 





    def apply_genenick(self,gID,hID):
        hShort=hID.split(".")[0].split("-")[0].split("_")[0]
        hLong = gID+'@'+hID

	
        if hLong in self.gene_nicks.keys(): return self.gene_nicks[hLong]
        elif hShort not in self.gene_nicks.values():
            self.gene_nicks[hLong] = hShort
            return hShort
        else:
            try:
                hTmp,hEnd,k = hShort,int(hShort[-1]),1
                myLetters = string.uppercase
                while k<len(myLetters) and hTmp in self.gene_nicks.values():
                    hTmp= hShort + myLetters[k]
                    k+=1
                hShort = hTmp
            except ValueError: END=False

            hTmp,k = hShort,2
            while hTmp in self.gene_nicks.values():
                hTmp = hShort+str(k)
                k+=1
            self.gene_nicks[hLong] = hTmp
            return hTmp

    def apply_trannick(self,gID,tID,hID):

        tShort =  hID.split("-")[-1]
        try:                tShort = 'T'+str(int(tShort))
        except ValueError:   tShort = 'T'+tShort
        tLong = tID+'@'+hID
        if tLong in self.tran_nicks[gID].keys():
            return self.tran_nicks[gID][tLong]
        elif tShort not in self.tran_nicks[gID].values():
            self.tran_nicks[gID][tLong] = tShort
            return tShort
        else:
            hTmp,k = tShort,2
            while hTmp in self.tran_nicks[gID].values():
                hTmp = tShort+'-'+str(k)
                k+=1
            self.tran_nicks[gID][tLong] = hTmp
            return hTmp




    def load(self):

        self.info, self.members = dd(list), [] 
        self.chr_jxns,self.starts,self.ends, self.jumps = dd(lambda: dd(list)),dd(lambda: dd(list)),dd(lambda: dd(list)),  dd(lambda: [dd(list),dd(list)])
        #self.breaks, self.jumps = dd(lambda: dd(list)), dd(lambda: [dd(list),dd(list)])
#        self.breaks, self.jumps = dd(lambda: dd(list)), dd(lambda: [dd(list),dd(list)])
	self.gene_nicks, self.tran_nicks = {},dd(lambda: {}) 
	if not self.open: 
		self.chr = None 
		return None 

	self.prev_chr = self.chr 
	while True:


	    gene = self.getGeneFromFile() 



	    if gene.empty:  break 
            starts,stops,jumps,valid_jxns = gene.update(self.readLen)
            for jxn,item in valid_jxns.items(): self.chr_jxns[jxn][gene.name].append(item) 
            for start,item in starts.items(): self.starts[start][gene.name].append(item)
            for end,item in stops.items(): self.ends[end][gene.name].append(item)

            for (a,b) in jumps:
                self.jumps[a][1][gene.name].append([b,jumps[(a,b)]])
                self.jumps[b][0][gene.name].append([a,jumps[(a,b)]])
            self.info[gene.name] = gene.info
	    self.members.append(gene)  	    
	    if not self.open: break  
	
        return self.prev_chr,self.jumps,self.starts,self.ends,self.chr_jxns,self.info








    def getGeneFromGencodeFile(self):	

	self.chr,feature,start,end,strand,gID = self.line[0],self.line[2],int(self.line[3]),int(self.line[4]),self.line[6],self.line[9].split('\"')[1]
	if self.prev_chr != self.chr: return GtfGene() 		
	if feature != 'gene': gtf_error('GTF Filetype Unfamiliar') 
	gene   = GtfGene(start,end,self.apply_genenick(gID,self.line[15].split('\"')[1]),strand)
	while True: 
		self.line = self.filehandle.readline().split()
		if len(self.line) == 0: 
			self.open = False
			return gene 
		elif self.line[2] == 'exon': 	
			start, end, geneName = int(self.line[3]), int(self.line[4]),  self.apply_genenick(self.line[9].split('\"')[1],self.line[17].split('\"')[1]) 
			gene.add_exon(start,end,self.apply_trannick(geneName,self.line[11].split('\"')[1],self.line[23].split('\"')[1]))
		elif self.line[2] == 'gene': return gene 


	
    def getGeneFromRhesusFile(self):	




	self.chr,feature,start,end,strand,gID = self.line[0],self.line[2],int(self.line[3]),int(self.line[4]),self.line[6],self.line[9].split('\"')[1]


	if self.prev_chr != self.chr: return GtfGene() 		

	if feature != 'aGene': gtf_error('GTF Filetype Unfamiliar') 
	gene   = GtfGene(start,end,self.apply_genenick(gID,self.line[13].split('\"')[1]),strand)
	while True: 
		self.line = self.filehandle.readline().split()
		if len(self.line) == 0: 
			self.open = False
			return gene 
		elif self.line[2] == 'exon': 	
			start, end, geneName = int(self.line[3]), int(self.line[4]),  self.apply_genenick(self.line[9].split('\"')[1],self.line[13].split('\"')[1]) 

			gene.add_exon(start,end,self.apply_trannick(geneName,self.line[9].split('\"')[1],self.line[11].split('\"')[1]))
		elif self.line[2] == 'aGene': return gene 


	

		


		




    def getline(self):
	if self.GTF_STYLE == 'RHESUS':
            try:
                mychr,feature,start,end,strand,gID = self.line[0],self.line[2],int(self.line[3]),int(self.line[4]),self.line[6],self.line[9].split('\"')[1]
                tID,tName,geneType = self.line[11].split('\"')[1],self.line[11].split('\"')[1],self.line[13].split('\"')[1]
                
                geneStatus,hugoName = self.line[15].split('\"')[1],self.line[17].split('\"')[1]
                #return mychr,feature,start,end,strand,gID,tID,geneType,geneStatus,hugoName
            except IndexError:
                self.open = False
                return False,False,False,False,False,False,False
        elif self.GTF_STYLE == 'HUMAN':
	
	    if self.current_gene == None:
		try:	mychr,feature,start,end,strand,gID = self.line[0],self.line[2],int(self.line[3]),int(self.line[4]),self.line[6],self.line[9].split('\"')[1]
		except IndexError: 
			self.open = False 
              		return False,False,False,False,False,False,False
		if feature != 'gene': sage_error('GTF Filetype Unfamiliar') 
		
		return mychr,feature,start,end,strand,self.apply_genenick(gID,self.line[15].split('\"')[1]),None		
	
	    else:

		
                if feature == 'gene':
                    tID,tName,hugoName = self.line[9].split('\"')[1],    self.line[15].split("\"")[1], self.line[15].split('\"')[1]
                else:
                    tID,tName,hugoName = self.line[11].split('\"')[1], self.line[23].split("\"")[1], self.line[17].split('\"')[1]
                    #return mychr,feature,start,end,strand,gID,tID,geneType,geneStatus,hugoName
            #except IndexError:
             #   self.open = False
              #  return False,False,False,False,False,False,False
        else:
		gtf_error('Unknown gtf type') 




        geneName = self.apply_genenick(gID,hugoName)

        if feature == 'gene':
            return mychr,feature,start,end,strand,geneName,geneName
        else:
            return mychr,feature,start,end,strand,geneName,self.apply_trannick(geneName,tID,tName)











































class GtfGene:
    def __init__(self, start=None,stop=None,gID=None,strand="+"):


        if gID!=None:
	    self.empty = False
            self.start,self.stop = int(start),int(stop) 
            self.name = gID
            self.strand = strand
            self.info = [start,stop,gID,strand]

        else:
	    self.empty = True
            self.id,self.name,self.strand = None,None,None


        self.transcripts = []
        self.exons = dd(list)


    def add_exon(self,start,stop,tID):
        self.exons[tID].append([start,stop])





    def update(self,readLen):
        readlen = readLen -1
        starts,stops,valid_jxns,jumps = dd(list),dd(list),dd(list),dd(list)
        for transcript in self.exons.keys():
            for i,(s,e) in enumerate(sorted(self.exons[transcript])):
                starts[s].append([transcript,i+1]); stops[e].append([transcript,i+1])
                if e-s < readlen:
                    dist,S,E = e-s,(s,e),(s,e)
                else:
                    dist,S,E = readlen,(s,s+readlen),(e-readlen,e)
                if i == 0:
                    trailing_jxns = [[E,dist]]

                else:
                    jumps[(trailing_jxns[-1][0][-1],s)].append([transcript,i+1,i+2])
                    remainder = readlen

                    if sum([t[-1] for t in trailing_jxns]) > readlen:
                        for j in range(len(trailing_jxns)-1,-1,-1):
                            if trailing_jxns[j][1] >= remainder:    trailing_jxns[j] = [(trailing_jxns[j][0][1]-remainder,trailing_jxns[j][0][1]),remainder]
                            remainder -= trailing_jxns[j][1]
                            if remainder == 0: break
                        trailing_jxns = trailing_jxns[j::]

                    if sum([t[-1] for t in trailing_jxns]) > readlen:
                        print 'wtf',n
                        sys.exit()
                    else:
                        #print trailing_jxns
                        if len(trailing_jxns) == 1:
                            valid_tuple = tuple([x for x in tuple([t[0] for t in trailing_jxns]+[S])])
                            valid_cnt   = "-".join([str(x) for x in range(i-len(valid_tuple)+2,i+2)])
                            valid_jxns[valid_tuple].append([transcript,valid_cnt])
                            trailing_jxns.append([E,dist])
                        else:

                            needlen = readlen-sum([t[-1] for t in trailing_jxns[1::]])-1
                            if dist >= needlen:
                                valid_tuple = tuple([x for x in tuple([t[0] for t in trailing_jxns]+[(S[0],S[0]+needlen)])])
                                valid_cnt   = "-".join([str(x) for x in range(i-len(valid_tuple)+2,i+2)])
                                valid_jxns[valid_tuple].append([transcript,valid_cnt])
                                trailing_jxns,trailing_len = [t for t in trailing_jxns[1::]],sum([t[1] for t in trailing_jxns[1::]])
                                if dist+trailing_len > readlen:
                                    new_tuple = tuple([t[0] for t in trailing_jxns]+[S])
                                    new_cnt   = "-".join([str(x) for x in range(i-len(new_tuple)+2,i+2)])
                                    valid_jxns[new_tuple].append([transcript,new_cnt])
                            else:
                                valid_tuple = tuple([x for x in tuple([t[0] for t in trailing_jxns]+[S])])
                                valid_cnt   = "-".join([str(x) for x in range(i-len(valid_tuple)+2,i+2)])
                                valid_jxns[valid_tuple].append([transcript,valid_cnt])
                            trailing_jxns.append([E,dist])


        return starts,stops,jumps,valid_jxns





