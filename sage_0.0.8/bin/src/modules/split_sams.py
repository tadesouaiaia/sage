#!/usr/bin/env python
import os
import sys
from collections import defaultdict as dd



def revcomp(seq):
    C = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join([C[x] for x in seq[-1::-1]])







class SamReads:
    def __init__(self,fileHandle,outprefix,fprefix,ftype='CHRS',options=None):


	if options == None: 
		self.readlen, self.subdist = 100, 0 
	else:
		self.readlen, self.subdist = options.readlen, options.subdist

        self.qual = "".join(['D' for i in range(self.readlen)])



	

        self.f = open(fileHandle)
        self.line = self.f.readline().strip()
        self.s    = self.line.split()
        while len(self.s)>0 and self.line[0][0] == "@":
            self.line = self.f.readline().strip()
            self.s = self.line.split()
        if len(self.line) == 0: self.open = False
        else:                   self.open = True

	self.outpath, self.prefix, self.ftype = outprefix , fprefix, ftype

	self.vispath = self.outpath+'/vis'
	self.cntpath = self.outpath+'/cnts'
	self.sampath = self.outpath+'/sams'
	self.visfiles = {}

 
	if self.ftype == 'CHRS': 
		self.next = self.next_chr_map
		self.create_lines = self.create_chr_lines
		self.process = self.process_chr_map
	elif self.ftype == 'NOV':
		self.next = self.next_novel_map
	else:
		self.next = self.next_jxn_map 
		self.create_lines = self.create_jxn_lines
		self.process = self.process_jxn_map
    		self.gene_cnts,self.iso_cnts = dd(int),dd(int)

   
    def close_vis_files(self):

	for k in self.visfiles.keys():

		self.visfiles[k].close() 
    
    def get_pts(self,pStart,jPairs):
            remLen,pts,jumps,pSkip= self.readlen,[],[],False
            for p,(p1,p2) in enumerate(jPairs):
                pGo = p1+pStart
                if pGo > p2:
                            pStart -= (p2-p1)+1
                            pSkip = True
                elif remLen > p2-pGo:
                            pts.append((pGo,p2))
                            remLen -= 1+p2-pGo
                            pStart  = 0 
                else:
                    pts.append((pGo,pGo+remLen-1))
                    break
            cigar = [str(pts[0][1]-pts[0][0])+"M"]
            for i in range(1,len(pts)):
                jumps.append(str(pts[i-1][1])+'-'+str(pts[i][0]))
                cigar.append(str(pts[i][0] - pts[i-1][1])+'N')
                cigar.append(str(pts[i][1]-pts[i][0])+'M')
            return pts[0][0],tuple(pts),jumps,"".join(cigar)




    def determine_locations(self,ref,loc):
        jxns,gene,isos = ref.split("=")[1].split(":")[0].split("-"),ref.split(":")[2].split("=")[0],ref.split("=")[-1]

        gene_data = [g.split("=") for g in ":".join(ref.split(":")[2::]).split("/")]

        jJumps = [jxns[i]+"-"+jxns[i+1] for i in range(1,len(jxns)-2,2)] 
        jxns = [(int(jxns[i-1]),int(jxns[i])) for i in range(1,len(jxns),2)] 
        cLoc,cPts,cJumps,cCigar = self.get_pts(loc,jxns)
        if cJumps != jJumps:
            foundJ,pass_data=[i for i in range(len(jJumps)) if jJumps[i]  in cJumps],[]
            for gene,isos in gene_data:
                iso,vIsos = isos.split("|"),[]
                for ob in iso:
                    iName,iJumps = ob.split(":")[0],ob.split(":")[1].split("-")
                    vIsos.append(iName+":"+"-".join([iJumps[j] for j in foundJ]))
                    if len(foundJ) > 0:
                        vIsos.append(iName+":"+"-".join([iJumps[j] for j in foundJ]+[iJumps[foundJ[-1]+1]]))
                    else:
                        ### ME ###
                        vIsos.append(iName+":"+"-".join([iJumps[j] for j in range(len(iJumps))]))
                pass_data.append([gene,"|".join(vIsos)])                    

            gene_data = pass_data

        return cLoc,cPts,tuple(cJumps),cCigar,gene_data






    def record_expression(self): 
	gene_file =  open(self.cntpath+'/'+self.prefix+'.genes','w')
	iso_file =  open(self.cntpath+'/'+self.prefix+'.isos','w')
	for g in self.gene_cnts.keys(): gene_file.write('%s %s\n' % (g,self.gene_cnts[g]))
	for g in self.iso_cnts.keys(): iso_file.write('%s %s\n' % (g,self.iso_cnts[g]))




    def	next_novel_map(self):
        self.fail,self.unique = False, True
        self.id,self.seq,self.ref_chr,self.ref,self.loc,self.line,self.data,self.ref_data = self.s[0],self.s[9],  self.s[2],self.s[2],int(self.s[3]),self.line,[],[]
        if self.ref_chr not in self.visfiles.keys(): self.visfiles[self.ref_chr] = open(self.vispath+'/'+self.prefix+'.'+self.ref_chr+'.vis','w')
	self.visfiles[self.ref_chr].write(self.line) 
        self.line = self.f.readline().strip()
        self.s = self.line.split()
        if len(self.s) == 0:
        	self.open = False




    def next_chr_map(self):
        self.fail,self.unique, self.data, self.ref_data = False, True,[], [] 
        self.id,self.seq,self.ref_chr,self.ref,self.loc,self.subs,self.line = self.s[0],self.s[9],  self.s[2],self.s[2],int(self.s[3]),int(self.s[-1].split(":")[-1]),self.line
        if self.s[1] == '0':self.rev_seq = revcomp(self.seq)
        else:               self.rev_seq,self.seq = self.seq,revcomp(self.seq)

        while self.open and self.id == self.s[0]:
            ref_chr,ref,loc,subs,cigar,strand =  self.s[2].split(":")[0],self.s[2],int(self.s[3]),int(self.s[-1].split(":")[-1]),self.s[5],self.s[1]
            self.data.append((subs,loc,cigar,strand,ref_chr))
            self.line = self.f.readline().strip()
            self.s = self.line.split()
            if len(self.s) == 0:	self.open = False

        if len(self.data) == 1:	self.chrs = [self.data[0][-1]]
        elif len(self.data)>1:
            self.data.sort()
            subs,subs_set = sorted([d[0] for d in self.data]),sorted(list(set([d[0] for d in self.data])))
            if len(subs_set)>1 and subs_set[0] + 1 < subs_set[1]:	
		self.data = [d for d in self.data if d[0] == min(subs)]
            self.chrs = list(set([d[-1] for d in self.data]))
            if len(list(set([d[-1] for d in self.data]))) > 1: self.fail = True
        if not self.fail:
            self.create_lines()
        return self




    def create_chr_lines(self):
        self.sam_lines = []
        for d in self.data:
            subs,loc,cigar,strand,chr = d
            NM,SN,SR,OR='NM:i:'+str(len(self.data)),'SN:i:'+str(subs),'SR:i:1','OR:Z:'+str(strand)
            self.sam_lines.append([self.id,strand,chr,str(loc),'255',cigar,'*','0','0',self.seq,self.qual,NM,SN,SR,OR])




    def process_chr_map(self): 
	if not self.fail and len(self.chrs) == 1:
        	if self.chrs[0] not in self.visfiles.keys(): self.visfiles[self.chrs[0]] = open(self.vispath+'/'+self.prefix+'.'+self.chrs[0]+'.vis','w')
                for sam_line in self.sam_lines:  self.visfiles[self.chrs[0]].write("\t".join(sam_line)+'\n')






















    def next_jxn_map(self):
        self.fail, self.data, self.ref_data = False, [], [] 
        self.id,self.seq,self.ref_chr,self.ref,self.loc,self.subs,self.line = self.s[0],self.s[9],self.s[2].split(":")[0],self.s[2],int(self.s[3]),int(self.s[-1].split(":")[-1]),self.line

	
        if self.s[1] == '0':	self.rev_seq = revcomp(self.seq)
        else:			self.rev_seq,self.seq = self.seq,revcomp(self.seq)

        while self.id == self.s[0]:
            ref_chr,ref,loc,subs =  self.s[2].split(":")[0],self.s[2],int(self.s[3]),int(self.s[-1].split(":")[-1])
            cLoc,cPts,cJumps,cCigar,gene_data = self.determine_locations(ref,loc)
            self.data.append((subs,(cLoc,cPts,cJumps,cCigar,self.s[1],ref_chr),gene_data))
            self.line = self.f.readline().strip()
            self.s = self.line.split()
            if len(self.s) == 0:
                self.open = False
                break
        if len(self.data)>1:
            self.data.sort()
            self.sub_dict,self.unique = dd(lambda: dd(list)),False
            for d in self.data: self.sub_dict[d[0]][d[1]].extend(d[2])

	    sub_keys = sorted(self.sub_dict.keys())
            if len(sub_keys) > 1: 
		if sub_keys[0] + self.subdist <= sub_keys[1]: 
			self.sub_dict = {sub_keys[0]: self.sub_dict[sub_keys[0]]} 

	

            if len(self.sub_dict.keys()) == 1:
                loc_dict,gene_dict = self.sub_dict[self.sub_dict.keys()[0]],dd(list)
                if len(loc_dict) == 1:
                    for v,t in loc_dict.values()[0]: gene_dict[v].append(t)
                    self.loc_genes = [[g,"|".join(t)] for g,t in gene_dict.items()]
                    self.data = [(self.sub_dict.keys()[0],loc_dict.keys()[0],self.loc_genes)]
                else:
                    for V in loc_dict.values():
                        for v,t in V:   gene_dict[v].append(t)
                    self.loc_genes = [[g,"|".join(t)] for g,t in gene_dict.items()]
            else:
		self.fail, gene_dict  = True, dd(list) 
		for loc_dict in self.sub_dict.values(): 
			for V in loc_dict.values(): 
				for v,t in V: gene_dict[v].append(t) 
                self.loc_genes = [[g,"|".join(t)] for g,t in gene_dict.items()]
        if not self.fail:
            self.create_lines()
        return self




        
    def create_jxn_lines(self):
        if len(self.data) == 1:
            self.unique = True
            self.subs = self.data[0][0]
            self.loc,self.pairs,self.jxns,self.cigar,self.strand,self.chr = self.data[0][1]
            self.genes = [g[0] for g in self.data[0][-1]]
            self.isos = [g[0]+"="+g[1] for g in self.data[0][-1]]
            self.chrs = [self.chr]
        else:
            self.unique=False
            self.subs = self.data[0][0]
            self.genes = sorted(list(set([lg[0] for lg in self.loc_genes])))
            self.chrs = sorted(list(set([d[1][-1] for d in self.data])))
            self.isos = [g[0]+"="+g[1] for g in self.loc_genes]

        self.sam_lines,self.v_lines = [],[]
        for i,d in enumerate(self.data):
            subs = d[0]
            loc,pairs,jxns,cigar,strand,chr =d[1]
            if strand == '0': seq = self.seq
            else:             seq = self.rev_seq
            Sgenes,Sisos="/".join([x[0] for x in d[2]]),"/".join([x[0]+"="+x[1] for x in d[2]])
            NM,SN,SR,OR='NM:i:'+str(len(self.data)),'SN:i:'+str(subs),'SR:i:1','OR:Z:'+str(strand)
            GE,TR="GE:Z:"+"/".join([x[0] for x in d[2]]),'TR:Z:'+"/".join([x[0]+"="+x[1] for x in d[2]])
            self.sam_lines.append([self.id,strand,chr,str(loc),'255',cigar,'*','0','0',seq,self.qual,NM,SN,SR,OR,"MT:Z:KJXN",GE,TR,'FT:Z:EXONIC','AD:Z:NA'])
            vLine = ["|".join([str(p[0])+'-'+str(p[0]) for p in pairs]),Sisos,'EXONIC','KJXN']
            self.v_lines.append(vLine+[vLine[0]]+[NM])





















    def process_jxn_map(self): 	
        #if len(self.chrs) == 1:
        if not self.fail and len(self.chrs) == 1:
        	if self.chrs[0] not in self.visfiles.keys(): 
			self.visfiles[self.chrs[0]] = open(self.sampath+'/'+self.prefix+'.'+self.chrs[0]+'.sam','w')
                	#vert_files[samread.chrs[0]] = open(prefix+".kjxns."+samread.chrs[0]+".vert",'w')



                for sam_line in self.sam_lines:  self.visfiles[self.chrs[0]].write("\t".join(sam_line)+'\n')



		if self.unique:
			sam_gene = self.chrs[0]+":KJXN:"+"/".join(self.genes)
			sam_isos = self.chrs[0]+":KJXN:"+"/".join(self.isos)
			self.gene_cnts[sam_gene]+=1
			self.iso_cnts[sam_isos]+=1

		else:
			if len(self.genes)<5:    self.gene_cnts[self.chrs[0]+":KJXN:"+"/".join(self.genes)]+=1
			if len(self.isos)<4:     self.iso_cnts[self.chrs[0]+":KJXN:"+"/".join(self.isos)]+=1






























































if __name__ == '__main__':

    from optparse import OptionParser
    usage = "usage: ./%prog [options] data_file"
    parser = OptionParser(usage=usage)


    parser.add_option("-r", "--readlen", default = 100, type=int, help="Output Filename Prefix")
    parser.add_option("-p", "--prefix", default = None, type='string', help="Output Filename Prefix")

    (options, args) = parser.parse_args()


    if not options.prefix:
        print 'need prefix'
        sys.exit()

    split_mapping(args[0],options.readlen,options.prefix)























