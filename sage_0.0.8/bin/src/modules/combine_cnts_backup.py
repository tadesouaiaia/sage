#!/usr/bin/env python
import os
import sys
from collections import defaultdict as dd



def revcomp(seq):
    C = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join([C[x] for x in seq[-1::-1]])







class CombineCnts:
    def __init__(self,outpath,feature_key):

	self.feature_key = feature_key
	self.outpath = outpath+'/'


	
	self.swap_key = {}

	self.swap_gene_type = {'EXONIC': ['CODING','TOTAL'],'KJXN': ['CODING','TOTAL'], 'INTRONIC': ['INTRONIC','TOTAL'],'INTERGENIC': ['INTERGENIC']}


	for c,f in feature_key.items():
		cr = f.split('.feature.key')[0].split('.')[-1]
		self.swap_key[cr+'~NA'] = cr+'~000~NA'
		for line in open(f):
			line = line.split() 
			if line[1] == 'GENE':
				gL = [line[0],line[3].split('@')[0].split('.')[0], line[3].split('@')[-1]]	
				self.swap_key[line[0]+'~'+line[2]] = "~".join(gL)
		
    def swap_name(self,cr,name):

	try:
		return self.swap_key[cr+'~'+name] 

	except KeyError:

		if 'copy' in name: return cr+'~999~'+name
		print name		
		print 'ERROR',cr,name
		sys.exit() 




	
	


    def translate_type(self,ref,TYPE='GENE'): 

	
	rs = ref.split(':') 
	cr,ct,names = ref.split(':')[0],ref.split(':')[1],":".join(ref.split(':')[2::]).split('/')

	swapped_types =   self.swap_gene_type[ct]
	
	
	if TYPE == 'GENE':
		
		if ct != 'INTERGENIC': return swapped_types,"/".join([self.swap_name(cr,gene_name) for gene_name in names])
		else:		       return [],None 


	


	else:
		print TYPE
		sys.exit()
	sys.exit()  
	gUpdate = [] 

	cr,ct = name.split(':')[0],name.split(':')[1]
	gn = ":".join(name.split(':')[2::]).split('/')
	gUpdate = [] 
	for gx in gn:

		if ct in ['EXONIC','KJXN','INTRONIC']:			
			if cr+'~'+gx in self.swap_key:
				gUpdate.append(self.swap_key[cr+'~'+gx])
			else:
				gUpdate.append(cr+'~000~'+gx)
		elif ct in ['INTERGENIC']:
			gout = '>'.join([self.swap_key[cr+'~'+gi] if cr+'~'+gi in self.swap_key else cr+'~000~'+gi for gi in gx.split('>')])
		else:
			print cr,ct,gn		

	return
	if TYPE == 'GENE':

		cr,ct,gn = name.split(':') 
#		print cr,ct,gn	


		gSlash = gn.split('/') 
		gCarrot = gn.split('>') 

		if len(gSlash) > 1:		
			for gx in gSlash:
				if cr+'~'+gx in self.swap_key: gUpdate.append(self.swap_key[cr+'~'+gx])
				else:			       gUpdate.append(cr+'~000~'+gx) 
			return '/'.join(gUpdate) 


		elif len(gCarrot) > 1:
			for gx in gCarrot:
				if cr+'~'+gx in self.swap_key: gUpdate.append(self.swap_key[cr+'~'+gx])
				else:			       gUpdate.append(cr+'~000~'+gx) 
			return '>'.join(gUpdate) 

		else:
			if cr+'~'+gn in self.swap_key:	return self.swap_key[cr+'~'+gn] 
			else:			        return cr+'~000~'+gn 
			




    def add_gene_cnts(self,files):
	ref_cnts = dd(int)

	self.prefix = 'sage'
#        self.coding = open(self.outpath+self.prefix+'.coding.cnts','w')
#        self.total = open(self.outpath+self.prefix+'.total.cnts','w')
#        self.intergenic = open(self.outpath+self.prefix+'.intergenic.cnts','w')
#        self.intronic = open(self.outpath+self.prefix+'.intronic.cnts','w')

	self.file_key = {} 
	self.file_key['CODING'] =  open(self.outpath+self.prefix+'.coding.cnts','w')
        self.file_key['TOTAL']  = open(self.outpath+self.prefix+'.total.cnts','w')
        self.file_key['INTRONIC']  = open(self.outpath+self.prefix+'.intronic.cnts','w')
        self.file_key['INTERGENIC']  = open(self.outpath+self.prefix+'.intronic.cnts','w')
        self.file_key['ISOS']  = open(self.outpath+self.prefix+'.feature.cnts','w')


	self.swap_gene_type = {'EXONIC': ['CODING','TOTAL'],'KJXN': ['CODING','TOTAL'], 'INTRONIC': ['INTRONIC','TOTAL'],'INTERGENIC': ['INTERGENIC']}
	outputs = dd(lambda: dd(int)) 
 
	for f in files:


		for line in open(f):
			ref,cnt = line.split()
			rs = ref.split(':') 

			keys,name = self.translate_type(ref)
			for k in keys:
				outputs[k][name] += int(cnt) 
			continue



			if rs[1] in ['EXONIC','KJXN']:
				ref_cnts[rs[0]+':CODING:'+rs[2]] += int(cnt) 
				ref_cnts[rs[0]+':TOTAL:'+rs[2]] += int(cnt) 
			elif rs[1] in ['INTRONIC']:
				ref_cnts[rs[0]+':INTRONIC:'+rs[2]] += int(cnt) 
				ref_cnts[rs[0]+':TOTAL:'+rs[2]] += int(cnt) 
			elif rs[1] in ['INTERGENIC']:
				ref_cnts[rs[0]+':INTERGENIC:'+rs[2]] += int(cnt) 
			else:
				 print rs 


	for fk in outputs.keys():
		for name,cnt in outputs[fk]:
			self.file_key[fk].write('%s %s\n' % (name,cnt))

	 

	for r,c in ref_cnts.items():

		r_chr = r.split(':')[0] 
		r_type = r.split(':')[1] 
		r_genes = r.split(':')[2].split('/') 
		r_regs = r.split(':')[2].split('>') 






		if r_type == 'CODING':
			try: r_swap = "/".join([self.swap_key[r_chr+'~'+rg] for rg in r_genes])
			except KeyError:	r_swap = "/".join([r_chr+'~'+rg for rg in r_genes])
			self.coding.write('%s %s\n' % (r_swap,c))
		elif r_type == 'TOTAL':
			try: r_swap = "/".join([self.swap_key[r_chr+'~'+rg] for rg in r_genes])
			except KeyError:	r_swap = "/".join([r_chr+'~'+rg for rg in r_genes])
			self.total.write('%s %s\n' % (r_swap,c))
		elif r_type == 'INTRONIC': 
			try:  r_swap = "/".join([self.swap_key[r_chr+'~'+rg] for rg in r_genes])
			except KeyError: 
				try:                    r_swap = "/".join([self.swap_key[r_chr+'~'+rg.split('@')[0]] for rg in r_genes])
				except KeyError:	r_swap = "/".join([r_chr+'~'+rg for rg in r_regs])
			self.intronic.write('%s %s\n' % (r_swap,c))
		elif r_type == 'INTERGENIC': 	

			try: 			r_swap = "/".join([self.swap_key[r_chr+'~'+rg] for rg in r_regs])
			except KeyError:	r_swap = "/".join([r_chr+'~'+rg for rg in r_regs])
				
			self.intergenic.write('%s %s\n' % (r_swap,c))
		else:
			print r_type,r_swap



    def add_iso_cnts(self,files):
	ref_cnts = dd(int)

	self.prefix = 'sage'
        self.isos = open(self.outpath+self.prefix+'.feature.cnts','w')
 
	for f in files:
		for line in open(f):
			ref,cnt = line.split()
			self.translate_name(ref) 
			rs = ref.split(':') 
			if rs[1] in ['INTERGENIC']: continue 
			
			elif rs[1] == 'KJXN': 
				if len(ref.split('=')) > 2: continue 
				try: r_swap = self.swap_key[rs[0]+"~"+rs[2].split('=')[0]]+'~KJXN='+ref.split('=')[1]
				except KeyError: r_swap = ref 
				ref_cnts[r_swap] += int(cnt) 
			elif rs[1] == 'EXONIC': 
				if len(ref.split('=')) > 2: continue 
				try:    r_swap = self.swap_key[rs[0]+"~"+rs[2].split('=')[0]]+'~EXON='+ref.split('=')[1]
				except KeyError: r_swap = ref 
				ref_cnts[r_swap] += int(cnt) 
			elif rs[1] in ['INTRONIC']:
				g1,g2 = rs[2].split('=')[0], rs[2].split('>')[-1].split('=')[0] 
				if g1 != g2: continue 
				try: 
					ref_cnts[self.swap_key[rs[0]+'~'+g1]+'~INTRONIC='+rs[2].split('>')[0]+'>'+rs[2].split('>')[-1].split('=')[1]] += int(cnt) 
				except KeyError:
					ref_cnts[ref] += int(cnt) 
			elif rs[1] == 'NOVEL':
				try:             
					rc = rs[0].split('>')[-1] 
					rs_fin = [self.swap_key[rc+'~'+x] if x != 'NA' else 'NA' for x in rs[-1].split('|')]
					ref_cnts[rc+'~NJXN~'+rs[2]+':'+rs[3]+':'+rs[4]+':'+"|".join(rs_fin)] += int(cnt)
				except KeyError:
					ref_cnts[ref] += int(cnt) 

			else:
				print 'hmmm'


	for r,c in ref_cnts.items():
		self.isos.write('%s %s\n' % (r,c))
