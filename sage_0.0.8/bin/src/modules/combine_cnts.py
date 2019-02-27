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

	self.prefix,self.file_key,self.swap_key = 'sage',{},{}

	
	self.swap_gene_type = {'EXONIC': ['CODING','TOTAL'],'KJXN': ['CODING','TOTAL'], 'INTRONIC': ['INTRONIC','TOTAL'],'INTERGENIC': ['INTERGENIC']}
	for k in ['coding','total','intronic','intergenic','feature']:	self.file_key[k.upper()] = open(self.outpath+self.prefix+'.'+k+'.cnts','w')


	self.swap_key['chrX~NA'] = 'chrX~000~NA'

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
		elif len(name.split('|')) > 1:	
			try: return "|".join([self.swap_key[cr+'~'+n] for n in name.split('|')])
			except KeyError: return cr+'~999~'+name

		else:				return cr+'~999~'+name 





	
	


    def translate_type(self,ref,TYPE='GENE'): 

	
	rs = ref.split(':') 


	cr,ct,names = ref.split(':')[0],ref.split(':')[1],":".join(ref.split(':')[2::]).split('/')
	
	if TYPE == 'GENE':
		swapped_types =   self.swap_gene_type[ct]
		
		if ct != 'INTERGENIC': return swapped_types,"/".join([self.swap_name(cr,gene_name) for gene_name in names])
		else:		       return swapped_types,"/".join([">".join([self.swap_name(cr,n) for n in N.split('>')]) for N in names])

	else:
		try:	
			if ct in ['KJXN','EXONIC','INTRONIC']:	
				return "/".join([self.swap_name(cr,nd.split('=')[0])+':'+ct+'='+"=".join(nd.split('=')[1::]) for nd in names])

			elif ct in ['NOVEL']:
				return cr+':NJXN:'+names[0].split(':')[0]+'='+"|".join([self.swap_name(cr.split('>')[-1],nd) for nd in names[0].split(':')[-1].split('|')])	
			else:
				return None
		except:
			sys.stderr.write('warning: translate error\n')
			return None



    def add_gene_cnts(self,files):

	outputs = dd(lambda: dd(int)) 
	for f in files:
		for line in open(f):
			ref,cnt = line.split()
			keys,name = self.translate_type(ref)
			for k in keys:	outputs[k][name] += int(cnt) 



	for fk in outputs.keys():
		for name,cnt in outputs[fk].items():
			self.file_key[fk].write('%s %s\n' % (name,cnt))

	 

    def add_iso_cnts(self,files):

 
	iso_output = dd(int) 
	for f in files:
		for line in open(f):
			ref,cnt = line.split()
			f_name = self.translate_type(ref,'FEATURE') 
			if f_name != None:
				iso_output[f_name] +=1 	


	for name,cnt in iso_output.items():
		self.file_key['FEATURE'].write('%s %s\n' % (name,cnt))
