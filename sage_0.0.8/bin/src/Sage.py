#!/usr/bin/env python2.7

import sys
import os
from collections import defaultdict as dd
from subprocess import call
import subprocess
import multiprocessing


from modules.GtfAnnotate import GtfFile 
from modules.ReadFiles   import ReadFile			
from modules.track_progress import dot_log
from modules.parse_vis import Annotation, NovelAnnotation
from modules.split_sams import SamReads
from modules.combine_cnts import CombineCnts

from contextlib import contextmanager

@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)




def sage_error(eString):
    sys.stderr.write('\nSageError: '+eString+'\n')
    sys.exit()









def read_config_file(config):
	config_key = dd(bool) 
      	parser.add_argument('--seqs',type=str,default=None,help='Key File (Reqired)')
	if type(config) == str:
		for line in open(config):

			line = line.strip().split('=')
			config_key[line[0]] = line[1]
		return config_key
	for line in config:
		line = line.strip().split('=')
		try: config_key[line[0]] = line[1]
		except IndexError: continue

	config.close()
	return config_key 
		

def add_dirs(parent,children,grandchildren=[]):

	if not os.path.exists(parent):
        	os.makedirs(parent)

	for c in children:
		if not os.path.exists(parent+'/'+c):
        		os.makedirs(parent+'/'+c)
		for g in grandchildren:
			if not os.path.exists(parent+'/'+c+'/'+g):
        			os.makedirs(parent+'/'+c+'/'+g)






class Sage:
        def __init__(self,args,sagedir,rundir):


		self.progress = dot_log(args) 	
		self.args = args
		self.sagedir = sagedir
		self.rundir = rundir
		self.pp, self.cp = self.sagedir +'/src/compiled/perm', self.sagedir +'/src/compiled/clipR'

		if args.local: 
			self.pp, self.cp = self.sagedir +'/src/compiled/perm_local', self.sagedir +'/src/compiled/clipR' 
		else:
			self.pp, self.cp = self.sagedir +'/src/compiled/perm_server', self.sagedir +'/src/compiled/clipR'

		if args.cores == 0:
			args.cores = int(multiprocessing.cpu_count())-1

		if args.option == 'annotate':

 
			if args.outpath == None: args.outpath = self.rundir+'/sage_anno'
			if args.prefix  == None: args.prefix   = 'sage' 

			args.chrs = os.path.abspath(args.chrs)
			args.outpath = os.path.abspath(args.outpath) 

			add_dirs(args.outpath,['anno','lists','refs','indexes','filters','stats'])



			self.progress.start_minor('Loading GTF File') 
			anno = GtfFile(args.gtf, args.chrs,args) 


			self.progress.start_minor('Indexing GTF File',1) 

   	 		while anno.open:
        			anno.printData()
				self.progress.mark() 
    			anno.print_unannotated_chrs(args.outpath+'/filters')



			cList = anno.listPrefix+'/CHRS.txt'
			jList = anno.listPrefix+'/JXNS.txt'


			for i in [35,70,100]:
				if args.print_gene_seqs: continue 
				try:		
					if i == 35: seeds = ['F1']
					else:       seeds = ['F3','F4'] 
					for seed in seeds:
						seed_map = {'F0': '0','F4': '8', 'F2': '4', 'F1': '2', 'F3': '6'}			
						log = args.outpath+'/build.log'			
						sn = args.outpath+'/indexes/'+'CHRS_'+str(i)+"_"+seed_map[seed]+'.index'
						call([self.pp,cList,str(i),'--readFormat','fastq','--seed',seed,'-s',sn],stdout=open(log,'w'))
						sn = args.outpath+'/indexes/'+'JXNS_'+str(i)+"_"+seed_map[seed]+'.index'
						call([self.pp,jList,str(i),'--readFormat','fastq','--seed',seed,'-s',sn],stdout=open(log,'w'))
				except OSError:
					sage_error('Version Incompatibility') 


			
			self.progress.start_minor('Saving Gtf Configuration')

			w=open(args.outpath+'/'+args.prefix+'.config','w') 
			w.write('%s' % 'CHRS='+anno.listPrefix+'/CHRS.txt\n')
			w.write('%s' % 'JXNS='+anno.listPrefix+'/JXNS.txt\n')
			w.write('%s' %  'ANNOS='+args.outpath+'/anno'+'\n')
			w.write('%s' %  'IDXS='+args.outpath+'/indexes'+'\n')
			self.progress.end() 


		elif args.option == 'qc': 
			if args.outpath == None: args.outpath = self.rundir+'/sage_reads'
			self.qc(args,args.outpath,args.prefix)

		elif args.option == 'align':
			self.progress.start_major('Alignment') 
			if args.outpath == None: args.outpath = self.rundir+'/sage_map'
			add_dirs(args.outpath,['reads','trims','maps'])
			self.link_read_config(args.readconfig.name,args.outpath+'/reads/')
			self.load_ref_indexes(args.refconfig.name) 
			self.align(args,args.outpath,args.prefix)

		elif args.option == 'quantify':

			self.load_ref_config(args.refconfig.name)
			if args.outpath == None: args.outpath = self.rundir+'/sage_output'
			mapdir = args.mapdir 
			self.parse(mapdir,args.outpath,args.prefix)


		elif args.option == 'easyrun':

			if not args.outpath: args.outpath = self.rundir+'/easyrun'	
			qcdir = args.outpath+'/qc'
			mapdir = args.outpath+'/align'
			parsedir = args.outpath+'/output'

			self.progress.start_major('Easyrun') 
			self.progress.start_minor('Read QC') 	
			
			self.qc(args,qcdir,'sage') 
			add_dirs(mapdir,['reads','trims','maps'])
			self.progress.start_major('Alignment') 
			self.link_read_config(qcdir+'/sage.reads',mapdir+'/reads/')
			self.load_ref_indexes(args.refconfig.name) 
			self.align(args,mapdir,'sage')
		
			self.parse(mapdir+'/maps/',parsedir,'sage')




	def parse(self,mapdir,outpath,outprefix): 


		self.progress.start_major('Quantifying Expression') 
		add_dirs(outpath,['vis','cnts','sams','verts','cands','output'])
		vis_novels,vis_jxns,vis_knowns = dd(list), dd(list) ,dd(list) 

		self.progress.start_minor('Quantifying Gene Expression') 
	
		TEST=False
		for filename in os.listdir(mapdir):

			if TEST: break


			if filename.endswith(".sam"):
				fprefix = filename.split('.sam')[0] 
				if fprefix.endswith('_miss_CHRS'):
					samReads = SamReads(mapdir+'/'+filename,outpath,fprefix,'JXNS', self.args) 
					while samReads.open:	
						samread = samReads.next().process()
						self.progress.mark() 
					samReads.record_expression() 	


				elif filename.split('.sam')[0].endswith('_miss_JXNS'):
					samReads = SamReads(mapdir+'/'+filename,outpath,fprefix,'NOV',  self.args) 
					while samReads.open:	
						self.progress.mark() 
						samread = samReads.next()
				else:



					samReads = SamReads(mapdir+'/'+filename,outpath,fprefix,'CHRS', self.args) 
					while samReads.open:	
						self.progress.mark() 
						samread = samReads.next().process()

			
				samReads.close_vis_files() 



					
		for filename in os.listdir(outpath+'/vis'):
			fprefix, fchr = filename.split('.vis')[0], filename.split('.vis')[0].split('.')[-1]
			


			if ".".join(fprefix.split('.'+fchr)[0:-1]).endswith('miss_CHRS'): 	continue 
			elif ".".join(fprefix.split('.'+fchr)[0:-1]).endswith('miss_JXNS'): 	vis_novels[fchr].append(outpath+'/vis/'+filename)
			else:									vis_knowns[fchr].append(outpath+'/vis/'+filename)


		self.progress.start_minor('Quantifying Isoform Expression') 



		chr_key = {chrfile.split('/')[-1].split('.')[0]  : chrfile.strip() for chrfile in open(self.ref_config['CHRS'])}
		anno_key = {x.split('.')[-2]: self.ref_config['ANNOS']+'/'+x for x in [X for X in os.listdir(self.ref_config['ANNOS']) if X.endswith('annotation')]}
		feature_key = {x.split('.')[-3]: self.ref_config['ANNOS']+'/'+x for x in [X for X in os.listdir(self.ref_config['ANNOS']) if X.endswith('feature.key')]}


		
#		TEST=True



		for k,cr in chr_key.items():
			self.progress.mark(1)

			if TEST:	break
			

			anno = Annotation(k,anno_key[k],outpath).process_mappings(vis_knowns[k]) 
			novel_anno = NovelAnnotation(k,anno_key[k],outpath, readlen = self.args.readlen).add_ref(chr_key[k]).process_mappings(vis_novels[k]) 
		
		self.progress.start_minor('Combining Data Across Genome') 

		
		combos = CombineCnts(outpath+'/output',feature_key)
		combos.add_gene_cnts([outpath+'/cnts/'+x for x in os.listdir(outpath+'/cnts') if x.endswith('genes')])
		combos.add_iso_cnts([outpath+'/cnts/'+x for x in os.listdir(outpath+'/cnts') if x.endswith('isos')])

		self.progress.end() 



































	def qc(self,args,outpath,outprefix): 
		self.progress.start_major('Read QC') 
		add_dirs(outpath,['reads'])
		rpath = outpath+'/reads/reads'

		self.progress.start_minor('Reading Data') 
		self.progress.start_minor('Calculating Quality')
		reads = ReadFile(args.reads,args,outpath,outprefix,rpath) 
		reads.addTrimLengths(args.trimLengths)
		reads.addAdaptors(args.threePrimeAdaptors,args.fivePrimeAdaptors,3)
		reads.loadEntropyDist()
		self.progress.start_minor('Trimming Adaptors')

		reads.filter() 
		rpath = outpath+'/reads/reads'
#		rpath = args.outpath+'/reads/'+args.prefix
#		fw=open(args.outpath+'/'+args.prefix+'.reads','w') 
		fw=open(outpath+'/'+outprefix+'.reads','w') 
		for i in args.trimLengths:
#			fstr = args.outpath+'/reads/reads_'+str(i)+'.txt'
			w=open(outpath+'/reads/'+outprefix+'_'+str(i)+'.txt','w')

			 
			for j in range(1,reads.writeIdxs[i]+1):	w.write(rpath+'_'+str(j)+'_'+str(i)+'.'+reads.readstream.ftype+'\n')
			fw.write(str(i)+'='+outpath+'/reads/'+outprefix+'_'+str(i)+'.txt\n')

		
		self.progress.end()






	def run_perm(self,mpath,rpath,RF,r,rLen,log,K={}):

		ref = self.ref_config[RF]
		if self.index_key[(RF,rLen,8)]:	ref = self.index_key[(RF,rLen,8)]
		elif self.index_key[(RF,rLen,6)]:	ref = self.index_key[(RF,rLen,6)]

		if RF == 'CHRS':	self.progress.start_minor('Aligning to Contiguous Genome...')
		else:			self.progress.start_minor('Aligning to Known Splice Junctions')

	
		with cd(mpath):
			#print 'yo'	
			call([self.pp,self.ref_config[RF],r,'-u','-A','--outputFormat','sam','--printNM','-v','6','--seed','F3','-T',str(rLen),'--noSamHeader'],stdout=open(log,'w'))

		#sys.exit() 

		if rLen == self.read_len:
			if RF == 'CHRS':
				call('ls '+rpath+'*'+str(self.read_len)+'_miss_CHRS.fastq', shell=True,stdout=open(rpath+'JXN.txt','w'))
				return rpath+'JXN.txt'
			else:
				call('ls '+rpath+'*miss_CHRS_miss_JXNS.fastq', shell=True,stdout=open(rpath+'NOV.txt','w'))
				return rpath+'NOV.txt'
		elif RF == 'CHRS':
			call('ls '+rpath+'*_miss_CHRS.fastq', shell=True,stdout=open(rpath+'T_JXN.txt','w'))
			return rpath+'T_JXN.txt'

		else:
			call('ls '+rpath+'*miss_CHRS_miss_JXNS.fastq', shell=True,stdout=open(rpath+'UNMAPPED.txt','w'))
			return rpath+'UNMAPPED.txt'



	def run_clip(self,mpath,rpath,tpath,r,rLen,log,k={}):


		if self.args.local:
			self.progress.start_minor('Skipping Splice Site Detection For Local Version')
 
			tw = open(tpath+'T_CHRS.txt','w')
			for i,f in enumerate(open(r)):
				w = open(tpath+'trim_'+str(i+1)+'.fastq','w')
				for line in open(f.strip()):
					line = line.strip()[0:self.trim_len]
					w.write('%s\n' % line) 	
				w.close()
				tw.write(tpath+'trim_'+str(i+1)+'.fastq\n')
			for j,f in enumerate(open(self.read_config[self.trim_len])):
				w = open(tpath+'trim_'+str(i+j+2)+'.fastq','w')
				for line in open(f.strip()):
					line = line.strip()[0:self.trim_len]
					w.write('%s\n' % line) 	
				w.close()
				tw.write(tpath+'trim_'+str(i+1)+'.fastq\n')
			tw.close()	
			return tpath+'T_CHRS.txt'

		self.progress.start_minor('Novel Splice Junctions')
		ref = self.ref_config['CHRS']
		if self.index_key[('CHRS',35,2)]: ref = self.index_key[('CHRS',35,2)]
		if rLen == self.read_len:
			with cd(mpath):
				call([self.cp,ref,r,'--seed','F1','--anchorL','35','-e','-v','1','-s','-u','--noSamHeader','--ignoreDummyR','40','--ignoreRepeatR','10'],stdout=open(log,'w'))
			call('ls '+rpath+'*miss_CHRS_miss_JXNS_miss_CHRS.fastq | cat - '+self.read_config[self.trim_len], shell=True,stdout=open(rpath+'TRIM.txt','w'))	
			tw = open(tpath+'T_CHRS.txt','w')
			for i,f in enumerate(open(rpath+'TRIM.txt')):
				w = open(tpath+'trim_'+str(i+1)+'.fastq','w')
				for line in open(f.strip()):
					line = line.strip()[0:self.trim_len]
					w.write('%s\n' % line) 	
				w.close()
				tw.write(tpath+'trim_'+str(i+1)+'.fastq\n')
			tw.close()	
			return tpath+'T_CHRS.txt'

		sage_error('CLIP NOT SUPPORTED FOR READLENGTH') 

			




	def align(self,args,outpath,outprefix):


	#	add_dirs(outpath,['reads','trims','maps'])
		rpath, mpath, tpath = outpath+'/reads/', outpath+'/maps/', outpath+'/trims/'



		miss_chrs = self.run_perm(mpath,rpath,'CHRS',self.read_config[self.read_len],self.read_len,mpath+'F_CHR.log')
		miss_jxns = self.run_perm(mpath,rpath,'JXNS',miss_chrs,self.read_len,mpath+'F_JXN.log')
		miss_novs = self.run_clip(mpath,rpath,tpath,miss_jxns,self.read_len,mpath+'F_NOV.log')	

		self.progress.start_major('Trimmed Re-Alignment') 	

		

		miss_trim_chrs = self.run_perm(mpath,tpath,'CHRS',miss_novs,self.trim_len,mpath+'T_CHR.log')
		unmapped_reads = self.run_perm(mpath,tpath,'JXNS',miss_trim_chrs,self.trim_len,mpath+'T_JXN.log')
		lc = 0 
		for filename in os.listdir(mpath):
			if filename.endswith('.sam'):
				ll = ''
				for line in open(mpath+filename):
					line = line.split()
					if line[0] != ll: 
						ll = line[0] 
						lc += 1 			 

		self.total_mapped = lc 
		self.progress.end()




	def load_ref_config(self,config):
		self.ref_config = {} 
		try:		  config = open(config)
		except TypeError: config = config
		for line in config: 
			line = line.strip().split('=') 
			rkey,rname = line[0], line[1] 
			self.ref_config[rkey] = rname 
		config.close()

	def load_ref_indexes(self,config):

		self.progress.start_minor('Checking References...')
		self.ref_config = {} 
		try:		  config = open(config)
		except TypeError: config = config
		for line in config:
			line = line.strip().split('=') 
			rkey,rname = line[0], line[1] 

			self.ref_config[rkey] = rname

		self.index_key = dd(bool) 
		if 'IDXS' in self.ref_config:
			for f in os.listdir(self.ref_config['IDXS']):
				R,L,S = f.split('.index')[0].split('_')
				self.index_key[(R,int(L),int(S))] = self.ref_config['IDXS']+'/'+f

		config.close()	










		
	def link_read_config(self,config,rpath):

		self.progress.start_minor('Checking Reads........')
		self.read_config = {} 

		try:		  config = open(config) 
		except TypeError: config = config 


		for line in config: 
			line = line.strip().split('=') 
			rlen,rlist = int(line[0]), line[1] 

			local_list = open(rpath+rlist.split('/')[-1],'w')
			for line in open(rlist):
				line = line.strip()
				call(['ln','-sf',line,rpath+line.split('/')[-1]])
				local_list.write(rpath+line.split('/')[-1]+'\n')
			self.read_config[rlen] = rpath+rlist.split('/')[-1]

		self.trim_len, self.read_len = min(self.read_config.keys()), max(self.read_config.keys()) 



		
			

























	def load_read_config(self,config):


		self.read_config = {} 
		if type(config) == str: config = open(config) 
		for line in config: 
			line = line.strip().split('=') 
			rlen,rname = int(line[0]), line[1] 
			self.read_config[rlen] = rname 
			

		config.close() 

		self.trim_len, self.read_len = min(self.read_config.keys()), max(self.read_config.keys()) 
		
		










