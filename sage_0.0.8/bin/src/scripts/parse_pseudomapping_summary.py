#!/usr/bin/env python
from collections import defaultdict as dd

from subprocess import call


def input_error(msg):
	sys.stderr.write(msg+'\n')
	sys.exit() 




def get_line(line):

	rID = line[0] 
	ref = line[2] 
	if line[1] == '0': strand = '+'
	else:		   strand = '-'
	subs = int(line[11].split(':')[-1])

	rFeat = "~".join(rID.split('~')[1:5])
	fFeat = "~".join(ref.split("~")[1:5])
	if rFeat == fFeat: MATCH=True
	else:		   MATCH=False
	return MATCH,rID,strand,subs,rFeat,fFeat



	

########################################################################################################################################################################

if __name__ == '__main__':

	import sys,os,argparse,shutil
				
	class LineWrapRawTextHelpFormatter(argparse.RawDescriptionHelpFormatter):
    		def _split_lines(self, text, width):
        		text = self._whitespace_matcher.sub(' ', text).strip()

        		return _textwrap.wrap(text, width)	
        		#return _textwrap.wrap(text, width)	



	class MyParser(argparse.ArgumentParser):
		#formatter_class=argparse.RawTextHelpFormatteri
		formatter_class=LineWrapRawTextHelpFormatter
    		def error(self, detail=None,values=None,choice=None):
			
			if not detail: self.print_help()

			elif detail and values and choice:
				sys.stderr.write(str('InvalidOption:\n'))
				sys.stderr.write(str(values+" "+" ".join(choice)+'\n'))
				sys.stderr.write(str(detail)+'\n')

        			self.print_help()
			elif detail:
				sys.stderr.write(str(detail)+'\n')
				self.print_help()

        		sys.exit(2)
	
	parser=MyParser(formatter_class=argparse.RawTextHelpFormatter)
#	parser=MyParser(formatter_class=LineWrapRawTextHelpFormatter)
	help_format="%-40s %40s"


########################################################################################################################################################################
########################################################################################################################################################################
# GLOBAL PARENT OPTIONS ###############################################################################################################################################
	
	#Run = MyParser(add_help=False)
	#Paths.add_argument('--counts',type=argparse.FileType('r'),default=None,help='foo') 
	#Run.add_argument('--runpath', dest='runpath',action='store',type=str,default=None,help='Min and Max group size for k-cluster (eg. 3,10 )')
	#from argparse import RawTextHelpFormatter
	#parser = ArgumentParser(description='test', formatter_class=RawTextHelpFormatter)

	parser.add_argument('--outpath', dest='outpath',action='store',type=str,default=None,help='Output path for sage to run in (default local folder created)')
	parser.add_argument('--verbose',action='store_true', dest='verbose',default=False,help='Talkative?')
	parser.add_argument('--local',action='store_true', dest='local',default=False,help='Local Enviroment')
	parser.add_argument('--prefix', dest='prefix',action='store',type=str,default='sage',help='Output Prefix for Files')
#	Paths.add_argument('--test',action='store_true', dest='test',default=False,help='Talkative?')
	parser.add_argument('--cores', dest='cores',action='store',type=int,default=0,help='By default sage will run in parralel using (n-1) cores.  Set value to 1 if not desired')
	parser.add_argument( '--readlen', dest='readlen',action='store',type=int,default=100,help='Read length')
	parser.add_argument('file',type=argparse.FileType('r'),default=None,help='foo') 

	command_parsers = {} 
	













########################################################################################################################################################################
########################################################################################################################################################################
# RUN                          #########################################################################################################################################

#	os.path.abspath() 

	args = parser.parse_args()
	
	refLists = dd(lambda: dd(list)) 
	refCnts = dd(lambda: dd(int)) 
	refNeighbors = dd(lambda: dd(int))
	refSteal = dd(lambda: dd(int))
	refGive = dd(lambda: dd(int))

	for li,line in enumerate(args.file):
		
		lp = line.strip()
		if li == 0: 
			print lp 
			continue 
		line = line.split('|') 
		init =  "|".join(line[0:2])

		NAS =  [('NA','NA'),('NA','NA'),('NA','NA')]
		print init,
		for N in line[3::]:
			nS = N.split()
			nLen = len(nS)/2
			nDict = dd(int) 
			nData = []  
			nSum = 0 
			if nLen > 0:
				for i in range(0,len(nS),2):	
					if i < 5: nData.append((nS[i],nS[i+1]))
					nDict[nS[i]] = int(nS[i+1])
				nM,nSum = max(nDict.values()), sum(nDict.values())	

			nGOOD = nData + NAS
			print '|',nLen,nSum," ".join([" ".join(x) for x in nGOOD[0:2]]),	
		print 








		
