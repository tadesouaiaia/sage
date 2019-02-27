#!/usr/bin/env python
from collections import defaultdict as dd
from collections import Counter as cc

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


def rswap(r):
	rsplit = r.split('~') 

	rnick =  ";".join([rsplit[2].split('.')[0],rsplit[0],rsplit[-1]])
	return rnick

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


	refDouble = dd(lambda: dd(list))
	X= args.file.readline().split()

	MATCH,rID,strand,subs,rFeat,fFeat = get_line(X) 
	myID,myRef = rID , rFeat
	myData = [] 
	ITER=True

	pp = 0 
	while ITER:
		pp+=1
		while myID == rID:
			myData.append([MATCH,strand,subs,rFeat,fFeat])	
			try: 
				MATCH,rID,strand,subs,rFeat,fFeat = get_line(args.file.readline().split())
			except IndexError: 
				ITER=False
				break		

		res =  [m[0] for m in myData] 
		if True in res: 
			refCnts[myRef]['total'] +=1
			rDub = {z: 0 for z in range(0,11)}
			if not False in res:
				refCnts[myRef]['uniq']+=1 
				refLists[myRef]['dists'].append(11) 
				refLists[myRef]['locs'].append(0)
			else:

				nData,nDict = sorted([(m[2],m[4]) for m in myData if not m[0]]) , {} 
				for nd in nData:
					if nd[1] not in nDict: nDict[nd[1]] = nd[0] 

				refLists[myRef]['dists'].append(nData[0][0])
				refLists[myRef]['locs'].append(len(nDict.keys()))

				for x,y in nDict.items():	
					rDub[y] += 1 
					refNeighbors[myRef][x] += 1
					refNeighbors[x][myRef] += 1
					refGive[myRef][x] += 1
					refSteal[x][myRef] += 1

			for z in range(11):
				refDouble[myRef][z].append(sum([rDub[w] for w in range(0,z+1)]))


		myID,myRef,myData = rID, rFeat, [] 

#		if pp > 1000: break 


	for r in refCnts.keys():



#		print r,refCnts[r]['total'],refCnts[r]['uniq'],refLists[r]['dists'] 
		tC,uC,dists,locs,neighbors =   refCnts[r]['total'],refCnts[r]['uniq'],refLists[r]['dists'],refLists[r]['locs'],refNeighbors[r]

		nbSrt = sorted(refNeighbors[r].items(),key=lambda x: x[1],reverse=True)
		nGSrt = sorted(refGive[r].items(),key=lambda x: x[1],reverse=True)
		nSSrt = sorted(refSteal[r].items(),key=lambda x: x[1],reverse=True)

		UNQ_DICT = {} 
		NB_DICT = {} 
		DCX = dd(int) 
		for (a,b) in cc(dists).items(): DCX[a] = b 
		for i in range(11): 
			#print i,DCX[i]
			#print i,sum([DCX[j] for j in range(0,i+1)]),DCX[i],len(dists) 
			#print i,1-sum([DCX[j] for j in range(0,i+1)])/float(len(dists))
			UNQ_DICT[i] = 1-sum([DCX[j] for j in range(0,i+1)])/float(len(dists))
			NB_DICT[i] =  sum(refDouble[r][i])/float(len(dists))


		UNQ_STRINGS = " ".join([str(round(UNQ_DICT[i],4)) for i in range(11)])
		NB_STRINGS =  " ".join([str(round(NB_DICT[i],4)) for i in range(11)])


		rnick = rswap(r) 


		print rnick,'STATS',tC,uC,round(uC/float(tC),5),round(sum(dists)/float(len(dists)),5),'|',UNQ_STRINGS,'|',NB_STRINGS
		print rnick,'NBS',len(dists),
		nbd = [] 
		for X,Y in zip([nbSrt,nGSrt,nSSrt],['ALL','GIVE',"TAKE"]):

			xSum =  sum([x[1] for x in X]) / float(len(dists))
			xLen =  len(X)/float(len(dists)) 
			xTmp =  [a for b in  [(rswap(x[0]),str(round(x[1]/float(len(dists)),2))) for i,x in enumerate(X) if i < 2] for a in b]+['NA','NA','NA','NA']

			print Y,round(xLen,4),round(xSum,4)," ".join(xTmp[0:2]),

		print ""
		continue
		sys.exit() 



		for X,Y in nbSrt:
			print X.split('~')[0]+'~'+X.split('~')[-1],Y,
		print '|',

		for X,Y in nSSrt:
			print X.split('~')[0]+'~'+X.split('~')[-1],Y,
		print '|',
		
		for X,Y in nGSrt:
			print X.split('~')[0]+'~'+X.split('~')[-1],Y,
		print '|',
		print ""










