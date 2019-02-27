#!/usr/bin/env python
import pandas as pd
import sys
import os
import random
from collections import defaultdict as dd
from collections import Counter as cc
import scipy.stats as stats
from scipy.stats import variation as coVar 

from random import random
import numpy as np

import random
from math import fabs
from scipy.stats import pearsonr as pearsonr
#from scipy.stats import spearmanr as spearmanr
import seaborn
from math import log
import math
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from random import shuffle
from sklearn.cluster import KMeans	
from sklearn.cluster import KMeans
from sklearn.neighbors import KernelDensity
from sklearn.preprocessing import MinMaxScaler
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle as Rect
sns.set(color_codes=True)
import matplotlib as mpl
label_size = 5
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size 
import warnings
warnings.filterwarnings("ignore")

def get_colors(inp, colormap, vmin=None, vmax=None):
    norm = plt.Normalize(vmin, vmax)
    return colormap(norm(inp))




def mad_based_outlier(points, thresh=3.5):
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

#    print get_colors(modified_z_score, plt.cm.jet)
    return modified_z_score > thresh, modified_z_score, get_colors(modified_z_score, plt.cm.seismic) 
#    return modified_z_score > thresh, modified_z_score, get_colors(modified_z_score, plt.cm.jet) 



def plot_stats(sdict):


	sp = subplot(2,2)



 
	#[[155.0, 154.0, 0.99355, 10.99355], [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.9935], [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.9935]]
	names = sdict.keys()
	basic = [x[0] for x in sdict.values()]
	urates = [x[1] for x in sdict.values()]
	nrates = [x[2] for x in sdict.values()]


	sp.add_hist({n:b[2] for n,b in zip(names,basic)}).update({'title': 'Percent of Unique Reads (no neighbors within distance 10)'})
	sp.add_hist({n:b[3] for n,b in zip(names,basic)}).update({'title': 'Avg Distance to Nearest Neighbor','xticklabels': [[0,4,8,10.5],['0','4','8','>10']]})

	X = range(len(urates[0]))
	myX,myY,myN=[],[],[]
	myX_bad,myX_reg,myX_good = [],[],[]
	myY_bad,myY_reg,myY_good = [],[],[] 
	for p,ur in enumerate(urates):

		yd = [np.random.normal(z,0.0000001+z*0.001) for z in ur] 
		xd = [np.random.normal(z,0.5+z*0.0000001) for z in X]
		if ur[0] == 0: 
			myY_bad.extend(yd)
			myX_bad.extend(xd)
		elif ur[-1] == 1.0: 
			myY_good.extend(yd)
			myX_good.extend(xd)
	
		else:
			myY_reg.extend(yd)
			myX_reg.extend(xd) 


		myY.extend(xd)
		myX.extend(yd)

	sp.add_scatter_trend(myX_reg,myY_reg,clr='blue',alpha=0.025) #.update({'title': 'uniq-rates at different mismatch thresholds'})
	sp.add_scatter_pts(myX_good,myY_good,'GREEN')
	sp.add_scatter_pts(myX_bad,myY_bad,'RED').update({'title': 'uniq-rates at different mismatch thresholds'})





	myX_bad,myX_reg,myX_good = [],[],[]
	myN_bad,myN_reg,myN_good = [],[],[] 
	for p,ur in enumerate(nrates):
#		myN.extend([np.random.normal(z,0.0000001+z*0.001) for z in ur]) 
		yd = [np.random.normal(z,0.0000001+z*0.001) for z in ur] 
		xd = [np.random.normal(z,0.5+z*0.0000001) for z in X]

		if ur[-1] == 0: 
			myN_good.extend(yd)
			myX_good.extend(xd)
			#myN_good.extend([np.random.normal(z,0.0000001+z*0.001) for z in ur]) 
			#myX_good.extend([np.random.normal(z,0.00001+z*0.02) for z in X]) 
		else:
			myN_reg.extend(yd)
			myX_reg.extend(xd) 
			#myN_reg.extend([np.random.normal(z,0.0000001+z*0.001) for z in ur]) 
			#myX_reg.extend([np.random.normal(z,0.00001+z*0.02) for z in X]) 

	sp.add_scatter_trend(myX_reg,myN_reg,clr='blue',alpha=0.25)#.update({'title': 'avg-collisions at different mismatch thresholds'})
	sp.add_scatter_pts(myX_good,myN_good,'GREEN').update({'title': 'avg-collisions at different mismatch thresholds'})

	plt.suptitle('Analysis of Collisions for Simulated Reads',fontsize=15,fontweight='bold') 
	plt.show() 










class subplot:
        def __init__(self,xLen,yLen,options=None,key={}):

		label_size = 9.5
		matplotlib.rcParams['xtick.labelsize'] = label_size 
		matplotlib.rcParams['ytick.labelsize'] = label_size 
	
                self.fig = matplotlib.pyplot.gcf()
                self.fig.set_size_inches(18.5, 9.5)
		self.xLen, self.yLen = xLen,yLen
		self.xLoc, self.yLoc = 0,0 
		self.options = options
		self.trend = False


		if 'titlepos' in key: 
			self.titleX,self.titleY = key['titlepos']
		else:
			self.titleX,self.titleY = 0.0,1.05


		self.ax_key = {} 
		self.ax = plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 1)
		self.ax_key[(self.xLoc,self.yLoc)] = self.ax 

	def update(self,key={}):

		if 'clear_axes' in key: 
			self.ax.set_xticks([]) 
			self.ax.set_yticks([]) 

		elif 'xticklabels' in key:
			ticks,labs = key['xticklabels']
			self.ax.set_xticks(ticks) 
			self.ax.set_xticklabels(labs,horizontalalignment='left') 


		if 'xlab' in key: 	self.ax.set_xlabel(key['xlab'])
		if 'ylab' in key: 	self.ax.set_ylabel(key['ylab'])

		if not self.trend:
			if 'title' in key: 	self.ax.set_title(key['title'],fontweight='bold',x=self.titleX,y=self.titleY,horizontalalignment='left',fontsize=12.5)


		self.yLoc += 1 
		if self.yLoc == self.yLen: 
			self.xLoc,self.yLoc = self.xLoc+1, 0 
		if self.xLoc >= self.xLen: return False
		
		self.ax = plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 1)
		self.ax_key[(self.xLoc,self.yLoc)] = self.ax 
		self.trend=False





		return True 



	def skip_row(self):

		#if self.yLoc != 0: 
		self.xLoc +=1 
		self.yLoc = 0 

		if self.xLoc >= self.xLen: return False
		self.ax = plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 1)
		self.ax_key[(self.xLoc,self.yLoc)] = self.ax 
		return True



	def save(self,f_name,key={}):

	        #plt.subplots_adjust(left=0.02, bottom=0.02, right=0.80, top=0.90,wspace=0.15,hspace=0.5)


		if 'axis_off' in key: plt.axis('off') 
		if 'title' in key: plt.suptitle(key['title'],fontsize=20,fontweight='bold')
		elif self.options and 'title' in vars(self.options).keys() and self.options.title != None:	plt.suptitle(self.options.title)
               	self.fig.savefig(f_name, dpi=100) 

		if self.options.show: 
			plt.show() 



	def add_hist(self,h_data):

		x,y = h_data.keys(), sorted(h_data.values())

	#	a = np.hstack(h_data)
		if len(self.ax.hist(y,bins='auto')[0]) < 2:

			self.ax.clear() 
			self.ax.hist(np.hstack(y),bins=min(4,len(cc(y)))) 


		yMin, yLim = self.ax.get_ylim()
		xMin, xLim = self.ax.get_xlim()
		HI,LO=False,False
		if yMin == 0: 
			h_srt = sorted(y) 
			out_bool, out_scrs, out_colors =  mad_based_outlier(np.array(h_srt)) 
			g_color = out_colors[sorted([(out_scrs[oi],oi) for oi in range(len(out_scrs))])[0][1]]
			for i,h in enumerate(h_srt): 
				if not out_bool[i]:
					self.ax.scatter(h,yLim*1.025,alpha=0.7,color=g_color,clip_on=False) 
				else: 
					self.ax.scatter(h,yLim*1.025,alpha=0.7,color=out_colors[i],clip_on=False) 
				#	if i > 5 and (out_scrs[i]/(out_scrs[i-1]+out_scrs[i-2])) > 0.66:  HI=True 
				#	if HI: 			self.ax.text(h,yLim*1.030,x[i].name.split(";")[-1])


			self.ax.set_ylim(yMin,yLim) 
			self.ax.set_xlim(xMin,xLim) 


		return self




	def add_pca_data(self,pca_pts,key={}):

		#sns.set(rc={'axes.facecolor':'black', 'figure.facecolor':'cornflowerblue'})
		if 'components' in key: comps = key['components']
		else:			comps = [0,1] 
		if 'type' in key:
			if key['type'] == 'tsne':	
				axis_labs,title = ['TS'+str(c+1) for c in comps], 'TSNE'
			elif key['type'] == 'kca': 
				axis_labs,title = ['KC'+str(c+1) for c in comps], 'KCA'
			elif key['type'] == 'mds': 
				axis_labs,title = ['MD'+str(c+1) for c in comps], 'MDS'
			else:
				axis_labs,title = ['X'+str(c+1) for c in comps], key['type'].upper() 
		else:
			title = 'PCA'
			if 'vars' in key: 				axis_labs = ['PC'+str(c+1)+':  '+str(round(100*key['vars'][c],2))+'%' for c in comps]
			else:						axis_labs = ['PC'+str(c+1) for c in comps]

		for i,p in enumerate(pca_pts):
			if 'colors' in key and key['colors'][i] != 'white': 
				if 'sizes' in key: self.ax.scatter(p[0],p[1],s = key['sizes'][i],color=key['colors'][i],alpha=0.6)
				else: 		   self.ax.scatter(p[0],p[1],color=key['colors'][i])
			else: 		    self.ax.scatter(p[0],p[1])

		self.ax.set_xlabel(axis_labs[0]) 	
		self.ax.set_ylabel(axis_labs[1])

		if 'zoom' in key: 
			pX = sorted([p[0] for p in pca_pts])
			pY = sorted([p[1] for p in pca_pts])
			xQ1,xQ3 = int(len(pX)*0.25), int(len(pX)*0.75)
			yQ1,yQ3 = int(len(pY)*0.25), int(len(pY)*0.75)
			iqX = pX[xQ3]-pX[xQ1]
			iqY = pY[yQ3]-pY[yQ1]  
			self.ax.set_xlim(pX[xQ1]-(iqX*2.5),pX[xQ3]+(iqX*2.5))
			self.ax.set_ylim(pY[yQ1]-(iqY*2.5),pY[yQ3]+(iqY*2.5))



		if 'title' in key: self.ax.set_title(key['title']) 
		else:		   self.ax.set_title(title) 
		return self		
		


	def change_limits(self,key):
		x0,x1 = self.ax.get_xlim()
		y0,y1 = self.ax.get_ylim()
		if 'x0' in key: x0 = key['x0'] 
		if 'x1' in key: x1 = key['x1'] 
		if 'y0' in key: y0 = key['y0'] 
		if 'y1' in key: y1 = key['y1'] 
		self.ax.set_xlim(x0,x1)
		self.ax.set_ylim(y0,y1)


	def add_lines(self,x,y,xlabel=None,ylabel=None,clr=None):
		if not clr: self.ax.plot(x,y,zorder=1)
		else: 	    self.ax.plot(x,y,color=clr,zorder=1) 
		if xlabel: self.ax.set_xlabel(xlabel) 
		if ylabel: self.ax.set_ylabel(ylabel) 
		return self

	def add_line(self,x,y,key={}):
		lw,alpha = 1.0, 1.0 
		if 'lw' in key.keys():    lw = key['lw'] 
		if 'alpha' in key: alpha = key['alpha'] 
		if 'color' in key: self.ax.plot(x,y,linewidth=lw,color=key['color'],zorder=1,alpha=alpha) 
		else:  self.ax.plot(x,y,linewidth=lw,zorder=1,alpha=alpha)
		if 'xlab' in key: 	self.ax.set_xlabel(key['xlab'])
		if 'ylab' in key: 	self.ax.set_ylabel(key['ylab'])

		return self

	def scatter_pts(self,X,Y,key={}):
		alpha,size,mark =  1.0, 20,'o'
		if 'size' in key.keys():    size = key['size'] 
		if 'alpha' in key: alpha = key['alpha'] 
		if 'mark' in key: mark= key['mark'] 
		if 'yjitter' in key: Y = [np.random.normal(y,(0.0195)) for y in Y]
		

		if 'color' in key: self.ax.scatter(X,Y,s=size,marker = mark, color=key['color'],zorder=1,alpha=alpha) 
		else:  self.ax.scatter(X,Y,zorder=1,alpha=alpha,s=size)
		return self

	def add_scatter_trend(self,x,y,xlabel=None,ylabel=None,clr=None,alpha=0.99):



		if clr and len(clr.split(',')) == 2: 
			self.ax = sns.regplot(x=np.array(x),y=np.array(y),scatter_kws={"color": clr.split(',')[0]}, line_kws={"color": clr.split(',')[0]})
		else:
			if clr == None: clr = 'green'
			self.ax = sns.regplot(x=np.array(x),y=np.array(y),scatter_kws={"color": clr,"alpha": alpha}, line_kws={"color": clr,'lw': 1})






#		 pd.DataFrame(np.array([np.array([values[j][i] for j in range(len(values))]) for i in range(len(idxs))]), columns = labels)

#		self.ax = sns.lmplot(x="x", y="y", data=data);

		return self


	def add_scatter_data(self,x,y,xlabel=None,ylabel=None,clr=None):




		if clr and len(clr.split(',')) == 2: 
			self.ax = sns.regplot(x=np.array(x),y=np.array(y),scatter_kws={"color": clr.split(',')[0]}, line_kws={"color": clr.split(',')[1]})
		else:
			if clr == None: clr = 'k'
			self.ax = sns.regplot(x=np.array(x),y=np.array(y),color=clr)
		if xlabel: self.ax.set_xlabel(xlabel) 
		if ylabel: self.ax.set_ylabel(ylabel) 

		return self
	
	def add_scatter_pts(self,X,Y,notes='None'):
		notes = notes.split(',') 
		if 'RED' in notes: clr = 'red' 
		elif 'PURPLE' in notes: clr = 'purple' 
		elif 'BLUE' in notes: clr = 'blue'
		elif 'GREEN' in notes: clr = 'green'
		else:		      clr = 'k'  
		if 'XY_JITTERS' in notes:
			xJ = [np.random.normal(x,(0.0000001+x*0.1)) for x in X]
			yJ = [np.random.normal(y,(0.0000001+y*0.2)) for y in Y]
			self.ax.scatter(xJ,yJ,color=clr)

		else:
			self.ax.scatter(X,Y,color=clr,alpha=0.20) 

		return self
	




	def add_labels(self,title,xlabel,ylabel):
		self.ax.set_title(title) 
		if xlabel: self.ax.set_xlabel(xlabel) 
		if ylabel: self.ax.set_ylabel(ylabel) 

	def add_legend(self,labels,colors):
		labs,items = [],[] 
		for a,b in zip(labels,colors):
			labs.append(a) 
			items.append(Rect((0,0),1,1,fc=b))

		self.ax_key[(0,self.yLen-1)].legend(items,labs,ncol=1, loc='upper right',   bbox_to_anchor=(2.0,1.0),fontsize=12)


#		plt.legend(items,labs)
		#plt.legend(items,labs,loc = 'upper left',ncol = len(labs)/1,bbox_to_anchor=(-2.5,-.16),fontsize=10)
		
	def add_skree(self,pca_vars,key={}):

		if 'color' in key: 	
			self.ax.plot([sum(pca_vars[0:i]) for i in range(len(pca_vars))],color=key['color'])
		else:
			self.ax.plot([sum(pca_vars[0:i]) for i in range(len(pca_vars))])

		self.ax = plt.subplot2grid((self.xLen,self.yLen), (self.xLoc,self.yLoc), rowspan = 1, colspan = 1)
		return self	



		
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
	
	myDict = {} 
	for li,line in enumerate(args.file):
		
		lp = line.strip()


		if li == 0: 
			ls = line.split() 
			if ls[1] == 'STATS': STATS=True
			else:		     STATS=False 
			continue 


		if STATS:
			uDict,nDict = {},{} 
			line = line.split('|') 
			gene = line[0].split()[0] 
			reads,uniq,uniqRate,avgDist = [float(x) for x in line[0].split()[2::]]
			
			uRates = [float(x) for x in line[1].split()]
			avgNBS = [float(x) for x in line[2].split()]

			
			myDict[gene] = [[reads,uniq,uniqRate,avgDist],uRates,avgNBS]



#		if li > 1000: break 


	if STATS:
		plot_stats(myDict) 







		
