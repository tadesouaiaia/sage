#!/usr/bin/env python

import random
from collections import defaultdict as dd
from collections import Counter as cc
import sys
import os
from math import log 




class dot_log:
	#def __init__(self,VERBOSE,msg=' New_Module ',PARENT=False):
	def __init__(self,args):

		if args.verbose:	self.active = True 
		else:			self.active = False 

#		self.active = True

		self.intro = 'Module: '
		self.out = sys.stderr
		self.log = [[]] 


		if self.active: 	self.out.write('Sage Begins: '+'\n\n')

		self.blank = ' '
		self.sub_blank = ' '
		self.topic, self.subtopic = None, None 
		self.topics = [] 
		self.events = {} 

		self.status = None 


	def start_major(self,topic,msg=None,subtopics=[]): 


		if self.status == 'minor': self.end() 

		self.status = 'major'
		if self.active: 
			self.out.write(self.intro+topic+'\n') 
			self.blank = ''.join([' ' for x in range(len(self.intro))])

			if msg: 
				self.out.write(self.blank+msg+' '+', '.join([str(s) for s in subtopics])+'\n') 
				self.sub_blank = ''.join([' ' for x in range(len(self.intro+msg))])
			else:
				self.sub_blank = self.blank	

		self.log.append([topic])
		self.topic, self.subtopic = topic, None 
		self.topics.append(topic) 
		self.events[self.topic] = dd(list) 
		self.topic_counter = 0 
		


	def start_minor(self,topic,block_len=10,ENUMERATE=False):

		if self.status == 'minor': self.end() 

		self.status = 'minor'
		self.subtopic = topic
		self.log[-1].append(topic)
		self.counter = 0 
		self.ENUMERATE=ENUMERATE
		self.block_len,self.mp,self.dots,self.numbers = block_len,block_len,0,0 
		if block_len > 100: self.mp = 100 
		if self.active: self.out.write(self.sub_blank+topic+'...')



	def mark(self,dots=0):
		self.counter +=1 
		if self.active: 
			if dots > 0: self.out.write('.'.join(['' for x in range(dots)]))
			else:
				if self.counter % self.mp == 0: 
					self.out.write('.') 
					self.dots +=1
					if self.dots == 100:  self.mp = self.mp * self.mp 
					elif self.dots == 50: self.mp *= 4
					elif self.dots == 25: self.mp *= 2  
					elif self.dots == 10: self.mp += self.mp 

					if self.dots > 10 and self.dots % 10 == 0: self.mp += self.dots 

				

				if self.counter % self.block_len == 0 and self.ENUMERATE: 

					self.out.write('.'+str(self.counter)+'.')
					self.numbers += 1 
					if self.numbers == 4: self.block_len *= 2 
					elif self.numbers == 8: self.block_len *= 2 
					elif self.numbers % 10 == 0: self.block_len *= 2 
						
			
					


	def end(self): 



		if self.active and self.subtopic != None: self.out.write('Complete\n') 
		try:	
			self.events[self.topic][self.subtopic] = [self.counter]
			self.subtopic = None 
		except AttributeError:
			self.subtopic = None 
		except KeyError:
			self.subtopic = None
