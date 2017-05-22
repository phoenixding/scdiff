#!/usr/bin/env python
import pdb,sys,os
from scipy.stats import spearmanr

class Distance:
	def __init__(self,fromx,toy):
		if len(fromx)!=len(toy):
			print('error! the length of x and y must match')
			sys.exit(0)
			
		self.fromx=fromx
		self.toy=toy
		
	def euclidean_distance(self):
		S=0
		for i in range(len(self.fromx)):
			S+=(self.fromx[i]-self.toy[i])**2
		return S
		
	def spearmanr_distance(self):
		S=1-spearmanr(self.fromx,self.toy)[0]
		return S
		
	def manhattan_distance(self):
		S=0
		for i in range(len(self.fromx)):
			S+=abs(self.toy[i]-self.fromx[i])
		return S
	
		
