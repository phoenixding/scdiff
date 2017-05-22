#!/usr/bin/env python
"""
author: jun ding
date: June 8th,2016
function: calculate BIC score for clustering results
"""

import pdb,sys,os,math
from scipy.stats import spearmanr
from scipy.spatial.distance import *
from Distance import *

# get the center of a cluster
def getAvgEx(X):
		# get average experssion 
		# X: List of all Expression 		
		L=len(X[0])
		n=len(X)
		AE=[]
		for i in range(L):
			iAvg=sum([item[i] for item in X])/n
			AE.append(iAvg)
		return AE

#Inner-distance within cluster
def InnerDistance(A,dtype):
	mu=getAvgEx(A)
	S=0
	ct=0
	for i in A:
		di=eval(dtype)(i,mu)	
		S+=di
		ct+=1
	S=float(S)/ct
	return S

#Inter-distance between clusters 	
def InterDistance(A,B,dtype):
	mu1=getAvgEx(A)
	mu2=getAvgEx(B)
	di=eval(dtype)(mu1,mu2)
	return di
	
	
#Davies-Bouldin index
def DBI(X,Y,dtype="euclidean"):
	# X: Data
	# Y: Label
	YC=set(Y)
	K=len(YC)
	XC=[]
	for i in YC:
		Xi=[X[j] for j in range(len(Y)) if Y[j]==i]
		XC.append(Xi)
	
	DB=[]
	
	for i in range(len(XC)):
		x=XC[i]
		Sx=InnerDistance(x,dtype)
		DBx=[]
		for j in range(len(XC)):
			if i!=j:
				y=XC[j]
				Sy=InnerDistance(y,dtype)
				Mxy=InterDistance(x,y,dtype)
				DBx.append((Sx+Sy)/Mxy)
		DBx=max(DBx)
		DB.append(DBx)
	DB=sum(DB)/len(DB)
	return DB
	
						
#AIC
def AIC(X,Y,dtype="default"):
	YC=set(Y)
	K=len(YC)
	S=0
	for i in YC:
		CI=[X[j] for j in range(len(Y)) if Y[j]==i]
		mui=getAvgEx(CI)
		for j in CI:
			if dtype=="default":
				di=Distance(j,mui)
				dij=di.euclidean_distance()
			else:
				dij=eval(dtype)(j,mui)
			S+=dij
	FT=S
	R=len(X)
	D=len(X[0])
	AC=FT+4*K*D
	AC=math.log(AC)
	return AC	

#pham
def PhamAK(K,D):
	if K<1:
		print("error,K must larger than 1")
	elif K==1:
		return 1
	elif K==2:
		return 1-3.0/(4*D)
	else:
		return (5*PhamAK(K-1,D)+1)/6
	
def PhamSK(X,Y):
	YC=set(Y)
	K=len(YC)
	XC=[]
	S=0
	for i in YC:
		CI=[X[j] for j in range(len(Y)) if Y[j]==i]
		mui=getAvgEx(CI)
		for j in CI:
			di=euclidean(j,mui)
			S+=di
	return S

def PhamFK(X,Y,SKP):
	# X: [n_samples,n_features[
	# Y: labels_
	# SKP: measure of distance between points to its centroid for previous K
	#pdb.set_trace()
	SK=PhamSK(X,Y)
	YC=set(Y)
	K=len(YC)
	D=len(X[0])
	FK=1
	#pdb.set_trace()
	if (SKP!=0) and (K>1):
		ak=PhamAK(K,D)
		FK=SK/(SKP*ak)
	return [FK,SK]




	
