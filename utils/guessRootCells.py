#!/usr/bin/env python
import pdb,sys,os
import argparse
import csv
import numpy as np
from sklearn.cluster import KMeans
import math

def Eentropy(Eg):
	P=np.histogram(Eg,bins=10)[0]
	P=[float(item)/sum(P) for item in P]
	P=[item for item in P if item>0]
	ee=0
	for i in P:
		ee+=-1*i*math.log(i)
	return ee
	
	
parser=argparse.ArgumentParser(description="Guess root cells for the model based on Entropy")
parser.add_argument('-i','--input',required=True,help='input single cell RNA-seq expression data')
parser.add_argument('-k','--ncluster', help="# of clusters for the first time point",default=2)
args = parser.parse_args()


with open(args.input,'r') as f:
	lf=f.readlines()
	ex=[item.strip().split("\t") for item in lf]
	
kT={}
FR=ex[0]
ex=ex[1:]
for i in ex:
	ti=float(i[1])
	if ti not in kT:
		kT[ti]=[i]
	else:
		kT[ti].append(i)

FT=sorted(kT.keys())[0]

ExFT=kT[FT]
rExFT=[[float(k) for k in item[3:]] for item in ExFT]
kmeans=KMeans(n_clusters=int(args.ncluster),random_state=0).fit(rExFT)
clusters=[]


dC={}
for i in range(len(ExFT)):
	ci=kmeans.labels_[i]
	if ci not in dC:
		dC[ci]=[ExFT[i]]
	else:
		dC[ci].append(ExFT[i])

LE=[]
for i in dC:
	Ei=dC[i]
	Ei=[item[3:] for item in Ei]
	Ei=[[float(k) for k in item] for item in Ei]
	EE=[]
	for g in range(len(Ei)):
		Eg=[item[g] for item in Ei]
		ee=Eentropy(Eg)
		EE.append(ee)
	EE=np.mean(EE)
	LE.append([i,EE])
LE=sorted(LE,reverse=True)
MLEC=[item[0] for item in dC[LE[0][0]]]


for i in ex:
	if i[0] in MLEC:
		i[1]=float(i[1])-1
	
ex=sorted(ex,key=lambda x: x[1])
ex=[FR]+ex	
ex=[[str(k) for k in item] for item in ex]
ex=["\t".join(item) for item in ex]
ex="\n".join(ex)

with open("%s.rootGuesses.txt"%(args.input),'w') as f:
	f.write(ex)





	
