#!/usr/bin/env python
"""
author: Jun Ding
usage: this moduleTestExample script demonstrates the usages of the modules of scdiff package
--scdiff.Graph
--scdiff.Cell
--scdiff.Clustering
--scdiff.viz

 
"""

import pdb,sys,os
from File import *
# please install scdiff package to use all modules given in this example test script. 
import scdiff
from scdiff import *

# reading example cells ...
AllCells=[]
print("reading cells...")
with open("example.E",'r') as f:
	line_ct=0
	for line in f:
		if line_ct==0:
			GL=line.strip().split("\t")[3:]
		else:
			line=line.strip().split("\t")
			iid=line[0]
			ti=float(line[1])
			li=line[2]
			ei=[round(float(item),2) for item in line[3:]]
			ci=scdiff.Cell(iid,ti,ei,li,GL)
			AllCells.append(ci)
		line_ct+=1
		print('cell:'+str(line_ct))

# clustering cells using scdiff.Clustering module 
print("testing scdiff.Clustering  ...")
ClusteringExample=scdiff.Clustering(AllCells,'auto',None)
clusters=ClusteringExample.performClustering()


print("testing scdiff.Graph  ...")
# creating graph using scdiff.Graph module. 
g1=scdiff.Graph(AllCells,"example.tf_dna",'auto')

print ("testing scdiff.Cluster...")
cluster1=scdiff.Cluster([item for item in AllCells if item.T==14],14,'C1')

print ("testing scdiff.Path ...")
p1=scdiff.Path(g1.Nodes[0],g1.Nodes[1],g1.Nodes,g1.dTD,g1.dTG,g1.dMb)

print ("testing scdiff.viz ...")
# visualizing graph using scdiff.viz module 
os.mkdir("e1_out")
scdiff.viz("example",g1,"e1_out")

print("module testing done!")
