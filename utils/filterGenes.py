#!/usr/bin/env python
import pdb,sys,os
import argparse
import csv
import numpy as np


parser=argparse.ArgumentParser(description="Filter non-informative genes for large datasets")
parser.add_argument('-i','--input',required=True,help='input single cell RNA-seq expression data')
parser.add_argument('-s','--std',help="remove genes with standard deviation smaller than the specified cutoff",default=0.5)
parser.add_argument('-n','--ngenes', help="keep the top n genes with the largest varience",default=5000)
parser.add_argument('--setime', help="set all time point of all cells to a given number")
parser.add_argument('-o','--output', help='output for the filtered single cell RNA-seq expression data')
args = parser.parse_args()
try:
	ng=int(args.ngenes)
	ns=float(args.std)
except:
	print("check your input! -s -n")
	sys.exit(0)
setime=args.setime 
output=args.output 

with open(args.input) as f:
	keptcols=[]
	reader=csv.reader(f,delimiter="\t")
	ncol=len(next(reader))
	stdlist=[]
	#==========
	span=5000
	st=3
	ed=st+span
	
	while (st<ncol):
		f.seek(0)
		cols=[]
		for row in reader:
			cols.append(row[st:ed])
		cols=cols[1:]
		rows=[[float(k) for k in item] for item in cols]
		cols=np.array(rows).T.tolist()
		stds=[np.std(item) for item in cols]
		stdlist+=stds
		st=ed
		ed=st+span
		#print(ed)
	stdcut=sorted(stdlist,reverse=True)[ng]
	scols=range(3)+[3+k for k in range(len(stdlist)) if stdlist[k]>stdcut]
	
	f.seek(0)
	
	if output == None:
		for row in reader:
			sline=[row[item] for item in scols]
			sline[1]=sline[1] if setime==None else setime 
			sline="\t".join(sline)+'\n'
			print(sline)
	else:
		with open(output, "w") as outputfile:
			for row in reader:
				sline=[row[item] for item in scols]
				sline[1]=sline[1] if setime==None else setime 
				sline="\t".join(sline)+'\n'
				outputfile.write(sline)
