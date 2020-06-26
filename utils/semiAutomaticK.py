import pdb,sys,os
import argparse
import csv
import numpy as np
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import matplotlib.colors as colors
#colors_list = list(colors._colors_full_map.values())
colors_list=list(colors.CSS4_COLORS.values())

parser=argparse.ArgumentParser(description="determine K using semi-automatic way")
parser.add_argument('-i','--input',required=True,help='input single cell RNA-seq expression data')
args = parser.parse_args()
f=open(args.input,'r')
lf=f.readlines()
f.close()
lf=[item.strip().split("\t") for item in lf]
times=[item[1] for item in lf[1:]]
dtimes=list(set(times))
ex=[[float(k) for k in item[3:]] for item in lf[1:]]
pca=PCA(n_components=50)
pex=pca.fit_transform(ex)
tsne=TSNE(n_components=2)
tex=tsne.fit_transform(pex)
colors=[colors_list[dtimes.index(item)%len(colors_list)] for item in times]
X=[item[0] for item in tex]
Y=[item[1] for item in tex]


#pdb.set_trace()
for i in dtimes:
	xi=[X[j] for j in range(len(X)) if times[j]==i]
	yi=[Y[j] for j in range(len(Y)) if times[j]==i]
	ci=[colors[j] for j in range(len(X)) if times[j]==i]
	plt.scatter(xi,yi,c=ci,label=i,s=1,alpha=0.5)
plt.legend()
plt.savefig("%s.tsne.pdf"%(args.input),dpi=300)



