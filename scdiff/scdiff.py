#!/usr/bin/env python
# author: Jun Ding
# Date: Feb.  10th, 2018
# import system modules
import pdb,sys,os,random
import numpy as np
np.seterr(divide='ignore', invalid='ignore')
import warnings
warnings.filterwarnings("ignore")
import math
import copy
import pickle
import argparse
import functools

# import Statistics modules
from scipy.stats import spearmanr
from scipy.stats import pearsonr
from scipy.stats import ttest_ind
from StatTest import pbinom
from scipy.stats import ranksums
from scipy.stats import wilcoxon
from StatTest import HyperGeometricTest as HyperTest

# import File modules
from File import TabFile
from File import MotifFile
from File import FastaFile

# import Machine learning modules
from sklearn.cluster import SpectralClustering
from sklearn.cluster import MiniBatchKMeans
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import KMeans
from sklearn.cluster import Birch
from imblearn.over_sampling import SMOTE

from sklearn.metrics import silhouette_score
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import LogisticRegressionCV
from sklearn.decomposition import PCA
from KF2 import smootherMultiple # import the Kalman Smoother (handle multiple observations)
from KF2 import KalmanFilterInd  # import the kalman filter class
from ClusteringMetric import DBI,AIC,PhamFK

# import spatial 
from scipy.spatial.distance import pdist,squareform

# import graph drawing module
#import matplotlib.pyplot as plt
from viz import viz

import gc
import pkg_resources 

#-----------------------------------------------------------------------------------------------------------------------
# cell
class Cell:
        def __init__(self, Cell_ID, TimePoint, Expression,typeLabel,GL):
                self.ID=Cell_ID
                self.T=TimePoint
                self.E=Expression
                self.Label=None # Label for clustering purpose
                self.typeLabel=typeLabel
                self.GeneList=GL
        def distanceTo(self,Ancestor):
                dta=spearmanr(self.E, Ancestor.E)[0]
                return dta

#-----------------------------------------------------------------------------------------------------------------------
# Clustering
class Clustering:
	def __init__(self,AllCells,kc,largeType=None):
		self.cells=AllCells                               # all cells for clustering
		self.dET=self.getTimeEx()                         # cells for each time point
		self.KET=sorted(self.dET.keys())                  # time points sorted
		self.largeType=largeType                          # whether to use the large dataset model for clustering
		self.kc=kc                                        # clusters
		
	# get clustering parameters
	def getClusteringPars(self):
		if self.largeType==None:
			self.affMatrix=self.getAffMatrix(self.dET)            # get affinity matries (map: time->matrix) normal model 
		else:
			self.affMatrix=self.getAffMatrixLarge(self.dET)
			
		dCK=self.getdCK(self.kc)                                      # number of clusters for each time point, director
		dBS=self.determineSeed(dCK)                                   # best seed for each clustering (1 clustering for each time point)
		return [dCK,dBS]

#---------------------------------------------------------------
	# affinity matrix
	def getAffMatrix(self,det):
		#Cells->Affinity
		def calculateAffinity(X):
			# X: expression
			XZ=[[0 for i in range(len(X))] for j in range(len(X))]
			for i in range(len(X)):
				for j in range(len(X)):
					if i==j:
						XZ[i][j]=1
					elif i<j:
						XZ[i][j]=spearmanr(X[i],X[j])[0]
					else:
						XZ[i][j]=XZ[j][i]


				print("cell:%s"%(i))

			XZ=np.array(XZ)
			return XZ

		# calculat affinity matrix
		print('calculating affinity matrix ...')
		dTM={}
		for T in self.KET:
			CT = det[T]
			dTM[T] = calculateAffinity([item.E for item in CT])
		return dTM
	#---------------------------------------------------------------
	# affinity matrix (largeType)
	def getAffMatrixLarge(self,det):
		print('calculating affinity matrix ...')
		#cells->affinity
		dTM={}
		for T in self.KET:
			CT=det[T]
			X=np.array([item.E for item in CT])
			pca=PCA(n_components=10)
			X=pca.fit_transform(X)
			#pdb.set_trace()
			dTM[T]=X
			print(T)
		return dTM
	
	# Get number of clusters for each point (dCK)
	def getdCK(self,kc):
		#return {14.0:1,16.0:2,18.0:4}
		try:
				KCF = TabFile(kc).read('\t')
				dCK={float(item[0]):int(item[1]) for item in KCF}
		except:
				dCK = self.determineK(K0=1)
		return dCK

	
	# Affinity-> Distance
	def affinity2Distance(self,X):
		XX=[]
		for i in range(len(X)):
			XX.append([1-item for item in X[i]])
		XX=np.array(XX)
		return XX

	
	def getTimeEx(self):
                dET = {}
                for i in self.cells:
                        if i.T not in dET:
                                dET[i.T] = [i]
                        else:
                                dET[i.T].append(i)
                return dET

	# determine clustering parameter K for each time point.  
	def determineK(self,K0=1):
		# --------------------------------------------------------------
		# best K based on Silhouette score or other metrics
		# nested function: which can access variable inside the outside function
		#---------------------------------------------------------------
		def bestK(T):
			CS = detPeak(K, ST, 0.1)
			CM = detPeak(K, MT, 0.1, type='min')
			AK=CS+CM
			if len(AK)>0:
				CAK={item:AK.count(item) for item in set(AK)}
				SortKey=sorted(CAK.keys(),key=lambda item:(CAK[item],-1*float(item)), reverse=True)
				SortKey=[item for item in SortKey if CAK[item]>=CAK[SortKey[0]]-1]
				#SortKey=[item for item in SortKey if CAK[item]>=CAK[SortKey[0]]-0]
				ATP=[AT[K.index(item)] for item in SortKey]
				bk=AT.index(min(ATP))
			else:
				try:
					MTC=[MT[item] for item in range(len(K)) if K[item]>=pn]
					minMTC=min(MTC)
					ATC=[AT[item] for item in range(len(K)) if K[item]>=pn]
					minATC=min(ATC)
					SK=[MT.index(minMTC), AT.index(minATC)]
					bk=max(SK)
				except:
					return 1
			return K[bk]
		#----------------------------------------------------------------
		KET=self.KET
		dCK = {}
		dCK[KET[0]] = K0
		K = range(2, 10) if self.largeType==None else range(2,8)
		print("learning K...")
		
		#-----------------------------------------------------------------
		# try many iteration, choose the one with most votes
		
		Loop=50 if self.largeType==None else 5
		CK=[]
		L=len(self.cells[0].E)
		
		for li in range(Loop):
			pn = K0
			CKT=[]
			time_counter=0
			for T in KET[1:]:
				K = range(2, 10) if self.largeType==None else range(2,7)
				CT=self.dET[T]
				ST = []		# silhouette score
				MT=[]		# DBI score
				AT=[]		# AIC score
				FT=[]		# Pham Score
				if self.largeType=='1' or self.largeType=='True':
					X=copy.deepcopy(self.affMatrix[T])
					K=[item for item in K if item<len(X)]
					FT.append(PhamFK(X,[0 for i in range(len(X))],0))
					for n in K:
						gc.collect()
						SC = KMeans(n_clusters=n)
						SC.fit(X)	
						Y = SC.labels_
						sscore=silhouette_score(X,Y,metric="cosine")
						sscore = max(sscore, 0)
						mscore=DBI(X,Y,"cosine")
						ascore=AIC([item.E for item in CT],Y)
						[fk,sk]=PhamFK(X,Y,FT[-1][-1])
						ST.append(sscore)
						MT.append(mscore)
						AT.append(ascore)
						FT.append([fk,sk])
						print(n)

				else:
					X=copy.deepcopy(self.affMatrix[T])
					K=[item for item in K if item<len(X)]
					DX=self.affinity2Distance(X)
					FT.append(PhamFK(X,[0 for i in range(len(X))],0))
					for n in K:
						gc.collect()
						SC = SpectralClustering(n_clusters=n)
						if len(X)<=n:
							break
						SC.fit(X)
						Y = SC.labels_
						sscore = silhouette_score(DX, Y, metric="precomputed")
						sscore = max(sscore, 0)
						mscore=DBI(X,Y)
						ascore=AIC([item.E for item in CT],Y)
						[fk,sk]=PhamFK(X,Y,FT[-1][-1])
						ST.append(sscore)
						MT.append(mscore)
						AT.append(ascore)
						FT.append([fk,sk])
						print(n)
				#pdb.set_trace()
				# of clusters in previous time point
				pn = bestK(T)
				"""
				fk2=FT[1][0]
				#pdb.set_trace()
				if (time_counter==0 and pn==2 and fk2>1):
					pn=1
				"""
				CKT.append(pn)
				time_counter+=1
			CK.append(CKT)
			print('%s'%((li+1.0)/Loop*100)+'%')		
		CK=sorted([[CK.count(item),item] for item in CK],reverse=True)
		CK=[item[1] for item in CK]
		CKK=CK[0]
		for i in range(len(KET[1:])):
			T=KET[1:][i]
			dCK[T]=CKK[i]
		#pdb.set_trace()
		return dCK

	def determineSeed(self,dCK):
		#return {14.0:0,16.0:0,18.0:0}
		print("learning clustering seeds...")
		dBS = {}  # Best seeds
		KET=self.KET
		NSEEDS=100 if self.largeType ==None else 1 #100
		SPECTRALIMIT=100
		
		for T in KET[1:]:
			try:
				CT = self.dET[T]
				CKi = dCK[T]
				SS=[]
				if self.largeType=='1' or self.largeType=='True':
					X=copy.deepcopy(self.affMatrix[T])
					SEEDS = range(NSEEDS)
					for s in SEEDS:
						SC = KMeans(n_clusters=CKi)
						SC.fit(X)
						Y = SC.labels_
						sscore = silhouette_score(X, Y)
						SS.append(sscore)
						print("seeds:"+str(s))
					sbest = SEEDS[SS.index(max(SS))]
					dBS[T] = sbest
				else:
					X=copy.deepcopy(self.affMatrix[T])
					DX=self.affinity2Distance(X)
					SEEDS = range(NSEEDS)
					for s in SEEDS:
						SC = SpectralClustering(n_clusters=CKi, random_state=s)
						SC.fit(X)
						Y = SC.labels_
						sscore = silhouette_score(DX, Y, metric="precomputed")
						SS.append(sscore)
						print("seeds:"+str(s))
					sbest = SEEDS[SS.index(max(SS))]
					dBS[T] = sbest
			except:
				dBS[T]=0
		dBS[KET[0]] = 0
		return dBS

	def performClustering(self):
		print('start clustering...')
		KET=self.KET
		# default clustering model
		[dCK,dBS]=self.getClusteringPars()
		#pdb.set_trace()
		AC=[]
		gc.collect()
		for i in range(len(KET)):
			print("clustering for time: "+str(KET[i]))
			ti=KET[i]
			CT = self.dET[ti]
			CKT=dCK[ti]
			BST=dBS[ti]
			
			if CKT > 1:
				if (self.largeType=='1' or self.largeType=='True'):
					X=copy.deepcopy(self.affMatrix[ti])
					SC = KMeans(n_clusters=CKT, random_state=BST)
				else:
					X=copy.deepcopy(self.affMatrix[ti])
					SC = SpectralClustering(n_clusters=CKT, random_state=BST)
				
				
				SC.fit(X)
				Y = SC.labels_
				
				for j in range(len(CT)):
					CT[j].Label = Y[j]
				CC = [Cluster([item for item in CT if item.Label == j], ti, str(ti) + '_' + str(j)) for j in range(CKT)]
				AC += CC
			else:
				for j in range(len(CT)):
					CT[j].Label = 0
				CC = [Cluster([item for item in CT if item.Label == 0], ti, str(ti)+'_'+str(0))]
				AC += CC
		return AC

# cluster
class Cluster:
	def __init__(self,cells,timepoint,ID):
		self.cells=cells           # cells (data poitns) for the cluster
		self.GL=self.cells[0].GeneList # Cell Gene list
		self.T=timepoint           # time points observed
		self.ST=None               # temp time, save the old time based on observation
		self.ID=ID                 # ID

		self.P=None                # parent cluster
		self.C=[]                  # child clusters
		self.E=self.getAvgEx()     # initial mean expression
		self.R=self.getVariance()  # initial observation variance
		self.mT=None               # initial mean Time
		self.rT=None               # intial Time variance
		self.DTA=[]                # distance to ancestor List
		self.PR=0                  # probability of cluster


	def updateCluster(self):
		self.E=self.getAvgEx()
		self.R=self.getVariance()

	def getSimilarity(self,B):
		# get similariy score between current cluster and cluster B
		#SF=1-abs(sum(self.DTA)/len(self.DTA)-sum(B.DTA)/len(B.DTA))
		SF=(1-sum(B.DTA)/len(B.DTA))/(1-sum(self.DTA)/len(self.DTA))
		SAB=getDistance(self,B)
		SAB=SF*SAB
		return SAB

	def getAvgEx(self):
		# get average experssion of current cluster
		L=len(self.cells[0].E)         # dim of cell expression
		n=len(self.cells)              # number of cells for current cluster
		AE=[]
		for i in range(L):
			iAvg=sum([item.E[i] for item in self.cells])/n
			AE.append(iAvg)
		return AE

	def getVariance(self,mu=None):
		# get the variance of genes within this clusters
		L=len(self.cells[0].E)
		X=[item.E for item in self.cells]
		R=[]
		for i in range(L):
			gi=[item[i] for item in X]
			if mu==None:
				vi=getVarianceVector(gi)
			else:
				vi=getVarianceVector(gi,mu[i])
			R.append(vi)
		pcount=0.01
		R=[item+pcount for item in R]
		return R
    #-------------------------------------------------------------------
	# get expressed TFs for given cluster
	# if expressed in at least 20% of genes => expressed
	def getExpressedTF(self,TFList):
		PTF=[]
		EXCUT=0.85		# at least 1-EXCUT non-zero => expressed
		HGL=[item.upper() for item in self.GL]
		for i in TFList:
			if i in HGL:
				ix=HGL.index(i)
				x=[item.E[ix] for item in self.cells]
				x.sort()
				if x[int(EXCUT*len(x))]>0:
					PTF.append(i)
		return PTF

	#-------------------------------------------------------------------
        #get the differential expressed TFs at given time point for each cluster
	def getDiffTF(self,TFList,AC):
		#------------------------------------------------------------------
		# check Diff TF using student t-test
		def checkDiff1(X):
			diff_cut = 0.05		# student t-test cutoff
			pv=1
			for i in range(len(X)-1):
				for j in range(i+1,len(X)):
					#pdb.set_trace()
					A=X[i]
					B=X[j]
					pv=min(pv,ttest_ind(A,B)[-1])

			return pv<=diff_cut

		#-------------------------------------------------------------------
		# check Diff TF using Fold change
		def checkDiff(X):
			LogFC=0
			FCUT=1.5                    # Fold change cutoff 1.5
			for i in range(len(X)-1):
				for j in range(i+1,len(X)):
					A=X[i]
					B=X[j]
					AA=sum(A)/len(A)
					AB=sum(B)/len(B)
					LogFC=max(LogFC,abs(AA-AB))

			return LogFC>FCUT

		#------------------------------------------------------------------
		tfs=TFList
		PT=self.T
		PL=[] # parent cluster list for given path (parent cluster means the from node of the path)
		for i in AC:
			if i.T==PT:
				PL.append(i)
		DTF={}
		HGL=[item.upper() for item in self.GL]

		for i in tfs:
			if i in HGL:
				ix=HGL.index(i)
				Ei=[]
				for j in PL:
					Eij=[item.E[ix] for item in j.cells]
					Ei.append(Eij)
				if checkDiff(Ei):
					diff=sum(Ei[0])/len(Ei[0])-sum(Ei[1])/len(Ei[1])
					diff='>' if diff>0 else '<'
					DTF[i]=diff
		return DTF

	#-------------------------------------------------------------------
	#get Pcluster, get the brother clusters and estimate transition probability
	def getPcluster(self):
		BrotherCluster=[self] if self.P==None else self.P.C
		BrotherCells=[]
		for i in BrotherCluster:
			BrotherCells+=i.cells
		return 1.0*len(self.cells)/len(BrotherCells)

	#------------------------------------------------------------------------------
	# calculate the probability of a cell belonging to given cluster using guassian distribution (expression)
	def getCellProbability(self,cell):
		# cell: a cell
		# mu,sm
		mu=self.E
		sm=self.R

		P = []
		for i in range(len(cell.E)):
			p = ((-0.5) *math.log(2*math.pi*sm[i])) - ((cell.E[i] - mu[i]) ** 2 / (2 * sm[i]))
			P.append(p)
		return P

	def getCellProbabilityT(self,cell):
		mut=self.mT[0]
		smt=self.rT[0]
		P=((-0.5) *math.log(2*math.pi*smt)) - ((cell.dta - mut) ** 2 / (2 * smt))
		#P=P*len(self.GL)
		return P

	# ------------------------------------------------------------------------
	# calculate the weight w for each gene
	def getW(self):
		# T: time point
		W = []
		BrotherCluster = [self] if self.P == None else self.P.C
		CT = reduce(lambda x, y: x + y, [item.cells for item in BrotherCluster])
		pscount=1
		for i in range(len(self.GL)):
			ei = [item.E[i] for item in CT]
			ei_nonzero = [item for item in ei if item != 0]
			wi = float(len(ei_nonzero)+pscount) / (len(ei)+pscount)
			W.append(wi)
		return W

	#---------------------------------------------------------------------------------
	# calculate the probability of a cell belonging to given cluster using mixture model (expression)
	def getCellProbabilityMixture(self,cell,W,K):
		# W: dropout weight
		X=cell.E
		p1 = self.getCellProbability(cell)
		p = [math.log(W[k] * math.exp(p1[k]) + (1 - W[k]) * K) if X[k] == 0 else math.log(W[k]) + p1[k] for k in range(len(X))]
		return p  # vector

	#---------------------------------------------------------------------------
	# calculate the total probability of a cell belonging to given cluster
	def getAssignProbability(self,cell,W,K):
		# W: dropout weight
		# K: dropout prob
		PG=sum(self.getCellProbabilityMixture(cell,W,K)) # ex likelihood
		#PG=self.getCellProbability(self,cell)
		PT=self.getCellProbabilityT(cell)             # time likelihood
		PR = math.log(self.PR)                        # probability of the cluster
		P=PG+PR+PT
		#P=PG+PR
		return P
		#-------------------------------------------------------------------

class Path:
        def __init__(self,fromNode,toNode,Nodes,dTD,dTG,dMb,fChangeCut=1):
                self.fromNode=fromNode                                                  # from Node
                self.toNode=toNode                                                      # to Node
                self.AllNodes=Nodes
                self.GL=fromNode.GL
           
				
                self.diffF=[item[0] for item in self.getDiffGene(FCUT=fChangeCut)]      # get differnetial genes based on fold change 
                self.diffT=[item[0] for item in self.getDiffGeneTTest()]                # get differential genes based on t-test
                self.diffG=[item for item in self.diffF if item in self.diffT]          # get differnetial genes based on fold change and student t-test
                self.FC=self.getFC()                                                    # fold change
		
		#pdb.set_trace()
                [self.etf,self.dtf]=self.getetf(dTD,dTG,dMb)                                       # transcription factors and diff TFs
                self.atf=self.getActiveTF(dTD,dTG,dMb)                                             # TF targets are significantly different between fromNode and toNode

		#---------------------------------------------------------------------- 
                self.B=self.getTransition(dTD,dTG,dMb,fChangeCut)                       # transition offset
                self.Q=self.getProcessVariance(MU=self.fromNode.E)                      # initial process variance

                #self.fulltext = '' # for drawing purpose
                #self.showtext = '' # for drawing purpose

        #-------------------------------------------------------------------
        # calculate fold change between fromNode (cluster) and toNode (cluster)
        def getFC(self):
                def logfc(x, y):
                        return y - x

                A=self.fromNode
                B=self.toNode

                AE=A.getAvgEx()
                BE=B.getAvgEx()
                FC=[[abs(logfc(AE[i],BE[i])),logfc(AE[i],BE[i]),i,AE[i],BE[i]] for i in range(len(AE))]

                FC.sort(reverse=True)
                FC=[[self.GL[item[2]],item[1],item[3],item[4]] for item in FC]
                return FC

        # get differential genes between clusters along the path
        def getDiffGene(self,FCUT):
                FC=self.getFC()
                #pdb.set_trace()
                DG=[item for item in FC if abs(item[1])>FCUT]
                return DG

        #-------------------------------------------------------------------
        # get differential genes between clusters along the path
        # using student t-test
        def getDiffGeneTTest(self):
                cut=5e-2
                AE=[item.E for item in self.fromNode.cells]
                BE=[item.E for item in self.toNode.cells]
                G=range(len(AE[0]))
                #pdb.set_trace()
                TT=[]
                for i in G:
                        X=[item[i] for item in AE]
                        Y=[item[i] for item in BE]
                        pxy=ttest_ind(X,Y)[-1]
                        TT.append([pxy,i])
                TT.sort()
                TT=[item for item in TT if item[0]<cut]
                DG=[[self.GL[item[1]],item[0]] for item in TT if self.GL[item[1]]]
                return DG

        #-------------------------------------------------------------------
        def getetf(self,dTD,dTG,dMb):
		# get enriched TFs based on significantly diff genes
		#---------------------------------------------------------------
		# dMi: input sequence scanning result
		# dMb: background sequence scanning result
		# n: number of sequences in input
		# dTD: dictionary TF->DNA
		# dMb: TF binding for background
		# review erniched TF
                def getEnrichTF():
                        pcut=0.1
                        dMi=batchScanPrior([item.upper() for item in self.diffG],dTD)
                        K=[item for item in dMi.keys() if item in dMb.keys()]
                        K.sort()
                        n=len(self.diffG)		# number of diff genes
                        N=len(self.GL)				# N: number of sequences in background (all)
                        entf=[]
                        #pdb.set_trace()
                        for i in K:
                                Ti=len(dMi[i])
                                Tb=len(dMb[i])
                                pr=float(Tb)/N
                                pvi=1-pbinom(Ti-1,n,pr)
                                if pvi<pcut:
                                        entf.append([pvi,i])
                    
			#pdb.set_trace()                     
                        entf.sort()
                        return entf
                #-------------------------------------------------------------
                etf=getEnrichTF()
                ptf=self.fromNode.getExpressedTF(dTD.keys())
                dtf=self.fromNode.getDiffTF(dTD.keys(),self.AllNodes)
                etf=[item for item in etf if item[1] in ptf]
                return [etf,dtf]

        #-------------------------------------------------------------------
        def getActiveTF(self,dTD,dTG,dMb):
                # the idea : the targets for given TFs are signiciantly differenent between the fromNode and toNode
                HGL=[item.upper() for item in self.GL]
                pv_cut=1e-1
                ATF=[]

                FE=self.fromNode.getAvgEx()
                TE=self.toNode.getAvgEx()
                for i in dMb:
                        tfi=dMb[i]
                        tfi=[HGL.index(item) for item in tfi]

                        fe=[FE[item] for item in tfi]
                        te=[TE[item] for item in tfi]

                        #fe=[self.fromNode.E[item] for item in tfi]
                        #te=[self.toNode.E[item] for item in tfi]
                        #fe=reduce(lambda x,y:x+y, [[item.E[j] for item in self.fromNode.cells] for j in tfi])
                        #te=reduce(lambda x,y:x+y, [[item.E[j] for item in self.toNode.cells] for j in tfi])
                        #pdb.set_trace()
                        pv=wilcoxon(fe,te)[-1]
                        if pv<pv_cut:
                                ATF.append([pv,i])
                ATF.sort()
                return [item[1] for item in ATF]

        #-------------------------------------------------------------------
        # regresion model for each path
        def getTransition(self,dTD,dTG,dMb,FCUT=1):
                G = self.getFC()
                dFC = {item[0].upper(): item[1] for item in G}
                etfID = [item[1] for item in self.etf]
                [X, Y,U,D] = buildTrain(G, dTG, etfID,self.GL,FCUT)
                LR = LogisticRegressionCV(penalty='l1', Cs=[1.5, 2, 3, 4, 5], solver='liblinear', multi_class='ovr')
                dR = {0: U, 1: D, 2: 0}
                HGL = [item.upper() for item in self.GL]
                try:
                        LR.fit(X, Y)
                        CE = LR.coef_
                        petf = parseLR(self.etf, CE)
                        # ---------------------------------------------------------
                        XX = []
                        for i in HGL:
                                if i in dTG:
                                        tfi = dTG[i]
                                        xi = [1 if item in tfi else 0 for item in etfID]
                                else:
                                        xi = [0] * len(etfID)
                                XX.append(xi)
                        YY = LR.predict(XX)
                        self.etf = petf
                except:
                        YY = [0 if dFC[item] > FCUT else 1 if dFC[item] < -1 * FCUT else 2 for item in HGL]
                YY = [dR[item] for item in YY]
                return YY

        #-------------------------------------------------------------------
        def getProcessVariance(self,MU=None):
                # MU : average at time t-1, vector
                # X1: all observation at time point t
                # X2:  all observations at time point t
                X1=[item.E for item in self.fromNode.cells]
                X2=[item.E for item in self.toNode.cells]
                dim=len(self.GL)
                Q=[]
                if MU==None:
                        for i in range(dim):
                                x1=[item[i] for item in X1] # all observation for gene i at time t-1
                                mui=sum(x1)/len(x1)+self.B[i]
                                x2=[item[i] for item in X2] # all observation for gene i at time t

                                v=getVarianceVector(x2,mui)
                                Q.append(v)
                else:
                        for i in range(dim):
                                x2=[item[i] for item in X2]
                                mui=MU[i]+self.B[i]
                                v=getVarianceVector(x2,mui)
                                Q.append(v)
                pcount=0.01
                Q=[item+pcount for item in Q]
                return Q

class Graph:
	def __init__(self,Cells,tfdna,kc,largeType=None,dsync=None,virtualAncestor=None,fChangeCut=1,etfile=None):
		self.Cells=Cells
		self.largeType=largeType
		self.dsync=dsync
		self.virtualAncestor=virtualAncestor
		self.fChangeCut=fChangeCut
		self.etfile=etfile
		
		
		#Current mode requires an ancestor node (only 1 cluster). If this is not case, virtual ancestor node is needed (The mean expression of the first time point)
		if (virtualAncestor=='True' or virtualAncestor=='1'):
			va=buildVirtualAncestor(self.Cells)
			self.Cells.insert(0,va)
	
		self.clustering=Clustering(self.Cells,kc,largeType)#
		self.Nodes=self.buildNodes()							# build nodees based on cells
		self.GL=self.Nodes[0].GL 
		
		#----------------------------------------------------------------------------------------------
		[self.dTD,self.dTG,self.dMb]=self.parseTFDNA(tfdna)							# Gene List									
		self.updateTime()
		print("building graph...")
		if (self.dsync==None):
			self.splitTime()							# update Nodes with time splitting
			
		self.W=self.getW()								# get the dropout rate wight W in the single cell experiment
		#---------------------------------------------------------------------------------

		self.connectNodes()  # find connecting relationship between nodes
		self.getNodePR()     # calculate the probability of Node
		#pdb.set_trace()
		self.Edges = self.buildEdges()  # build edges based on nodes
		self.Paths = self.buildPaths()  # build paths based on edges
		
		# get representing TFs for each of the clusters
		self.adjustRTFs(self.fChangeCut,self.etfile)
		
		#pdb.set_trace()
		[self.Q, self.R] = self.estimateQR()  # estiamte Q,R for Kalman Filter
		self.rEstimateEx()		# estimate the expression for nodes (clusters) using Kalman Filter
		self.rEstimateT()		# estimate time for each node (clustering) using Kalman Filter
		
      
      #---------------------------------------------------------------------------------------------------
	def adjustRTFs(self,fcut=1,tflist=None):
		# get RTFs (representating TFs) based on its own expression for each node (the TF expression is different to both parent and siblings)
		tflistpath=pkg_resources.resource_filename(__name__,"tfdata/HumanTFList.txt") if tflist==None else tflist
		
		try:
			with open(tflistpath,'r') as f:
				TFs=f.readlines()
				TFs=[item.strip().split()[0] for item in TFs]	
		except:
			print("error! Please check your input TF List file")
			sys.exit(0)
					
		eTFs=[item for item in [item.upper() for item in self.GL] if item in TFs]
		
		for Node in self.Nodes:
			if Node.P:
				NodeParent=Node.P
				NodeSib=[item for item in Node.P.C if item!=Node]
				NodeSibCells=[] if NodeSib==[] else functools.reduce(lambda x,y:x+y,[item.cells for item in NodeSib])
				NodeParentCells=NodeParent.cells
				peTFs=[]
				for j in eTFs:
					jdex=[item.upper() for item in self.GL].index(j)
					[jflag1,pvp,fcp]=tellDifference(Node.cells,NodeParentCells,jdex,fcut)
					[jflag2,pvs,fcs]=tellDifference(Node.cells,NodeSibCells,jdex,fcut)
					if (jflag1*jflag2>0) or (jflag1!=0 and len(NodeSibCells)==0):
						peTFs.append([pvp,j,fcp,fcs])
				peTFs.sort()
				Node.eTF=peTFs
			else:
				Node.eTF=[]
		   
        #----------------------------------------------------------------------------------
	
	def parseTFDNA(self,tfdna):
		
		RTD=TabFile(tfdna).read('\t') # dictionary to store the TF-DNA info
		try:
			TD=[[item[0].upper(),item[1].upper()] for item in RTD]
		except:
			print("check the format of input TF-DNA interaction file")
			sys.exit(0)
			
		[dTD,dTG]=getTFDNAInteraction(TD)
		# TF binding in all input sequences (background)
		dMb=batchScanPrior([item.upper() for item in self.GL],dTD)
		return [dTD,dTG,dMb]
		
	# calculate the weight w for each gene
	def getW(self):
		# T: time point
		pscount=1
		W=[]
		for i in range(len(self.GL)):
			ei = [item.E[i] for item in self.Cells]
			ei_nonzero = [item for item in ei if item != 0]
			wi = float(len(ei_nonzero)+pscount) / (len(ei)+pscount)
			W.append(wi)
		return W

	def updateGraph(self):
		# nodes already there, no need to re-build
		self.connectNodes()                         # re-connect
		self.getNodePR()
		self.Edges=self.buildEdges()
		self.Paths=self.buildPaths()
		self.adjustRTFs(self.fChangeCut,self.etfile)
		[self.Q, self.R] = self.estimateQR()
		self.rEstimateEx()
		self.rEstimateT()

	def buildNodes(self):
		return(self.clustering.performClustering())
		
	def updateTime(self):
		# calculate dta for each cell
		Ancestor = self.Nodes[0]
		for i in self.Cells:
			i.dta=i.distanceTo(Ancestor)
		KET=self.clustering.KET
		# calculate level for each cluster
		for i in self.Nodes:
			ti=KET.index(i.T)
			i.ST=i.T
			i.T=ti
		#pdb.set_trace()
			
	def splitTime(self):
		# -----------------------------------------
		print("time adjustment...")
		def checkNeighbors1(XX):
			cut = 0.1
			BreakP = []
			for i in range(len(XX) - 1):
				pvr = ranksums(XX[i].DTA, XX[i + 1].DTA)[-1]
				if pvr < cut:
					BreakP.append(i + 1)

			CC = []
			st = 0
			for i in BreakP:
				end = i
				CC.append(XX[st:end])
				st = end
			CC.append(XX[st:])
			return CC

		def checkNeighbors(XX,ANS):
			pv=[]
			df=[]
			ANSDTA=sum(ANS.DTA)/len(ANS.DTA)
			pvcut=0.05  #p-value cutoff 0.5
			dfcut=0.1  # 10% difference of STA average (normalized)
			ASTA=[sum(item.DTA)/len(item.DTA) for item in XX]

			ASTA=[abs(ANSDTA-item)/(ANSDTA-min(ASTA)) for item in ASTA] # normalization
			for i in range(len(XX)-1):
				pvi=ranksums(XX[i].DTA,XX[i+1].DTA)[-1]
				pv.append(pvi)
				dfi=ASTA[i+1]-ASTA[i]
				df.append(dfi)
			rpv=[sorted(pv).index(item) for item in pv]
			rdf=[sorted(df,reverse=True).index(item) for item in df]
			tr=[rpv[i]+rdf[i] for i in range(len(rpv))]
			tr=[[tr[i],i] for i in range(len(tr))]
			tr.sort()
			KET=self.clustering.KET
			[NCL,NCH]=[len(KET)-2,(len(KET)-1)*2-1]                     # candidates for number of cuts
			
			NCL=max(0,len(KET)-1)
			NCH=int(math.log(len(XX)+1,2))                               # max # of level (determined by the maximal # of binary splits)
			
			breakP=[tr[k][1] for k in range(0,NCL)]
			for i in range(NCL,NCH):
				if pv[tr[i][1]]<=pvcut or df[tr[i][1]]>=dfcut:
					breakP.append(tr[i][1])

			if len(breakP)==NCL:
				dfcut=dfcut*0.75
			breakP=[tr[k][1] for k in range(0,NCL)]
			for i in range(NCL,NCH):
				if pv[tr[i][1]]<=pvcut or df[tr[i][1]]>=dfcut:
					breakP.append(tr[i][1])

			i=NCH

			#-----------------------------------------------
			try:
				if pv[tr[i][1]] <= pvcut or df[tr[i][1]] >= dfcut:
					if abs(tr[i-1][1]-tr[i][1])<=1 and abs(tr[i-1][0]-tr[i][0])<=1:
						if tr[i-1][1] in breakP:
							breakP.remove(tr[i-1][1])
							breakP.append(tr[i][1])
			except:
				pass

			breakP.sort()
			#pdb.set_trace()
			CC = []
			st = 0
			for i in breakP:
				end = i+1
				CC.append(XX[st:end])
				st = end
			CC.append(XX[st:])
			return CC
			
		# -------------------------------------------------------------
		# splitTime main 
		"""
		Ancestor = self.Nodes[0]
		for i in self.Cells:
			i.dta=i.distanceTo(Ancestor)
		"""
		
		for i in self.Nodes:
			i.DTA = [item.dta for item in i.cells]    # distance to ancestor vector
		XX = sorted(self.Nodes[1:], key=lambda item: sum(item.DTA)/len(item.DTA), reverse=True)
		CC = checkNeighbors(XX,self.Nodes[0])

		for i in range(len(CC)):
			for j in CC[i]:
				j.T = i + 1

		self.Nodes[0].ST =self.Nodes[0].T
		self.Nodes[0].T = 0

		self.Nodes[1:] = XX

		return self.Nodes


	# get Parent cluster for each node 
	def getParent(self,A,TP):
		#A is the cluster,we are trying to find the parent for cluster A.
		#TP is the parent level
		#PL is the cluste list in parent level
	
		PL=[item for item in self.Nodes if (item.T==TP and item.ST<=A.ST)]
		PL=sorted(PL,key=lambda item:sum([getDistance(item,ck) for ck in A.cells])/len(A.cells))
		#pdb.set_trace()
		#---------------------------------------------------------------
		# if time sync is disabled, output the closest node in parent level
		
		if self.dsync=='True' or self.dsync=='1':
			if len(PL)>0:
				return PL[0]
			return None
			
		#--------------------------------------------------------------
		
		pvcut=0.1
		if len(PL)>1:
			X=[getDistance(PL[0],item) for item in A.cells]
			Y=[getDistance(PL[-1],item)  for item in A.cells]
			# Length adjustment if Vector is too short
			SizeFactor=50 # Used for length adjustment
			X=X*int(SizeFactor/len(X)) if len(X)<SizeFactor else X
			Y=Y*int(SizeFactor/len(Y)) if len(Y)<SizeFactor else Y
			pv=ranksums(X,Y)[-1]
			if pv<pvcut:
				return PL[0]
			else:
			   return self.getParent(A,TP-1)
		elif len(PL)==1:
			if PL[0]==self.Nodes[0]:
				return PL[0]
			else:
				PLL=[item for item in self.Nodes if (item.T==TP-1 and item.ST<=A.ST)]
				if len(PLL)==1:
					return PL[0]
				else:
					return self.getParent(A,TP-1)
		else:
			if A!=self.Nodes[0]:
				return self.getParent(A,TP-1)

	def connectNodes(self):
		print("connecting nodes ....")
		# get Previous level for parent candidate
		KET = sorted(list(set([item.T for item in self.Nodes])))

		# delete old parent-child relationship
		for i in self.Nodes:
			i.P=None
			i.C=[]

		for i in self.Nodes:
			i.P = self.getParent(i,i.T-1)
			if i.P:
				i.P.C += [i]
		#pdb.set_trace()


	def buildEdges(self):
		P = []
		for i in self.Nodes:
			if i.P:
				p1 = Path(i.P, i,self.Nodes,self.dTD,self.dTG,self.dMb,self.fChangeCut)
				P.append(p1)
		return P


	def buildPaths(self):
		def getCompletePath(en):
			# en: end node
			for i in self.Edges:
				if i.toNode == en:
					return getCompletePath(i.fromNode) + [i]
			return []

		CP = []  # complete path
		for i in self.Nodes:
			if not i.C:
				cp =getCompletePath(i)
				if cp!=[]:
					CP.append(cp)
		#pdb.set_trace()
		return CP

	def estimateQR(self):
		# -----------------------------------------------------------------------
		# estimate Q and R based on current assignment
		# AC: clusters
		# P: Paths
		AC=self.Nodes
		P=self.Edges

		Q = []
		R = []
		for i in AC:
			Ri = i.R
			R.append(Ri)

		for i in P:
			Qi = i.Q
			Q.append(Qi)

		Q = [sum([item[i] for item in Q]) / len(Q) for i in range(len(self.GL))]
		R = [sum([item[i] for item in R]) / len(R) for i in range(len(self.GL))]
		return Q, R


	def estimateQRT(self):
		#-----------------------------------------------------------------------
		# estimate Qt and Rt based on current assignemnt

		# ------------------------------------------------------------------
		for i in self.Nodes:
			i.mT= [sum(i.DTA)/len(i.DTA)]
			i.rT =[getVarianceVector(i.DTA)]  # get initial time variance for each time
		Rt=[item.rT[0] for item in self.Nodes]
		Rt=[max(Rt)]

		KET = sorted(list(set([item.T for item in self.Nodes])))
		Tstart = sum([item.dta for item in self.Nodes[0].cells]) / len(self.Nodes[0].cells)
		Tend = sum([item.dta for item in self.Nodes[-1].cells]) / len(self.Nodes[-1].cells)
		Tstep = (Tstart - Tend) / len(KET)

		for i in self.Edges:
			i.Qt=getVarianceVector([item.dta for item in i.toNode.cells],i.fromNode.mT[0]+Tstep)

		Qt=[item.Qt for item in self.Edges]
		Qt=[max(Qt)]
		return [Qt,Rt]


	def rEstimateEx(self):
		[Q, R] = self.estimateQR()
		mu0 = self.Nodes[0].E
		sm0 = self.Nodes[0].R
		pcount = 0.01  # pseudocount, handling zero variance

		for cp in self.Paths:
			an = []
			an.append(cp[0].fromNode)
			for j in cp:
				an.append(j.toNode)
			X = [[item.E for item in j.cells] for j in an]
			# -----------------------------------------------
			B = [0] * len(self.GL)
			B = [B] + [item.B for item in cp]

			A = [[1] * len(self.GL)] * len(X)
			H = [[1] * len(self.GL)] * len(X)
			I = [[0] * len(self.GL)] * len(X)

			"""
			FF=KalmanFilterInd(A,B,H,I,mu0,sm0,Q,R)
			[mu,sm]=FF.filterMultiple(X)
			 """
			[mu, sm] = smootherMultiple(A, B, H, I, mu0, sm0, Q, R, X)
			for j in range(len(an)):
				an[j].E = mu[j]
				an[j].R = [item + pcount for item in sm[j]]

	def rEstimateT(self):
		#todo re-write time estimte using kalman filter
		KET = sorted(list(set([item.T for item in self.Nodes])))
		pcount = 0.01  # pseudocount, handling zero variance

		#[Qt, Rt] = self.estimateQRT()
		#mut0=self.Nodes[0].mT
		#smt0=self.Nodes[0].rT
		Qt = [1]
		Rt = [0.5]
		#Rt = [0.5]

		Tstart = sum([item.dta for item in self.Nodes[0].cells]) / len(self.Nodes[0].cells)
		mut0 = [Tstart]
		smt0 = [0.1]
	
		#--------------------------------------------------------------------
		Tstart = sum([item.dta for item in self.Nodes[0].cells]) / len(self.Nodes[0].cells)
		Tend = sum([item.dta for item in self.Nodes[-1].cells]) / len(self.Nodes[-1].cells)
		Tstep = (Tstart - Tend) / len(KET)
		for cp in self.Paths:
			an = []
			an.append(cp[0].fromNode)
			for j in cp:
				an.append(j.toNode)
			TX = [[[item.dta] for item in j.cells] for j in an]
			A = [[1]] * len(KET)
			B = [[0]] + [[-1 * Tstep]] * (len(KET) - 1)
			H = [[1]] * len(KET)
			I = [[0]] * len(KET)
			[mut, smt] = smootherMultiple(A, B, H, I, mut0, smt0, Qt, Rt, TX)
			for j in range(len(an)):
				an[j].mT = mut[j]
				an[j].rT = [item + pcount for item in smt[j]]
			#pdb.set_trace()
		
	def rEstimateT1(self):
		#todo modify time estimate using Kalman Filter
		Qt = [1]
		Rt = [0.1]
		smt0 = [1]
		KET = sorted(list(set([item.T for item in self.Nodes])))
		Tstart=sum([item.dta for item in self.Nodes[0].cells])/len(self.Nodes[0].cells)
		Tend=sum([item.dta for item in self.Nodes[-1].cells])/len(self.Nodes[-1].cells)
		Tstep=(Tstart-Tend)/len(KET)
		#pdb.set_trace()
		#-------------------------------------
		A = [[1]] * len(KET)
		B = [[0]] + [[-1*Tstep]] * (len(KET) - 1)
		H = [[1]] * len(KET)
		I = [[0]] * len(KET)
		#------------------------------------
		TX=[]
		for i in KET:
			ci=[item for item in self.Nodes if item.T==i]
			TXI=reduce(lambda x,y: x+y, [[[item.dta] for item in j.cells] for j in ci])
			TX.append(TXI)
		#pdb.set_trace()
		[mut, smt] = smootherMultiple(A, B, H, I, mut0, smt0, Qt, Rt, TX)

		for i in self.Nodes:
			ivt=KET.index(i.T)
			i.mT=mut[ivt]
			i.rT=smt[ivt]
		
		
	def getNodePR(self):
		# -------------------------------------------------------------------------
		# calculate the probability for each cluster
		for i in self.Nodes:
			i.PR = i.getPcluster()

		for i in self.Nodes:
			t = i
			pr = 1
			while (t):
				pr = pr * t.PR
				t = t.P
			i.PR = pr

	def ReAssign(self):
		#-----------------------------------------
		# re-assign
		AS = []
		Tlli = 0
		#K=0.05 #  mixture probability constant.
		K=0.01 #  mixture probability constant.
		for i in self.Cells:
			pi=[j.getAssignProbability(i,self.W,K) for j in self.Nodes]
			Tlli += max(pi)
			AS.append(pi.index(max(pi)))
		RAC = [[self.Cells[i] for i in range(len(self.Cells)) if AS[i] == j] for j in range(len(self.Nodes))]
		#pdb.set_trace()
		#-----------------------------------------
		# update cluster cells --nodes
		for i in range(len(self.Nodes)):
			self.Nodes[i].cells = RAC[i]

		# remove empty cluster
		self.Nodes = [item for item in self.Nodes if item.cells != []]

		# update cluster E,R

		for i in range(len(self.Nodes)):
			self.Nodes[i].E=self.Nodes[i].getAvgEx()
			self.Nodes[i].R=self.Nodes[i].getVariance()

#=======================================================================
# global functions


def tellDifference(nodeCells,nodePCells,geneIndex,fcut=0.6):
	if len(nodePCells)==0:
		return [0,1,0]
	X=[item.E[geneIndex] for item in nodeCells]
	Y=[item.E[geneIndex] for item in nodePCells]
	fc=sum(X)/len(X)-sum(Y)/len(Y)
	pv=ttest_ind(X,Y)[1]
	
	pcut=0.05
	if pv<pcut and fc>fcut:
		return [1,pv,fc]
	if pv<pcut and fc<-1*fcut:
		return [-1,pv,fc]	
	return [0,pv,fc]
#-----------------------------------------------------------------------------
# filter X
def filterG(X,GL):
	#X :cells
	Y=[]
	for i in X:
			j=[i.E[k] for k in GL]
			Y.append(j)
	return Y

# calculate distance between 2 clusters
def getDistance(x,y,dtype='spearman',metric=None):
	#-----------------------------------------------------------------------
	# eucliden distance cluster-to-cluster
	# x: cluster, y:cluster
	def euD(x,y):
			s=0
			for i in range(len(x.E)):
					di=(x.E[i]-y.E[i])**2
					s+=di
			s=math.sqrt(s)
			return s

	# spearmanr distance
	def spD(x,y):
			return 1-spearmanr(x.E,y.E)[0]

	# pearsonr distance

	def peD(x,y):
			return 1-pearsonr(x.E,y.E)[0]

	# manhattan distance
	def maD(x,y):
			s=0
			for i in range(len(x.E)):
					di=abs(x.E[i]-y.E[i])
					s+=di
			return s

	#distance cell-to-cell
	def DISTC(x,y,DF):
			# x: cluster x
			# y: cluster y
			# DF: distance function

			S = []
			for i in x.cells:
					for j in y.cells:
							dij = DF(i, j)
							S.append(dij)
			SAB = sum(S) / len(S)
			return SAB

	if dtype=='spearman':
			DF=spD
	elif dtype=='manhattan':
			DF=maD
	elif dtype=='pearson':
			DF=peD
	else:
			DF=euD

	if metric=='cell-to-cell':
			disxy=DISTC(x,y,DF)
	else:
			disxy=DF(x,y)

	return disxy
#-----------------------------------------------------------------------
# Peak detection
def detPeak(K, A, deltaP, type='max'):
	# get Local optimal score
	LM = []
	LX = []
	if len(A)==0:
		return []
	minA = min(0, min(A))
	for i in range(len(A) - 1):
			delta = deltaP * (A[i] - minA)
			if i == 0:
					ad = A[i + 1] - A[i]  # after difference
					if ad > delta:
							LM.append(i)
					if ad < -1 * delta:
							LX.append(i)
			else:
					pd = A[i] - A[i - 1]
					ad = A[i + 1] - A[i]
					if (pd > delta) and (ad < -1 * delta):
							# pdb.set_trace()
							LX.append(i)
					if (pd < -1 * delta) and (ad > delta):
							LM.append(i)
	if type == 'min':
			LM = sorted(LM, key=lambda x: A[x])
			return [K[item] for item in LM]

	LX = sorted(LX, key=lambda x: A[x], reverse=True)
	return [K[item] for item in LX]
#----------------------------------------------------------------------

def getTFDNAInteraction(TD):
	dTD = {}  # TF->DNA
	dTG = {}  # DNA->TF
	for i in TD:
			if i[0] not in dTD:
					dTD[i[0]] = [i[1]]
			else:
					dTD[i[0]].append(i[1])
					
			if i[1] not in dTG:
					dTG[i[1]] = [i[0]]
			else:
					dTG[i[1]].append(i[0])
	return [dTD,dTG]

# scan motifs-----------------------------------------------------------
# log- motifs

def logMotif(M,cut):
	for i in range(len(M)):
			MID=M[i][0]
			M[i]=M[i][1:]
			for j in range(len(M[i])):
					M[i][j]=[float(item) for item in M[i][j]]
					M[i][j]=[item+0.01 for item in M[i][j]]
					M[i][j]=[item/sum(M[i][j]) for item in M[i][j]]
					M[i][j]=[math.log(item/0.25) for item in M[i][j]]
			#----------------------------------------------------------
			MX=0
			for j in range(len(M[i])):
					maxI=max(M[i][j])
					MX+=maxI
			SC=MX*cut
			M[i]=[MID+'\t'+str(SC)]+M[i]

	return M

#-----------------------------------------------------------------------
# motif scanning...
def batchScan(SS,MS,cut):
	# motif scanning...
	def scan(S,M):
			def matchscore(x,y):
					dN={'A':0,'C':1,'G':2,'T':3}
					return sum([0 if x[i] not in dN else y[i][dN[x[i]]] for i in range(len(x))])

			[MID,SC]=M[0].split('\t')
			SC=float(SC)
			M=M[1:]
			SID=S[0]
			S=S[1].upper()
			for i in range(len(S)-len(M)+1):
					si=matchscore(S[i:i+len(M)],M)
					if si>SC:
							return True
			return False

	#----------------------------
	dM={}
	ct=0
	for i in MS:
			iid=i[0].split()[0]
			for j in SS:
					jid=j[0]
					if scan(j,i):
							if iid not in dM:
									dM[iid]=[jid]
							else:
									dM[iid].append(jid)
			ct+=1
			print(ct)
	return dM

#-----------------------------------------------------------------------
# scanning TF-DNA interaction prior
def batchScanPrior(A,dTD):
	# dTD  -> dictionary of TF-DNA interaction
	# A -> Gene list
	K=list(dTD.keys())
	K.sort()
	dM={}
	dA={item:0 for item in A}
	for i in K:
			GI=dTD[i]
			GI=list(set([item for item in GI if item in dA]))
			if len(GI)>0:
					dM[i]=GI

	return dM

#-----------------------------------------------------------------------
# building traning dataset for regression
def buildTrain(G,dTG,etf,GL,Fcut=1):
	# G: differential genes for a given path
	# dTD: DNA->TF dictionary
	# TF candidate
	Ncut=Fcut/2.0
	
	#UP=[item[0].upper() for item in G if item[1]>Fcut]
	#DN=[item[0].upper() for item in G if item[1]<-1*Fcut]
	#NN=[item[0].upper() for item in G if abs(item[1])<Ncut]
	UP=[item for item in G if item[1]>Fcut]
	DN=[item for item in G if item[1]<-1*Fcut]
	NN=[item for item in G if abs(item[1])<Ncut]
	
	
	U=sum([item[1] for item in UP])/len(UP)
	D=sum([item[1] for item in DN])/len(DN)
	
	UP=[item[0].upper() for item in UP]
	DN=[item[0].upper() for item in DN]
	NN=[item[0].upper() for item in NN]
	
	
	XU=[]
	XD=[]
	XN=[]

	YU=[]
	YD=[]
	YN=[]

	HGL=[item.upper() for item in GL]
	for i in HGL:
			if i in dTG:
					tfi=dTG[i]
					xi=[1 if item in tfi else 0 for item in etf]
					if i in UP:
							yi=0
							XU.append(xi)
							YU.append(yi)
					elif i in DN:
							yi=1
							XD.append(xi)
							YD.append(yi)
					elif i in NN:
							yi=2
							XN.append(xi)
							YN.append(yi)

	X=XU+XD+XN
	Y=YU+YD+YN
	
	# to solve the imbalanced training set issue, use over-sampling techqniue- SMOTE
	sm=SMOTE(random_state=0)
	Xs,Ys=sm.fit_sample(X,Y)
	
	Xs=list(Xs)
	Ys=list(Ys)
	
	#pdb.set_trace()
	return [Xs,Ys,U,D]

# parse Logistic regression result
def parseLR(etf,LRC):
	out_etf=[]
	for i in range(len(LRC[0])):
			ct=0
			for j in LRC:
					ct+=1 if j[i]==0 else 0
			if ct!=len(LRC):
					out_etf.append(etf[i])
	return out_etf

#------------------------------------------------------------------------

#-----------------------------------------------------------------------
# calculate variance for given vector
def getVarianceVector(x,mu=None):
	# x:the observation of one gene at specific time t
	if mu==None:
			mu=float(sum(x))/len(x)
	v=[(item-mu)**2 for item in x]
	v=sum(v)*1.0/len(x)
	return v
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# build the virtual ancestor node
def buildVirtualAncestor(AllCells,VT=None):
	TT=[item.T for item in AllCells]
	FirstT=sorted(TT)[0]
	FirstTCells=[item for item in AllCells if item.T==FirstT]
	GL=AllCells[0].GeneList
	L=len(FirstTCells[0].E)         # dim of cell expression
	n=len(FirstTCells)              # number of cells for current cluster
	AE=[]
	for i in range(L):
		iAvg=sum([item.E[i] for item in FirstTCells])/n
		AE.append(iAvg)
	TA=float(VT) if VT!=None else FirstT-1
	virtual_ancestor=Cell('virtual_ancestor',TA,AE,'NA',GL)
	return virtual_ancestor

#-----------------------------------------------------------------------
# differenc between ACL and ACL_update

def ClusteringDifference(ACL,ACL_update):
	
	totaloverlap=0
	for i in ACL_update:
		ci=[]
		for j in ACL:
			ov=[item for item in i if item in j]
			ci.append(len(ov))
		maxci=max(ci)
		totaloverlap+=maxci
	pertotaloverlap=totaloverlap*1.0/sum([len(item) for item in ACL_update])
	pdiff=1-pertotaloverlap
	return pdiff
	
#----------------------------------------------------------------------
#=======================================================================
# MAIN program starts here!
def  main():
	parser=argparse.ArgumentParser()
	parser.add_argument('-i','--input',required=True,help='input single cell RNA-seq expression data')
	parser.add_argument('-t','--tf_dna',required=True,help='TF-DNA interactions used in the analysis')
	parser.add_argument('-k','--clusters',required=True,default='auto', help='how to learn the number of clusters for each time point? user-defined or auto?  if user-defined, please specify the configuration file')
	parser.add_argument('-o','--output',required=True, help='output folder to store all results')
	parser.add_argument('-l','--large',required=False, help='(1/None), Optional, specify whether the data is largeType, in which case PCA+KMeans will be used for clustering instead of spectral clustering ' + 
                                                                'to improve the time and space efficiency. It is recommended for large dataset with more than 2k cells.')
	parser.add_argument('-s','--speedup',required=False, help='(1/None), Optional, if set as 1, scdiff will speedup the running by reducing the iteration times.')
	parser.add_argument('-d','--dsync',required=False,help='(1/None), Optional, if set as 1, the cell synchronization will be disabled. The cell capture time will be used directly. ' +
                                                               'This option is recommended when the users believe that the cells captured at the same time are mostly at similar differentiation stage.')
	parser.add_argument('-a','--virtualAncestor',required=False,help='(1/None), Optional, By default, scidff uses the first time point as the ancestor for all following time points. ' +
                                                                         'It is recommended to use this option if users believe that the cells at the first time points are already well differentiated and there ' +
                                                                         'exits at least 2 clusters/sub-types at the first time point. To enable this option, set it as 1.')
	parser.add_argument("-f",'--log2foldchangecut',required=False, default=1, help='Float, Optional, by default, scdiff uses log2 Fold change 1(~2^1=2) as the cutoff for differential genes (together with t-test p-value cutoff 0.05). '+ 
											   'However, users are allowed to customize the cutoff based on their application scenario (e.g. log2 fold change 1.5).')
	
	parser.add_argument("-e",'--etfListFile',required=False,help='String, Optional, by default, scdiff recognizes 1.6k TFs (we collected in human and mouse).  Users are able to provide a customized list of TFs instead using this option. '+
												'It specifies the path to the TF list file, in which each line is a TF name.')
	
	args=parser.parse_args()

	scg=args.input
	kc=args.clusters
	output=args.output
	tfdna=args.tf_dna
	largeType=args.large
	dsync=args.dsync
	virtualAncestor=args.virtualAncestor
	etf=args.etfListFile
	
	#pdb.set_trace()
	
	try:
		fChangeCut=float(args.log2foldchangecut)
	except:
		print("Error! log2foldchangecut (-f) must be a float, please check your input.")
		sys.exit(0)
		
	#pdb.set_trace()
	#-----------------------------------------------------------------------
	# 1) : read in gene expression
	AllCells=[]
	print("reading cells...")
	with open(scg,'r') as f:
		line_ct=0
		for line in f:
			if (len(line.strip())==0):
				continue
				
			if line_ct==0:
				GL=line.strip().split("\t")[3:]
			else:
				line=line.strip().split("\t")
				iid=line[0]
				ti=float(line[1])
				li=line[2]
				ei=[round(float(item),2) for item in line[3:]]
				ci=Cell(iid,ti,ei,li,GL)
				AllCells.append(ci)
			line_ct+=1
			print('cell:'+str(line_ct))
			
	firstTime=min([float(item.T) for item in AllCells])
	
	#===================================================================
	# error processing ...
	if kc!="auto":
		kcLines = TabFile(kc).read("\t")
		firstConfigTime=float(kcLines[0][0])
		firstConfigCluster=float(kcLines[0][1])
		
		if firstConfigCluster>1:
			print("Error! the first line of the config file must have only 1 cluster. Please add another virtual ancestor time point as the first line if you have multiple clusters at the first time point in your config file. Virtual ancestor option also needs to be enabled (-a 1) ")
			sys.exit(0)	
		if (firstConfigTime<firstTime) and (virtualAncestor==None):
			print("Error! VirtualAncestor option deteced, but VirtualAncestor option (-a) is not enabled. Please enable the Virtual Ancestor option.")
			sys.exit(0)
		elif (firstConfigTime==firstTime) and (virtualAncestor!=None):
			print("Error! Virtual Ancestor option enabled, but no Virtual Ancestor time point information was given in the config file. Please add Virtual Ancestor time point in the config file or disable the opition.")
			sys.exit(0)
	#=======================================================================
	# 3): Clustering starts here!
	G1=Graph(AllCells,tfdna,kc,largeType,dsync,virtualAncestor,fChangeCut,etf)

	#========================================================================
	#drawing graphs
	if os.path.exists(output)==False:
		os.mkdir(output)

	scg_name=scg.split('/')[-1]
	viz(scg_name,G1,output)
	
	sflag=0
	if (args.speedup=='1') or (args.speedup=='True'):
		sflag=1
		pdiff=0.05 
		
	#=======================================================================
	# start cell-reassignment
	# starting Kalman Filter --Expression /Time
	#-----------------------------------------------------------------------
	#Initial parameters for Kalman Filter -Time
	condition=True
	maxLoop=10
	lct=0
	while condition:
		ACL=[sorted([item.ID for item in K.cells]) for K in G1.Nodes]
		G1.ReAssign()
		G1.updateGraph()
		ACL_update=[sorted([item.ID for item in K.cells]) for K in G1.Nodes]
		if sflag!=1:
			condition = ((ACL != ACL_update) and (lct < maxLoop))
		else:	
			condition=((ClusteringDifference(ACL,ACL_update)>pdiff) and (lct<maxLoop))
		lct+=1
		viz(scg_name,G1,output)
	print("done!")
	#======================================================================#

if __name__=='__main__':
	main()
