#/usr/bin/env python

#-----------------------------------------------------------------------
# author: 


import random
import numpy as np


#=======================================================================
# Kalman Filter implementation

# Implements a linear Kalman filter.
class KalmanFilterLinear:
	def __init__(self,A, B, H, x, P, Q, R):
		self.A = A                      # State transition matrix.
		self.B = B                      # Control matrix.
		self.H = H                      # Observation matrix.
		self.current_state_estimate = x # Initial state estimate.
		self.current_prob_estimate = P  # Initial covariance estimate.
		self.Q = Q                      # Estimated error in process.
		self.R = R                      # Estimated error in measurements.

	def GetCurrentState(self):
		return self.current_state_estimate

	def Step(self,control_vector,measurement_vector):
		#---------------------------Prediction step-----------------------------
		predicted_state_estimate = self.A * self.current_state_estimate + self.B * control_vector
		predicted_prob_estimate = (self.A * self.current_prob_estimate) * numpy.transpose(self.A) + self.Q
		#--------------------------Observation step-----------------------------
		innovation = measurement_vector - self.H*predicted_state_estimate
		innovation_covariance = self.H*predicted_prob_estimate*numpy.transpose(self.H) + self.R
		#-----------------------------Update step-------------------------------
		kalman_gain = predicted_prob_estimate * numpy.transpose(self.H) * numpy.linalg.inv(innovation_covariance)
		self.current_state_estimate = predicted_state_estimate + kalman_gain * innovation
		# We need the size of the matrix so we can make an identity matrix.
		size = self.current_prob_estimate.shape[0]
		# eye(n) = nxn identity matrix.
		self.current_prob_estimate = (numpy.eye(size)-kalman_gain*self.H)*predicted_prob_estimate
	

class KalmanFilterInd:
	#-------------------------------------------------------------------
	# in order to speedup the calculation and improve the memory usage
	# Here, we assume the expression of genes are independent
	# That means: the expression of gene i at time t g(i,t) only depends:
	# 1) g(i,t-1) the prediction of expression of gene i based on time t-1
	# 2) g_ob(i,t), the observation of expression of gene i on time t.
	# g(i,t) does not rely any other gene j (j!=i). 
	# Based on this assumption, the Kalman Filter can be much simplied) 
	#-------------------------------------------------------------------
	def __init__(self,AT,BT,HT,IT,x,P,Q,R):
		self.AT = AT                    #transition matrix  (simplied to a diagonal matrix each time t) -> can be represented by a vector (dim) for each time t
										# format: List of n vectors, n represents all time points. first time point, its the transition to itelf=>I
										# if only 1 vector is given, it means, it's uniform for all time points
										
										#-------------------------------------------------------
		self.BT = BT                    #transition offset  -> A vector (dim) for each time t
										#  format: List of n vectors, n represents all time points. 
										# if only 1 vector is given, it means, it's uniform for all time points
										#-----------------------------------------------------------------------


		self.HT = HT                     # observation matrix (simplied to a diagonal instead of matrix for each time)-> can be represented by a vector (dim) for each time t
										#  format: List of n vectors, n represents all time points. 
										# if only 1 vector is given, it means, it's uniform for all time points
										
										#-----------------------------------------------------------------------
		self.IT = IT                     # observation offset -> A vector (dim) for each time t
										#  format: List of n vectors, n represents all time points. 
										# if only 1 vector is given, it means, it's uniform for all time points
		
		self.current_state_estimate = x # Initial state estimate.  -> a vector of /mu_s 
		self.current_prob_estimate = P  # Initial covariance estimate. -> a vector of /sigma_s

		self.Q = Q                      # Estimated error in process. -> simplied to a diagonal matrix of vairance
		self.R = R                      # Estimated error in measurements -> simplied to a diagonal matrix of vairance

	def Step(self,A,B,H,I,x):
		# x: x is the one observation 
		dim=len(self.current_state_estimate)  # dim: # of genes in observation	
		# prediction step-----------------------------------------------------------------------------	
		#pdb.set_trace()
		predicted_state_estimate= [A[i]*self.current_state_estimate[i]+B[i] for i in range(dim)] 
		predicted_prob_estimate = [A[i]*self.current_prob_estimate[i]*A[i]+self.Q[i] for i in range(dim)]
		
		# observation step----------------------------------------------------------------------------
		innovation = [x[i]-H[i]*predicted_state_estimate[i]-I[i] for i in range(dim)]
		innovation_covariance =[H[i]*predicted_prob_estimate[i]*H[i]+ self.R[i] for i in range(dim)]

		#correction step-------------------------------------------------------------------------------
		kalman_gain=[predicted_prob_estimate[i]*H[i]*(1.0/innovation_covariance[i]) for i in range(dim)]
		self.current_state_estimate=[predicted_state_estimate[i]+kalman_gain[i]*innovation[i] for i in range(dim)]
		self.current_prob_estimate = [(1-kalman_gain[i]*H[i])*predicted_prob_estimate[i] for i in range(dim)]
		#pdb.set_trace()
		
	def StepMultiple(self,A,B,H,I,x):
		# x: predict and update for multiple observations 
		# Here, we assume, H , I are the same matrix for all observations. Otherwise, 
		# should modify this function. 
		#-----------------------------------------------------------------
		
		# x: x is the mutliple observation 
		dim=len(self.current_state_estimate)  # dim: # of genes in observation	
		# prediction step-----------------------------------------------------------------------------	
		# prediction is the same, 1 transition 
	
		predicted_state_estimate= [A[i]*self.current_state_estimate[i]+B[i] for i in range(dim)] 
		predicted_prob_estimate = [A[i]*self.current_prob_estimate[i]*A[i]+self.Q[i] for i in range(dim)]
		
			
		# observation step----------------------------------------------------------------------------
		# observation different, multiple observation points
		innovation=[(sum([x[j][i] for j in range(len(x))])-len(x)*(H[i]*predicted_state_estimate[i]+I[i]))/len(x) for i in range(dim)]
		innovation_covariance =[H[i]*predicted_prob_estimate[i]*H[i]+ self.R[i] for i in range(dim)]
		
		#correction step-------------------------------------------------------------------------------
		kalman_gain=[predicted_prob_estimate[i]*H[i]*(1.0/innovation_covariance[i]) for i in range(dim)]		
		self.current_state_estimate=[predicted_state_estimate[i]+kalman_gain[i]*innovation[i] for i in range(dim)]
		self.current_prob_estimate = [(1-kalman_gain[i]*H[i])*predicted_prob_estimate[i] for i in range(dim)]
		#pdb.set_trace()
	
	def filter(self,X):
		mu=[]
		sm=[]
		
		for i in range(len(X)):		
			A=self.AT[i]
			B=self.BT[i]
			H=self.HT[i]
			I=self.IT[i]
			self.Step(A,B,H,I,X[i])
			yi=self.current_state_estimate
			si=self.current_prob_estimate
			mu.append(yi)
			sm.append(si)
		return [mu,sm]
		
	def filterMultiple(self,X):
		# X: multiple observation for each time point
		mu=[]
		sm=[]
		
		for i in range(len(X)):		
			A=self.AT[i]
			B=self.BT[i]
			H=self.HT[i]
			I=self.IT[i]
			self.StepMultiple(A,B,H,I,X[i])
			yi=self.current_state_estimate
			si=self.current_prob_estimate
			mu.append(yi)
			sm.append(si)
		return [mu,sm]	
#------------------------------------------------------------------------
		
def smoother(A,B,H,I,xhat,P,Q,R,X):
	def Average(mu1i,sm1i,mu2i,sm2i):
		K=[float(sm1i[j])/(sm1i[j]+sm2i[j]) for j in range(len(mu1i))]
		mui=[mu1i[j]+K[j]*(mu2i[j]-mu1i[j]) for j in range(len(mu1i))]
		smi=[(1-K[j])*sm1i[j] for j in range(len(mu1[i]))]
		return [mui,smi]
		
	# Forward Filter
	FF=KalmanFilterInd(A,B,H,I,xhat,P,Q,R)
	
	[mu1,sm1]=FF.filter(X)
	
	#Backward Filter
	A=[[1.0/item for item in j] for j in A]
	B=[[-1*item for item in j] for j in B]
	xhat=mu1[-1]
	P=sm1[-1]
	
	BF=KalmanFilterInd(A,B,H,I,xhat,P,Q,R)
	[mu2,sm2]=BF.filter(X[::-1])
	mu2=mu2[::-1] # reverse
	sm2=sm2[::-1] # reverse
	#pylab.plot(range(len(X)),[item[0] for item in X],'b',range(len(X)),[item[0] for item in mu1],'r',range(len(X)),[item[0] for item in Y],'g',range(len(X)),[item[0] for item in mu2],'+r')
	#pylab.show()
	#pdb.set_trace()
	mL=[] # mu List
	sL=[] # sL List
	for i in range(len(mu1)):
		[m,s]=Average(mu1[i],sm1[i],mu2[i],sm2[i])
		mL.append(m)
		sL.append(s)
	
	return [mL,sL]
	
#----------------------------------------------------------------------
# this is the smoother for multiple observation for each time point

def smootherMultiple(A,B,H,I,xhat,P,Q,R,X):
	def Average(mu1i,sm1i,mu2i,sm2i):
		K=[float(sm1i[j])/(sm1i[j]+sm2i[j]) for j in range(len(mu1i))]
		mui=[mu1i[j]+K[j]*(mu2i[j]-mu1i[j]) for j in range(len(mu1i))]
		smi=[(1-K[j])*sm1i[j] for j in range(len(mu1[i]))]
		return [mui,smi]
		
	# Forward Filter
	FF=KalmanFilterInd(A,B,H,I,xhat,P,Q,R)
	
	[mu1,sm1]=FF.filterMultiple(X)
	
	#Backward Filter
	xhat=mu1[-1]
	P=sm1[-1]
	B=[[-1*item for item in j] for j in B]
	BF=KalmanFilterInd(A,B,H,I,xhat,P,Q,R)
	[mu2,sm2]=BF.filterMultiple(X[::-1])
	mu2=mu2[::-1] # revese
	sm2=sm2[::-1] # reverse	
	mL=[] # mu List
	sL=[] # sL List
	for i in range(len(mu1)):
		[m,s]=Average(mu1[i],sm1[i],mu2[i],sm2[i])
		mL.append(m)
		sL.append(s)
	#pdb.set_trace()
	return [mL,sL]

#=======================================================================
# other statistical functions 

#----------------------------------------------------------------------
# handling variance ------------------------------------------

def getVarianceVector(x,mu=None):
		# x:the observation of one gene at specific time t
		if mu==None:
			mu=float(sum(x))/len(x)
		v=[(item-mu)**2 for item in x]
		v=sum(v)*1.0/len(x)
		return v
		
def getVariance(X,mu=None):
	# X: observation at specific time point t.  	
	dim=len(X[0]) # number of genes at each time point t
	
	R=[]
	for i in range(dim):
		gi=[item[i] for item in X]
		if mu==None:
			vi=getVarianceVector(gi)
		else:
			vi=getVarianceVector(gi,mu[i])
		R.append(vi)
	return R

def getVarianceT(XT,muT=None):
	# XT: all X for all time T
	# muT: all mu for all time 
	S=[]
	for i in range(len(XT)):
		if muT==None:
			ri=getVariance(XT[i])
		else:
			ri=getVariance(XT[i],muT[i])
		S.append(ri)

	dim=len(XT[0][0])	
	SR=[]
	for i in range(dim):
		vr=[item[i] for item in S]
		vr=sum(vr)/len(vr)
		SR.append(vr)
	return SR
	
def getProcessVariance(X1,X2,A,MU=None):
	# A: vector
	# MU : average at time t-1, vector
	# X1: all observation at time point t
	# X2:  all observations at time point t
	dim=len(X1[0])
	Q=[]
	if MU==None:
		for i in range(dim):
			x1=[item[i] for item in X1] # all observation for gene i at time t-1
			mui=A[i]*sum(x1)/len(x1)
			x2=[item[i] for item in X2] # all observation for gene i at time t
			v=getVarianceVector(x2,mui)
			Q.append(v)
	else:
		for i in range(dim):
			x2=[item[i] for item in X2]
			mui=MU[i]
			v=getVarianceVector(x2,mui)
			Q.append(v)
	return Q
	
def getProcessVarianceT(XT,A,MUT=None):
	S=[]
	# A: List of vector (dim)
	if MUT==None:
		for i in range(len(XT)-1):
			qi=getProcessVariance(XT[i],XT[i+1],A[i])
			S.append(qi)
			
		dim=len(XT[0][0])	
		SR=[]
		for i in range(dim):
			vr=[item[i] for item in S]
			vr=sum(vr)/len(vr)
			SR.append(vr)
	else:
		for i in range(len(XT)-1):
			qi=getProcessVariance(XT[i],XT[i+1],A[i],MUT[i])
			S.append(qi)
			
		dim=len(XT[0][0])	
		SR=[]
		for i in range(dim):
			vr=[item[i] for item in S]
			vr=sum(vr)/len(vr)
			SR.append(vr)
	return SR
	
#=======================================================================

#=======================================================================

def main():   
	#-----------------------------------------------------------------------
	# testing multiple
	N=60
	Noise=0.2
	Y=[[0.39] for item in range(N)]
	X=[[[random.gauss(Y[i][0],Noise)] for i in range(len(Y))] for item in range(6000)] 
	#X=[[[0]]*45,[[1]]*25,[[2]]*5]
	#pdb.set_trace()
	A=[[1]]*len(X)
	B=[[0]]*len(X)
	#B=[[0]]+[[1]]*(len(X)-1)
	H=[[1]]*len(X)
	I=[[0]]*len(X)
	
	mu=[[0.39]]*len(X)
	
	Q=[0]
	R=[0.1]
	#Q=getProcessVarianceT(X,A)                         # process noise 
	#R=getVarianceT(X)                                  # observation noise 
	xhat=[0]
	P=[1]
	
	KF=KalmanFilterInd(A,B,H,I,xhat,P,Q,R)
	FR=KF.filterMultiple(X)
	FR1=smootherMultiple(A,B,H,I,xhat,P,Q,R,X)
# for testing purpose
if __name__ == "__main__":
    main()
    
