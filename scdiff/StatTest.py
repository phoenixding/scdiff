"""
@Name: StatTest
@Author:Jun Ding
@Date: Mar.3,2011
@Version: 1.0
module 1: HyperGeometricTest
function: perform hypergeometric testing to get the p-value
usage:
python HyperGeometricTest(N,M,n,m)

module 2: pbinom
function:compute the cumulative probability densiuty function of the binomial distribution up P(X<=x)

"""
import math,pdb,sys,os
from scipy.stats import binom_test

def logc(a,b):
	s=0
	for k in range(b):
		s=s+math.log(float(a-k)/float(b-k))
	return s

def p(N,M,n,m,t):
	bot=logc(N,M)
	top1=logc(n,t)
	top2=logc(N-n,M-t)
	r=top1+top2-bot
	r=math.exp(r)
	return r

def HyperGeometricTest(N,M,n,m):
		p_i=0
		stop=min(M,n)
		for t in range(m,stop+1):
			p_i=p_i+p(N,M,n,m,t)
		return min(p_i,1.0)


#=======================================================================
# sub-routine for pbinom 
def erf(z):
        t = 1.0 / (1.0 + 0.5 * abs(z))
        # use Horner's method
        ans = 1 - t * math.exp( -z*z -  1.26551223 +
                                                t * ( 1.00002368 +
                                                t * ( 0.37409196 + 
                                                t * ( 0.09678418 + 
                                                t * (-0.18628806 + 
                                                t * ( 0.27886807 + 
                                                t * (-1.13520398 + 
                                                t * ( 1.48851587 + 
                                                t * (-0.82215223 + 
                                                t * ( 0.17087277))))))))))
        if z >= 0.0:
                return ans
        else:
                return -ans

def normal_estimate(s, p, n):
    u = n * p
    o = (u * (1-p)) ** 0.5

    return 0.5 * (1 + erf((s-u)/(o*2**0.5)))
 
#======================================================================
 
def pbinom(x,n,p):
	
	# handling special cases
	if x<0:
		return 0
	if n<=0:
		return 0
	if x>n:
		return 1
		
	# use scipy.binom_test to calculate binomial test p-value
	pv=binom_test(x,n,p,alternative="less")
	if (1-pv)<=sys.float_info.epsilon/2:
		return 1
	else:
		return pv
 
def pbinom1(x,n,p):
	# this is approximation
	# if n is larger (<2000), approximation 1
	if n<2000:
		q=1.0-p
		pdf=cdf=q**n
		f=p/q
		for i in range(1,x+1):
			pdf*=((n-i+1.0)/i*f)
			cdf+=pdf
		return cdf
	else:
	# if n>=2000 (relatively large, approximiation 2
		return normal_estimate(x,p,n)
		
		

def binomial(G,k):
	if k==0 or k==G:
		return 1
	elif k>G:
		print('error!')
	else:
		fact=1
		for i in range(0,k):
			up=G-i
			down=k-i
			divid=float(up)/down
			fact=fact*divid
		return fact

def pc(module,lamb,w):
	pc=1
	for i in range(len(module)):
		lam=float(lamb[module[i]])
		pr=max(0.0,1-math.exp(-w*lam))
		pc=pc*pr
	return pc

