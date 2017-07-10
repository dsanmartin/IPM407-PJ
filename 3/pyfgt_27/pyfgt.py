#! /usr/bin/env python
# pyfgt.py
# functions for a 1d FGT

import math
from numpy import *
from pylab import *
from scipy.special import hermite


def fact(n):
	if (n <= 1):
		return 1
	else:
		return n*fact(n-1)

# A_alpha co-efficients
# input: s = array of input sources
#		 q = array of weights
#		 sb  = cluster center
#		 delta
#		 alpha
def A_alpha(s, q, sb, delta, alpha):
	#print len(s), alpha, delta
	A = 0
	sq_delta = math.sqrt(delta)
	# sum over the points
	for i in range(len(s)):
		A = A + q[i]* ((s[i]-sb) / sq_delta)**alpha
	
	return (A / (1.*fact(alpha)))


# h_alpha coefficients
# input: t, n
def h_alpha(t, n):
	if (n == 0):
		temp = 1 #math.exp(-(t**2))
	elif (n ==1):
		temp = 2*t # *math.exp(-(t**2))
	else:
		# recurse
		temp = 2*t*h_alpha(t,n-1) - 2*(n-1)*h_alpha(t,n-2)
	
	print('temp	', n, ' = ',temp)
	return math.exp(-(t**2))*temp
#
#	
#	sign = (-1)**(2*n)
#	H_alpha = sign * (2*t)**n
#
#	return math.exp(-(t**2))*H_alpha

def B_beta(s, q, sb, tc, delta, beta, p):
	sign = ((-1)**beta)/(1.*fact(beta))
	B = 0
	# p is range for alpha
	for i in range(p):
		B_temp = A_alpha(s, q, sb, delta, i)
		x = (sb-tc)/math.sqrt(delta)
		h = exp(-x**2)*hermite(beta+i)(x) #h_alpha((sb-tc)/math.sqrt(delta), i + beta)
		B = B + (B_temp*h)
	
	return (sign*B)

# C_beta == B_beta
def C_beta(s, q, sb, tc, delta, beta, p):
	sign = (-1)**beta / 1.*fact(beta)
	
	C = 0.
	for i in range(p):
		a = A_alpha(s,q,sb,delta,i)
		x = (sb-tc)/math.sqrt(delta)
		h = exp(-x**2)*hermite(i+beta)(x)
		C = C + a*h
	
	return sign*C
	
# take list of 1d points, source points, center, strength & width of gaussian
def Gaussian(pts, s, sb, q, delta):
	r = zeros(len(pts))
	for i in range(len(s)):
		val = []
		for j in range(len(pts)):
			val.append(q[i]*math.exp(-(abs(pts[j]-s[i])**2)/delta))
		val = array(val)
		r = r + val
	return r	

# Parameters:
# s     = list of source points
# q     = list of source weights
# delta = bandwidth of the gaussian
# p     = hermite terms
# p2    = taylor terms
# shift = where the expansion is evaluated
def Gaussian_Hermite(s, q, delta, p, p2, shift):
	sb = mean(s)      # center of source cluster
	tc = sb - shift   # 
	#q_arr = [q,q,q]
	x    = arange(-1.,1.01,0.02) # evaluation points
	real = Gaussian(x, s, sb, q, delta)  # direct evaluation
	# build A_alpha coefficients
	As = []
	for i in range(p):
		As.append(A_alpha(s,q,sb,delta,i))
	
	# approximate the Gaussian
	G = []
	for i in range(len(x)):
		temp = 0.
		for j in range(p):
			temp = temp + As[j]*exp(-(((x[i]-sb)/sqrt(delta))**2))*hermite(j)((x[i]-sb)/math.sqrt(delta)) #h_alpha((x[i]-sb)/math.sqrt(delta),j)
		G.append(temp)
		
	# approximate hermite expansion with taylor expansion
	# already have A_alpha terms
	# get C_alpha coefficients
    Bs = []
    for beta in range(p2):
 		Bs.append(B_beta(s, q, sb, tc, delta, beta, p))
 
	# approximate the Gaussian with the translated series
	G2 = []
	for i in range(len(x)):
		temp = 0.
		for j in range(p2):
			temp = temp + Bs[j]*((x[i]-tc)/sqrt(delta))**j
		G2.append(temp)
	

	real = array(real)
	G = array(G)
	G2 = array(G2)
	error = abs(real-G)+1e-20
	error2 = abs(real - G2) + 1e-20
	print 'max hermite error = ',max(log10(error))
	print 'best taylor error = ',min(log10(error2))
	plot(x,real,'r-')
	plot(x,G,'b-')
	plot(x,G2,'c-')
	#legend(('Direct Solution','Hermite series','Taylor Series'))
	figure()
	plot(x,log10(error),'r-')
	plot(x,log10(error2),'b-')
	#legend(('Error, p = %d' %(p)))
	show()
	
	return real, As, G, error 

s     = [0.01, -0.01]
q     = [1, 0.1]
delta = 0.0001
p     = 20
p2    = 20
evalC = 0.80

Gaussian_Hermite(s, q, delta, p, p2, evalC)
