#!/usr/bin/env python
"""
2D FGT from Greengard & Strain
"""

from numpy import *
from math import sqrt, exp
from point2d import *
from scipy.special import hermite
import sys
from string import atof
from pylab import *
from optparse import OptionParser

# read file into lists of point_2d objects
# input - fname - input filename
def read_file(fname):
	sArr = []
	f = open(fname,'r')
	for line in f:
		sArr.append(line)

	X = []
	temp = []
	for i, s in enumerate(sArr):
		vals = s.split()
		for j in range(0,len(vals)):
			temp.append(atof(vals[j]))
		if (len(temp) == 3):
			X.append(point2d(temp[0],temp[1],temp[2],i))
		else:
			X.append(point2d(temp[0],temp[1],0,i))
		temp = []
	f.close()
	return X

# write results into a CSV file
# input - r - results
#		  t - list of original target points
#		  fname - name of output file
def write_file(r,t,fname):
	outfile = open(fname,'w')
	for i in range(0,shape(r)[1]):
		outfile.write('%.20f,%.20f,%.20f\n' %(t[i].x,t[i].y,r[0,i]))           

	outfile.close()

# go from matrix (a,b) to array address
def mat2arr(y,x,yb,xb):
	return (y*xb + x)

# cluster points into uniform boxes
# input - pts - list of points
#		  xb, yb number of clusters in x,y axis
def cluster(pts, xb, yb):
	
	centers = [[] for _ in range(xb*yb)]
	# initialise list for storing points
	B = [[] for _ in range(xb*yb)]
	
	dx = 1. / xb
	dy = 1. / yb
	
	x_bounds = []
	y_bounds = []
	
	# calculate lists of x and y bounds
	# assume domain = [0,1]^2
	for i in range(xb + 1):
		x_bounds.append(i*dx)
	for i in range(yb + 1):
		y_bounds.append(i*dy)
	
	# assign the cluster centers
	for i in range(yb):
		for j in range(xb):
			y_center = (y_bounds[i] + y_bounds[i+1])/2.0
			x_center = (x_bounds[j] + x_bounds[j+1])/2.0

			centers[i*yb+j] = point2d(x_center,y_center)
	
	# now assign points into the appropriate boxes
	for i in range(len(pts)):
		X = pts[i].x
		Y = pts[i].y

		# get y box
		for j in range(yb):
			if (Y >= y_bounds[j] and Y > y_bounds[j+1]):
				# not in this box
				continue
			if (Y >= y_bounds[j] and Y <= y_bounds[j+1]):
				y_box = j
				break
		
		# get x box
		for k in range(xb):
			if (X >= x_bounds[k] and X > x_bounds[k+1]):
				continue
			if (X >= x_bounds[k] and X <= x_bounds[k+1]):
				x_box = k
				break

		B[y_box*yb + x_box].append(pts[i])
	
	# build index array
	A = zeros([xb,yb],dtype=int)
	for i in range(shape(A)[0]):
		for j in range(shape(A)[1]):
			A[i,j] = (i*xb)+j
	
	return A, B, centers

# |alpha| = sum(alpha[i])
def a_abs(alpha):
    a_abs = alpha[0] + alpha[1]
    return a_abs

def factorial(n):
    if n <= 0:
        return 1

    return n*factorial(n-1)

# alpha! = product(alpha[i])
def a_fact(alpha):
    a_fact = factorial(alpha[0])*factorial(alpha[1])
    return a_fact

# generate multi-index vectors
def gen_alpha(p):
    alpha = []
    for i in range(0,p):
        for j in range(0,p-i):
            alpha.append([i,j])

    return alpha

def build_fact_cache(q):
	f = []
	for i in range(q+1):
		f.append(factorial(i))
		
	return f

# build A_alpha coefficients - return array of coefficients
# FGT paper pg. 83, Eqn. 12
# input - s = source points
#		  q = source strengths
#		  sb = center of cluster
#		  delta
#		  alpha - list of alpha multi-index vectors
def A_alpha(s, sb, delta, alpha):
	sq_delta = sqrt(delta)
	A        = []
	# go over each alpha value
	for i in range(len(alpha)):
		alp_temp = alpha[i]
		A_tmp    = 0.
		# sum over all the source points
		for j in range(len(s)):
			x      = (s[j]-sb) / sq_delta
			A_tmp += s[j].gamma*(x)**alp_temp
		A.append(A_tmp / (1.*a_fact(alp_temp)))
	return A

# directly get taylor coefficients for target cluster,
# rather than translate from hermite series
# FGT paper pg. 85, Eqn. 20
# input - s = source points
#		  tc - center of target cluster
#		  delta
#		  beta - vector of multi-index coefficients
#		  H - list of pre-computed hermite polynomials
def B_beta(s,tc,delta,beta,H):
	sq_delta = sqrt(delta)
	B = []
	#print 'tc = ', tc
	for i in range(len(beta)):
		sign = ((-1)**(1.*a_abs(beta[i]))) / (1.*a_fact(beta[i]))
		
		B_tmp = 0.
		# for each point, work out the expansion
		for j in range(len(s)):
			q = s[j].gamma
			x = (s[j]-tc)/sq_delta
			#print 's[j] = ',s[j],'tc = ',tc,'x = ',x
			t = ((x.x)**2) + ((x.y)**2)
			#print 't = ',t
			#print 'x.x = ',x.x,'x.y = ',x.y
			#temp = exp(-t)*hermite(beta[i][0])(x.x)*hermite(beta[i][1])(x.y)
			temp = exp(-t)*(H[beta[i][0]](x.x))*(H[beta[i][1]](x.y))
			B_tmp += q*sign*temp
			#print 'beta[i] = ', beta[i]
			#print 'B_tmp = ',B_tmp
		B.append(B_tmp)
	B = array(B)
	B = B.reshape([1,len(beta)])
	#print shape(B)
	#print 'tc = ',tc
	#print 'B = ',B
	return B

# translate hermite series about sources into taylor series about
# target cluster center - doesn't depend on source pts
# FGT paper pg. 84, Eqn. 18
# return array of values
# input - sb - old cluster center
#		  tc - new center of series
#		  A  - array of A_alpha values for cluster w/center sb
#		  delta
#		  beta - array of multi-index values
#		  H - list of pre-computed hermite polynomials
def C_beta(sb, tc, A, delta, p, q, alpha, beta,H):
	C = []
	for i in range(len(beta)):
		sign  = ((-1.)**(1.*a_abs(beta[i]))) / (1. * a_fact(beta[i]))
		B_tmp = 0.
		# iterate over A_alpha values
		for j in range(len(A)):
			a = A[j]
			x = (sb-tc)/sqrt(delta)
			t = (x.x)**2 + (x.y)**2
			
			h1 = H[alpha[j][0]+beta[i][0]](x.x)
			h2 = H[alpha[j][1]+beta[i][1]](x.y)

			B_tmp += a*exp(-t)*h1*h2
		C.append(B_tmp*sign)
	C = array(C)
	C = C.reshape([1,len(beta)])
	
	return C

# get the interaction list from A about cluster (i,j)
# with a distance of d in x,y directions
# input - A - cluster index array
#		  i - cluster row
#		  j - cluster column
#		  d - number of clusters in each direction
def int_list(A,i,j,d):
	if (i-d < 0):
		y_min = 0
	else:
		y_min = i-d
	
	if (j-d < 0):
		x_min = 0
	else:
		x_min = j-d
	
	return A[y_min:i+d+1,x_min:j+d+1]

# direct calculation of gaussians
# input - s - List of source points
#		  t - list of target points
#		  delta
#		  results - array of results to update
def direct_gaussian(s,t,delta,results):
	for i in range(len(t)):
		r = 0
		# index of the target point - element of results to store in
		idx = t[i].idx
		for j in range(len(s)):
			dist = (s[j]-t[i]).norm()
			v    = dist**2 / delta
			#print 't(%d), s(%d), v = %f' %(i,j,v)
			r   += s[j].gamma*exp(-v)
		results[0,idx] += r
		
	return results
	
# evaluate hermite expansion at a target point
# input - A_alpha - list of series coefficients
#		  t - target point
#		  sb - center of source cluster (series center)
#		  delta
#		  alpha - array of multi-index values
# 		  H - list of pre-computed hermite polynomials
def eval_hermite(A_alpha,t,sb,delta,alpha,H):
	r = 0.
	for i in range(len(alpha)):
		x        = (t-sb)/sqrt(delta)
		pre_mult = A_alpha[i]*exp(-(x.x**2 + x.y**2))
		#r      += pre_mult*hermite(alpha[i][0])(x.x)*hermite(alpha[i][1])(x.y)
		r       += pre_mult*H[alpha[i][0]](x.x)*H[alpha[i][1]](x.y)
	
	return r

# evaluate taylor series at a target point
# input - t - target point
#		  tc - center of taylor series
#		  delta
#		  t_series - taylor series coefficients
#		  beta - array of multi-index values
def eval_taylor(t,tc,delta,t_series,beta):
	x = (tc-t)/sqrt(delta)

	r = 0.
	for i in range(shape(t_series)[1]):

		r += t_series[0,i]*(x**beta[i])
	
	return r

# pre-compute the hermite polynomials to be used throughout
# the calculations (for speed)
# input - p - number of terms kept in hermite / taylor series
def gen_hermite_polys(p):
	H = []
	for i in range(p):
		H.append(hermite(i))
	
	return H

# get nCp needed for number of coefficients for taylor series
# input - n
#		  p
def nCp(n,p):
	top    = factorial(n)
	bottom = factorial(p)*factorial(n-p)
	
	return 1.*(top/bottom)

# sum the gamma values to get Qb
def sum_gamma(s):
	q = 0.
	for i in s:
		q += i.gamma
	
	return q

def truncation_bound(Qb,r,eps):
	pre_mult = exp(-2.*(r**2))
	for n in range(0,10):
		e = pre_mult*exp(-n**2)
		print('%d,%e,%e' %(n,e, eps))
		if (e < eps):
			return n+1
	
	return 10

def hermite_bound(Qb,r,eps):
	K = (1.09**2)
	for p in range(0,50):
		e = K*Qb*(1./factorial(p))*((r**(p+1))/(1-r))
		if (e <= eps):
			return p
	
	return 50

def taylor_bound(Qb,r,eps):
	K = (1.09**2)
	for p in range(0,50):
		e = K*Qb*(1./factorial(p))*((r**(p+1))/(1-r))
		if e <= eps:
			return p
	
	return 50

parser = OptionParser()
parser.add_option("-s","--source_dataset", dest="dataset",type="string",help="input dataset")#,metavar="FILE")
parser.add_option("-D","--delta",type="float",dest="delta",default=0.0001)
parser.add_option("-N","--nside",type="int",dest="N_side",default=50)
parser.add_option("-e","--eps",type="float",dest="eps",help="Desired accuracy",default=1e-6)
parser.add_option("-p","--hermite_terms",type="int",dest="p",help="Number of Hermite series terms",default=None)
parser.add_option("-q","--taylor_terms",type="int",dest="q",help="Number of Taylor series terms",default=None)
parser.add_option("--NF",type=int,dest="NF",help="Direct Evaluation cut-off",default=None)
parser.add_option("--ML",type=int,dest="ML",help="Hermite Evaluation cut-off",default=None)

(options,args) = parser.parse_args()

dataset = options.dataset
delta = options.delta
N_side = options.N_side
eps = options.eps

source_file = 'datasets/'+dataset + '/pts.dat'
target_file = 'datasets/'+dataset + '/targets.dat'

# read files ino lists of points
s = read_file(source_file)
t = read_file(target_file)
	
# setup of parameters
r = 0.5

xb = yb = N_side

Qb = sum_gamma(s)
# terms for hermite expansion
if (options.p == None):
	p = hermite_bound(Qb,r,eps)
else:
	p = options.p
# terms for taylor expansion
if (options.q == None):
	q = taylor_bound(Qb,r,eps)
else:
	q = options.q

# build a factorial cache - only compute this crap once
fact_cache = build_fact_cache(p+q)
# generate table of hermite polynomials
H = gen_hermite_polys(p+q)

# cutoffs
# taylor
#NF = 100
#ML = 0
# direct
#NF = 100
#ML = 100
# hermite
NF = 0
ML = 100
# h2t
#NF = 16
#ML = 10

# cutoffs based on error bounds
if (options.NF == None):
	NF = p+1
else:
	NF = options.NF
if (options.ML == None):
	ML = q+1
else:
	ML = options.ML
d = truncation_bound(Qb,r,eps)

# cluster source & target points
# A is the index matrix, B = source clusters, C = target clusters,
# sb = centers of source clusters
# tb = centers of target clusters
A, B, sb = cluster(s,xb,yb)
A, C, tb = cluster(t,xb,yb)

# generate multi-index vectors
alpha = gen_alpha(p)
beta = gen_alpha(q)

# storage for taylor series & hermite series
h_series = [[] for _ in range(xb*yb)]
t_series = [[] for _ in range(xb*yb)]

# loop over source boxes -> NB > NF: generate hermite series
for i in range(xb):
	for j in range(yb):
		idx = A[i,j]

		NB  = len(B[idx])

		if (NB > NF):
			#print 'generating hermite series for box (%d,%d)' %(i,j)

			h_series[idx] = A_alpha(B[idx],sb[idx],delta,alpha) 
		else:
			continue
			#print 'NB < NF -> no hermite series generated'

# loop over target boxes -> MC > ML: assign storage for taylor series
for i in range(xb):
	for j in range(yb):
		idx = A[i,j]

		MC = len(C[idx])

		if (MC > ML):
			#print 'assigning storage for taylor series in box (%d,%d)' %(i,j)

			# no. terms from iFGT paper, pg. 7
			t_series[idx] = zeros([1,nCp(q+1,2)],dtype=double)
		else:
			continue
			#print 'MC < ML -> no storage assigned'

# array to hold results in
results = zeros([1,len(t)],dtype=double)

direct_eval = 0
hermite_eval = 0
taylor_convert = 0
hermite_taylor = 0
# loop over source boxes 
# - form interaction list
for i in range(yb):
	for j in range(xb):
		# get interaction list for source box (i,j)
		idx  = A[i,j]
		IL   = int_list(A,i,j,d)
		size = shape(IL)
		
		NB = len(B[idx])

		if (NB == 0):
			continue
		elif (NB <= NF):
			# now loop over boxes in interaction list
			for y in range(size[0]):
				for x in range(size[1]):
					t_idx = IL[y,x]
					MC = len(C[t_idx])
					if (MC == 0):
						continue
					elif (MC <= ML):
						# evaluate gaussians directly
						#print 'direct gaussian - s=%d to t=%d' %(idx,t_idx)
						results = direct_gaussian(B[idx],C[t_idx],delta,results)
						#print 'evaluate gaussians directly'
						direct_eval += 1
					else:
						#print 'convert each source to taylor series & add to box'

						t_series[t_idx] += B_beta(B[idx],tb[t_idx],delta,beta,H)
						taylor_convert += 1
		
		# else (NB > NF)
		else:
			# loop over IL
			for y in range(size[0]):
				for x in range(size[1]):
					# hermite expansions already formed
					t_idx = IL[y,x]
					MC    = len(C[t_idx])
					if (MC == 0):
						continue
					elif (MC <= ML):
						# evaluate hermite series at each target & add to accumulated potential
						#print 'Evaluate Hermite & add to accumulated potential'
						for k in range(len(C[t_idx])):
							target_idx = C[t_idx][k].idx
							# eval_hermite(A_alpha,t,sb,delta,alpha)
							results[0,target_idx] += eval_hermite(h_series[idx],C[t_idx][k],sb[idx],delta,alpha,H)
						hermite_eval += 1
					# else (MC > ML)
					else:
						#print 'convert hermite -> taylor & add to taylor for box C'
						# h2t(), C_beta(sb, tc, A, delta, p, q, alpha, beta)
						t_series[t_idx] += C_beta(sb[idx],tb[t_idx],h_series[idx],delta,p,q,alpha,beta,H)
						hermite_taylor += 1

# loop through target boxes evaluating taylor series if MC > ML
taylor_eval = 0
for i in range(yb):
	for j in range(xb):
		idx = A[i,j]
		MC = len(C[idx])

		if (MC == 0):
			continue
		elif (MC > ML):
			taylor_eval += 1
			#print 'evaluate taylor series for C at each target position'
			for k in range(MC):
				target_idx = C[idx][k].idx
				# eval_taylor(t,tc,delta,t_series,beta)
				results[0,target_idx] += eval_taylor(C[idx][k],tb[idx],delta,t_series[idx],beta)

#print 'STATS:'
#print '\tDIRECT EVALUTION - %d BOXES' %(direct_eval)
#print '\tHERMITE EVALUATION - %d BOXES' %(hermite_eval)
#print '\tTAYLOR CONVERT - %d BOXES' %(taylor_convert)
#print'\tHERMITE -> TAYLOR - %d BOXES' %(hermite_taylor)
#print '\tTAYLOR SERIES EVALUATED - %d BOXES' %(taylor_eval)

r = []
x = []
y = []
for i in range(0, shape(results)[1]):
	r.append(results[0,i])
	x.append(t[i].x)
	y.append(t[i].y)

#scatter(x, y, c = r, marker = 's', s = 30)#, edgecolors='none')
#colorbar()
#show()

# create filename
fname = '%s_d%d_p%d_q%d_NF%d_ML%d_N%d_de%.4f.dat' %(dataset,d,p,q,NF,ML,N_side,delta)
#print 'WRITING RESULTS'
write_file(results,t,fname)
