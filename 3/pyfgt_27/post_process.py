#!/usr/bin/env python 
import os
from numpy import *
from pylab import *
from string import *

EPS = 1e-20


def readError(fname):
	sArr = []
	f = open(fname,'r')
	for line in f:
		sArr.append(line)
	X = [[] for _ in range(0,len(sArr))]
	for i, s in enumerate(sArr):
		vals = s.split()
		for j in range(0,len(vals)):
			X[i].append(atof(vals[j]))
	return X


# reading function
def readArray(fname):
	sArr = []
	f = open(fname,'r')
	for line in f:
		sArr.append(line)
	X = [[] for _ in range(0,len(sArr))]
	for i, s in enumerate(sArr):
		vals = s.split(',')
		for j in range(0,len(vals)):
			X[i].append(atof(vals[j]))
	return X

#def readError(fname):
#	e = array(readArray(fname))
#
#	return e

def process_args(argv):
	options = []
	# take out filename
	fname = argv[0]
	arg = []
	if os.path.exists(fname):
		print 'filename = %s' %(fname)
	else:
		print 'filename %s not found! Exiting...' %(fname)
	for args in argv[1:]:
		arg = args.split(',')
	for i in arg:
		options.append([i])

	return fname, options

def contour_plot(pts, path, n, save):
	# read X and Y from path given
	X_pts = readArray(path + "X.dat")
	Y_pts = readArray(path + "Y.dat")

	contour(X_pts,Y_pts,pts[:,:,2],n)
	colorbar()
	show()

def plot_centers(result, xb, yb):
	max_x, min_x = max(result[:,0]), min(result[:,0])
	max_y, min_y = max(result[:,1]), min(result[:,1])
	
	centers = zeros([yb,xb],dtype=complex)
	
	dx = (max_x - min_x) / xb
	dy = (max_y - min_y) / yb
	
	x_bounds = []
	y_bounds = []
	
	for i in range(xb + 1):
		x_bounds.append(min_x + i*dx)
	for i in range(yb + 1):
		y_bounds.append(min_y + i*dy)
		
	#print 'xbounds = ', x_bounds
	#print 'ybounds = ', y_bounds
	# get the cluster centers
	for i in range(yb):
		for j in range(xb):
			y_center = (y_bounds[i] + y_bounds[i+1])/2.0
			x_center = (x_bounds[j] + x_bounds[j+1])/2.0
			centers[i,j] = complex(x_center,y_center)
	
	#print centers
	
	centers = centers.reshape([xb*yb,1])
	
	# plot the centers
	#scatter(centers[:].real, centers[:].imag, marker='x')
	#show()
	
	return centers
	

#def log_error_plot(result, naive, save, name, pts, ry, h, p, xb, yb, eps):
def log_error_plot(result,naive,save,name,pts,p,q,NF,ML,N_side,delta):
	#result = result.reshape([shape(result)[0]**2,3])
	if (N_side > 75):
		# skip this
		skip = True
	else:
		centers = plot_centers(result,N_side,N_side)
	sidelength = int(sqrt(shape(naive)[0]))

	scale = naive[:,2].max()
	naive = naive + EPS
	r = result + EPS
	
	error = (r[:,2] - naive[:,2]) / scale
	error = error + EPS

	minX, maxX = min(result[:,0]), max(result[:,0])
	minY, maxY = min(result[:,1]), max(result[:,1])
	boundaryBox = 1.01 * array([minX, maxX, minY, maxY])
	scatter(r[:,0], r[:,1], c = log10(abs(error)), marker = 's', s = 40, faceted = False, vmin=-16, vmax=0)
	colorbar()
	if (N_side < 75):
		scatter(centers[:].real,centers[:].imag,marker='x')
	axis(boundaryBox)
	#title('Spatial logarithmic error for %d points, ry = %.4f\n h = %.2f, %dx%d boxes, p = %d, eps = %.1e' %(pts,ry,h,yb,xb,p,eps))
	xlabel('x')
	ylabel('y')

	max_err = log10(abs(error).max())
	min_err = log10(abs(error).min())
	print 'Max relative logarithmic error (max f) = %.4f' %(max_err)
	print 'Min relative logarithmic error (max f) = %.4f' %(min_err)

	# work out L2 error
	L2 = 0
	for i in range(len(error)):
		L2 += error[i]**2
	L2 = sqrt(L2)
	print 'L2 Error = %e' %(L2)
	
	if (not save):
		show()
	else:
		print 'saving!'
		savefig(name+'.png')

# histogram of cluster numbers
def cluster_histogram(results, save):
	#print 'Clusters histogram'
	no_cl = results[:,:,3]
	clf()
	hist(no_cl,100)
	if (save):
		#print 'Saving cluster histogram'
		show()
	else:
		show()

if __name__ == "__main__":
	# get the results filename & sort everything
	fname, opts = process_args(sys.argv[1:])
	# read results file and reshape -- assume square data file
	results = array(readArray(fname))
	sidelength = sqrt(shape(results)[0])
	#[s,s,3] for naive, [s,s,4] for results
	#results = results.reshape([sidelength,sidelength,3])
	
	# get data from the filename, i.e.
	# fname = '../results/timeTest_10000_ry0.374196_h0.1_p20_x10_y10.dat'
	fname = os.path.basename(fname)
	# remove the .dat from the end
	name = fname.split('.dat')[0]
	out_name = name
	name = name.split('_')
	data_name = name[0]
	no_pts = atoi(name[1])
	d = atof(name[2][1:])#len(name[2])])
	p = atof(name[3][1:])#len(name[3])])
	q = atoi(name[4][1:])#len(name[4])])
	NF = atoi(name[5][2:])
	ML = atoi(name[6][2:])
	N_side = atoi(name[7][1:])
	delta = atof(name[8][2:])
	
	print 'pts = %d, d = %d, p = %d, q = %d, NF = %d, ML = %d, N_side = %d, delta = %.4f' %(no_pts,d,p,q,NF,ML,N_side,delta)
	h = sqrt(delta)
	folder = 'datasets/%s_%d/' %(data_name,no_pts)
	
	err_name = folder + 'h%.4f_naive.dat' %(h)
	source_name = folder + 'pts.dat'
	print 'error name = %s\n' %(err_name)

	# read error values
	# needs [..] to work
	save = 0
	if ["save"] in opts:
		save = 1
	if ["contour"] in opts:
		contour_plot(results,folder,300,save)
	elif ["log"] in opts:
		error = readError(err_name)
		log_error_plot(results,error,save,out_name,no_pts,ry,h,p,x_boxes,y_boxes,eps)
	elif ["cluster"] in opts:
		print 'cluster plot'
	elif ["hist"] in opts:
		cluster_histogram(results,save)
	elif ["centers"] in opts:
		sources = array(readArray(source_name))
		plot_centers(sources,x_boxes,y_boxes)
	else:
		error = array(readError(err_name))
		#log_error_plot(results,error,save,out_name,no_pts,ry,h,p,x_boxes,y_boxes,eps)
		log_error_plot(results,error,save,out_name,no_pts,p,q,NF,ML,N_side,delta)

