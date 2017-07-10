#!/usr/bin/env python
"""
read in a set of results and graph errors of similar parameters
"""

import os
from numpy import *
from pylab import *
from string import *
from matplotlib import *

EPS = 1e-20
def main():
	p = [3,5,9,12]

	# N = 10,20,50,100,200
	N = [10,20,50,100,200]
	p3_errors = array([1.8384,0.8764,-0.2037,-1.1950,-2.0204])
	p5_errors = array([2.9357,1.3644,-0.4734,-1.8797,-3.2835])
	p9_errors = array([4.9237,2.2361,-1.1913,-3.6581,-4.5559])
	p12_errors = array([6.0729,2.6364,-1.7340,-5.2974,-4.5550])

	scale_max = max(p3_errors.max(),max(p5_errors.max(),max(p9_errors.max(),p12_errors.max())))
	scale_min = min(p3_errors.min(),min(p5_errors.min(),min(p9_errors.min(),p12_errors.min())))
	print 'scale_max = ',scale_max,'scale_min = ',scale_min

	# plot the individual error plots
	ax=subplot(111)
	#ax.plot(N,p3_errors,'r-')
	ax.set_ylim((scale_min,scale_max))
	#figure()
	plot(N,p3_errors,'r-')
	plot(N,p5_errors,'b-')
	plot(N,p9_errors,'g-')
	plot(N,p12_errors,'k-')
	xlabel('Number of boxes per side')
	ylabel('abs(log10(error))')
	ax.legend(('p=3','p=5','p=9','p=12'),loc='upper right')
	show()


if __name__=="__main__":
	main()