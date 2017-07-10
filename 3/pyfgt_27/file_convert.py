#!/usr/bin/env python

import os
from string import *
import sys

def read_file(fname):
	sArr = []
	f = open(fname,'r')
	for line in f:
		sArr.append(line)
	X = [[] for _ in range(0,len(sArr))]
	for i, s in enumerate(sArr):
		vals = s.split()
		for j in range(0,len(vals)):
			X[i].append(atof(vals[j]))
	f.close()
	return X

def write_file(r,fname):
	outfile = open(fname,'w')
	for i in range(0,len(r)):
		outfile.write('%.20f,%.20f,%.20f\n' %(r[i][0],r[i][1],r[i][2]))

	outfile.close()
	

# main method
# first argument is old file, second argument is new file
fname_old = sys.argv[1]
fname_new = sys.argv[2]
print 'converting %s to %s' %('fname_old','fname_new')

X = read_file(fname_old)
print 'read input file'
write_file(X,fname_new)