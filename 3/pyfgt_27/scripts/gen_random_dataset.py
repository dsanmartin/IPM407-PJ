#!/usr/bin/env python

"""
Generate a random dataset in [0,1]^2 - generates sources & targets
"""

#from random import random
import sys
from optparse import OptionParser
from numpy import *

def write_file(s,fname):
	outfile = open(fname,'w')
	for i in s:
		if (len(i) == 2):
			outfile.write('%.20f %.20f\n' %(i[0],i[1]))
		elif (len(i) == 3):
			outfile.write('%.20f %.20f %.20f' %(i[0],i[1],i[2]))
	
	outfile.close()

def gen_sources(N):
	#r = zeros([1,N],dtype=array)
	r = [[] for _ in range(N)]
	#print shape(r)
	for i in range(len(r)):
		X = random.random()
		Y = random.random()
		gamma = random.random()
		r[i] = [X,Y,gamma]
	
	return r

def gen_targets(N):
	r = [[] for _ in range(N)]
	for i in range(N):
		X = random.random()
		Y = random.random
		r[i] = [X,Y]
	
	return r

N=100000
s = gen_sources(N)
t = gen_targets(N)

# choose dataset name
s_fname = 'pts.dat'
t_fname = 'targets.dat'

write_file(s,s_fname)
write_file(t,t_fname)
