#!/usr/bin/env python
"""
tests on delta / a for different alpha values
"""

import os

# values of Nb -> 10,20,50,100,150,200,500,1000
# values of alpha -> 3,5,9,12,15,20
# keep eps @ 1e-8
# NF = 0
# ML = 1000

#p_alpha = [3,5,9,12,15,20,25]
p_alpha = [15,20,25]
N = [10,20,50,100,150,200]#,500,1000]

command = 'python pyfgt_2d.py '

# options:
# -N <N_side>
# -p <alpha>
# --NF <NF>
# --ML <ML>

for i in p_alpha:
	for j in N:
		t_cmd = command + '-s %s -e %e -D %f -N %d -p %d --NF %d --ML %d' %('timeTest_1600',1e-8,0.0001,j,i,0,10000)
		print(t_cmd)
		os.popen(t_cmd)
