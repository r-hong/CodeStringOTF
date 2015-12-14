#! /usr/bin/env python
'''
'''
###################################################
# Rolando P. Hong Enriquez (rolando.hong@iit.it) ##
###################################################
import argparse
import sys
import csv
import math
import os
import pylab as pl
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from functions import dataGenerator 
import random
 
#######################################################
###### Parsing arguments from the command line ########
#######################################################
def _make_parser():
        parser = argparse.ArgumentParser(description="This tool creates N randomly distributed gaussian clouds of data. Each cloud has M points. The final collection of clouds is saved in the output file using 3 columns. This data generator is used in the tutorial to produce fake data for cluster analysis or to estimate the correct number of clusters using nclust.py.",epilog="Thanks for using fakeData.py!!")
        parser.add_argument('N', type=int, help="Number of clusters.")
        parser.add_argument('M', type=int, help="Number of points in each cluster.")
	parser.add_argument('OUTPUT', help="Output filename prefix")
        return parser
#--------------------------------
_p    = _make_parser()
__doc__ += _p.format_help()
######################################################
if __name__ == '__main__':
        args = _p.parse_args()

	points=dataGenerator(args.N,args.M)
	p0=zip(*points)

	# Plot the initial data
	fig = plt.figure()
	ax = fig.add_subplot(111,projection='3d')
	ax.scatter(p0[0],p0[1],p0[2],c='m',marker='o',s=0.5)
	plt.title("Generated 3D data clusters")
	ax.set_xlabel('x')
	ax.set_ylabel('y')
	ax.set_zlabel('z')
	plt.show()

	coords=[]
	coords.append(p0[0])
	coords.append(p0[1])
	coords.append(p0[2])
	coords=zip(*coords)

	print "saving generated data..."
	args.OUTPUT +=`args.N`
	with open(args.OUTPUT, 'w') as file:
		for q in range(len(coords)):
			ss=coords[q]
			ss1=str(ss)[2 : -1]
			#print 'vavava',q,ss1.split(',')
			ss2="".join(ss1.split(','))
			#print 'vavava',q,ss2
			file.write('{}\n'.format(ss2))

