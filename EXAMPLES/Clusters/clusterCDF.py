#! /usr/bin/env python
"""
"""
import argparse
import sys
import csv
import math
import os
import pylab as pl
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from functions import * 
import random
import pickle 
import scipy.stats as s

def _make_parser():
	parser = argparse.ArgumentParser(description="This tool estimates and plots the positions of the K clusters centers in the data. The cumulative distributions functions for the values associated to each cluster are also calculated to determine the temporal order of the clusters (this assumes that the indices of the values for the data points or colective variables corrispond to time as it is for instance implicit when we collect colective variable values in time during a molecular dynamics simulation. Resuming, from this calculation we obtain the reordering of the clusters in a sequence that should be follow for the calculation of principal curves between pairs of cluster. In this tool we first evolve an initial guess formed by K points to find the K center of the clusters. At each new iteration, the new points are estimated as the mean value of the DIM-dimensional voronoi cells, after this, cumulative distribution functions are calculated and ploted.",epilog="Thanks for using clustersCDF.py!!")
	parser.add_argument('IN_file', help="Input file")
	parser.add_argument('K', type=int, help="Number of clusters.")
	parser.add_argument('Niter', type=int, help="Number of iterations updating the strings.")
	parser.add_argument('-CVi', action='append', dest='CV2plot', default=[], help="Indexes between [0, (DIM-1)] of the CVs to be visualized (2D and 3D options are possible).")
	parser.add_argument('-plotdata', action='append', dest='VAL', default=[], help="Default is 'yes'. If the value 'no' is specified then the intermediate plots will be produced with the curves only.")
	return parser
#--------------------------------
_p    = _make_parser()
__doc__ += _p.format_help()
######################################################
if __name__ == '__main__':
        args = _p.parse_args()

	#####################################################################################################
	# Error handling for the visualization
	#######################################
	if args.CV2plot==[]:
		sys.exit('Error (-CVi): Please enter the CVs to be analysed.')
	if (len(args.CV2plot)==1) or (len(args.CV2plot)>3):
		sys.exit('Error (-CVi): only 2D or 3D visualizations are possible.')	

	#Getting axis labels for visualization
	######################################
	lab=[]
	if (len(args.CV2plot)==2):
		exec "lab.append('CV %s')" % (args.CV2plot[0])
		exec "lab.append('CV %s')" % (args.CV2plot[1])
		CVX=int(args.CV2plot[0])
		CVY=int(args.CV2plot[1])
		labx=lab[0]
		laby=lab[1]
	else:
		exec "lab.append('CV %s')" % (args.CV2plot[0])
		exec "lab.append('CV %s')" % (args.CV2plot[1])
		exec "lab.append('CV %s')" % (args.CV2plot[2])
		CVX=int(args.CV2plot[0])
		CVY=int(args.CV2plot[1])
		CVZ=int(args.CV2plot[2])
		labx=lab[0]
		laby=lab[1]
		labz=lab[2]
	# error handling for the optional -plotdata argument 
	plotdata='yes'

	if (args.VAL==[]):
		plotdata='yes'
	else:
		dd=args.VAL[0]
		if dd=='no':
			#print 'no data willbe plotted...'
			plotdata=dd
		elif (dd=='yes') or (args.VAL==[]):
			#print 'plot everything...'
			plotdata='yes'
		else:
			sys.exit("Error (-plotdata): Only 'yes' or 'no' values are admisible.")
	
	####### MAIN ###################################################################
	sys.stdout.write("\a")
	#loading data points ('N' points for each one of the 'DIM' CVs)
	N=file_lines(args.IN_file)
	DIM=file_cols(args.IN_file)

	P1=load_data(args.IN_file,N,DIM)
	dK1=[]
	dK2=[]
	conv_all=[]

	# estimate the initial string (P2) of size Min
	P2 = string_ini(P1,DIM,args.K)
	Pini=string_ini(P1,DIM,args.K)
	PiniT=zip(*Pini)

	# maximun distances of the data to the images of the string
	maxdistances=maxvec(P1,P2)

	# Estimate the mean of the voronoid cells
	P3, fail, free_energy, allpoints, closeP, closeI = voroid_mean(P1,P2,DIM,maxdistances,0)

	#######################################################
	if (len(args.CV2plot)==3):
		converg=[]
		for i in range(args.Niter):
			print 'Iteration ', (i+1), 'out of ', args.Niter
			# get voronoid centers
			maxdistances=maxvec(P1,P3)
			Pm, fail, free_energy, allpoints, closeP, closeI = voroid_mean(P1,P3,DIM,maxdistances,0)
			d=np.linalg.norm(np.asarray(Pm)-np.asarray(P3))
			converg.append(d)
			P3=Pm

	else:
       		converg=[]
       		for i in range(args.Niter):
               		print 'Iteration ', (i+1), 'out of ', args.Niter
               		# get voronoid centers
               		maxdistances=maxvec(P1,P3)
               		Pm, fail, free_energy, allpoints, closeP, closeI = voroid_mean(P1,P3,DIM,maxdistances,0)
			d=np.linalg.norm(np.asarray(Pm)-np.asarray(P3))	
			converg.append(d)
			P3=Pm
	
	############################################################################
	##### PLOTING of the final curve with data #####
	if fail==1:
		print 'Failed run... (no output)'
	else:
		PP1=zip(*P1)
		PP2=zip(*P2)
		PP3=zip(*P3)
		#ploting
	if (len(args.CV2plot)==3):

		#if args.DIM==3:
		fig = pl.figure()
		ax = fig.add_subplot(111,projection='3d')
        	ax.set_xlabel(labx)
       		ax.set_ylabel(laby)
       		ax.set_zlabel(labz)
       		ax.set_title('Cluster centers', alpha=0.6)
		ax.scatter(PP1[CVX],PP1[CVY],PP1[CVZ],c='m',marker='o',s=0.5,alpha=0.5)
		ax.scatter(PP3[CVX],PP3[CVY],PP3[CVZ],c='r',marker='o',s=100)
		pl.show()
	else:
	#	if args.DIM==2:
		pl.title('Cluster centers',alpha=0.6)
		pl.scatter(PP1[CVX],PP1[CVY],c='m',marker='o',s=1.5,alpha=0.5)
		pl.scatter(PP3[CVX],PP3[CVY],c='r',marker='o',s=100)
        	pl.xlabel(labx)
       		pl.ylabel(laby)
		pl.show()

	#plot
	pl.ylabel('Cummulative Prob.')
	pl.xlabel('MD time')
	pl.title('CDF')	
	for i in range(args.K):
		exec "data_sorted%s=np.sort(allpoints[%s])" % (i,i)
		exec "p%s = 1. * np.arange(len(allpoints[%s])) / (len(allpoints[%s]) - 1)" % (i,i,i)
		exec "pl.plot(data_sorted%s,p%s,linewidth=2,label='cluster%s')" % (i,i,i)
	pl.legend()
	pl.show()
