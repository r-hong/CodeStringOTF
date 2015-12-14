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
from functions import file_lines, file_cols, load_data, maxvec, string_ini, voroid_mean 
import random
import pickle 

#######################################################
###### Parsing arguments from the command line ########
#######################################################
def _make_parser():
	parser = argparse.ArgumentParser(description="Calculates the optimal number of clusters in a dataset formed by a group of 'M' collective variables (columns in the input file) measured 'N' times (rows in the input file). The procedure evolves an initial guess formed by K cluster centers. At each new iteration, the new centers are estimated as the mean value of the M-dimensional voronoi cells (k-means method). Next we calculate the 'distortion' of the K cluster partition (dKi). We repeat the same procedure for a range of cluster partitions Ki (2 < Ki < maxC). The highest 'jump' in the cluster partitions (Jk_i = dK_i**(-M/2) - dK_(i-1)**(-M/2)) represents the optimal partition and therefore the correct number of clusters in the data set.",epilog="Thanks for using nCluster.py!!")
	parser.add_argument('INPUT', help="Input file")
	parser.add_argument('maxC', type=int, help="Maximal mumber of clusters (2 < Ki < maxC).")
	parser.add_argument('Niter', type=int, help="Number of iterations updating the strings.")
	parser.add_argument('-CV', action='append', dest='CV2plot', default=[], help="Indexes between [0, (DIM-1)] of the CVs to be visualized (2D and 3D options are possible).")
	parser.add_argument('-plotdata', action='append', dest='VAL', default=[], help="Default is 'yes'. If the value 'no' is specified then the intermediate plots will be produced with the curves only.")
	return parser
#--------------------------------
_p    = _make_parser()
__doc__ += _p.format_help()
######################################################
if __name__ == '__main__':
        args = _p.parse_args()

	# Error handling for the visualization
	#######################################
	if args.CV2plot==[]:
		sys.exit('Error (-CV): Please enter the CVs to be analysed.')
	if (len(args.CV2plot)==1) or (len(args.CV2plot)>3):
		sys.exit('Error (-CV): only 2D or 3D visualizations are possible.')	

	# Getting axis labels for visualization
	#######################################
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

	# Error handling for the optional -plotdata argument 
	####################################################
	plotdata='yes'
	if (args.VAL==[]):
		plotdata='yes'
	else:
		dd=args.VAL[0]
		if dd=='no':
			plotdata=dd
		elif (dd=='yes') or (args.VAL==[]):
			plotdata='yes'
		else:
			sys.exit("Error (-plotdata): Only 'yes' or 'no' values are admisible.")
	
	##########################
	####### MAIN #############
	##########################
	sys.stdout.write("\a")
	#loading data points ('N' points for each one of the 'DIM' CVs)
	N=file_lines(args.INPUT)
	DIM=file_cols(args.INPUT)

	P1=load_data(args.INPUT,N,DIM)
	dK1=[]
	dK2=[]
	conv_all=[]
	for k in range(2,args.maxC):
		print 'Cluster size: ', k, 'out of ', args.maxC
		# estimate the initial string (P2) of size Min
		P2 = string_ini(P1,DIM,k)
		Pini=string_ini(P1,DIM,k)
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
	
		ss1=np.asarray(allpoints)
		ss2=np.reshape(ss1,(1,len(allpoints)))
		PPm=np.asarray(Pm)

		#Getting the 'distortion'
		msum=0
		for i in range (len(allpoints)):
			vorcell=ss2[0,i]
			dist=0
			for j in range (len(vorcell)):
				dist=dist+(np.linalg.norm(PPm[i,:]-P1[vorcell[j]]))
			msum=msum+(dist/len(vorcell))
		distortion1=((msum/len(allpoints))/DIM)**(-1*(DIM/2))
		distortion2=((msum/len(allpoints)))**(-1*(DIM/2))
		dK1.append(distortion1)
		dK2.append(distortion2)
		conv_all.append(converg)

	jump1=[]
	jump2=[]
	for i in range(len(dK1)-1):
		jj1=dK1[i+1]-dK1[i]
		jj2=dK2[i+1]-dK2[i]
		jump1.append(jj1)	
		jump2.append(jj2)

	#plot dK1
	x=[]
	for i in range(2,args.maxC):
	        x.append(i)
	pl.ylabel('dK1[i]')
	pl.xlabel('Number of clusters')
	pl.xlim(1,args.maxC)
	pl.ylim(0,max(dK1)+0.1)
	pl.title('Distortion')
	pl.scatter(x,dK1,c='k', marker='o')
	pl.plot(x,dK1,c='k', linewidth=3)
	pl.fill_between(x, dK1, facecolor='blue', alpha=0.1)
	pl.show()

	#plot jump1
	x=[]
	for i in range(3,args.maxC):
		x.append(i)
	pl.ylabel('dK1[i]-dK1[i-1]')
	pl.xlabel('Number of clusters')
	pl.xlim(1,args.maxC)
	pl.ylim(0,max(jump1)+0.1)
	pl.title('Jump in the distortion')
	pl.scatter(x,jump1,c='k', marker='o')
	pl.plot(x,jump1,c='k', linewidth=3)
	pl.fill_between(x, jump1, facecolor='blue', alpha=0.1)
	pl.show()

