#! /usr/bin/env python
'''
'''
###################################################
# Rolando P. Hong Enriquez (rolando.hong@iit.it) ##
###################################################
import argparse
import pickle 
import sys
import math
from functions import load_data, load_data_range,loadpcurve 
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from collections import defaultdict, Counter
from mpl_toolkits.mplot3d import Axes3D
from itertools import product, combinations

#######################################################
###### Parsing arguments from the command line ########
#######################################################
def _make_parser():
parser = argparse.ArgumentParser(description="Estimates the density of points on the principal curve in the collective variable space (only 3D case, e.g., coordinates of the center of mass). As the data to estimate the density on the curve is necesarily discrete, here we extrapolate the density pattern from larger distances until reaching the points on the principal curve. Shortly, we use an equation like ln(d) = ln(C)+ln(e(-lr)) = ln(C)-lr, where d=density, r=radius, l=decay constant. We estimate the slope of the model and extrapolate for the case in which the distance to the curve is zero, that is, the density 'on the curve'. The procedure is repeated for each point of the principal curve to render a density profile that is both, plotted and saved.",epilog="Thanks for using pointDensity3D.py!!")
parser.add_argument('file1', help="Input pcurve_state file.")
parser.add_argument('file2', help="Original CV file.")
parser.add_argument('file3', help="Output file.")	
parser.add_argument('choice', type=int, help="Use Min or Mout? choice '0' for Min and '1' for Mout.")
args = parser.parse_args()

#######################################
############ CHECK_PCURVE ################
####################################### 
#File with the pcurv estate (output from pcurve.py)
## loading pcurve state 
N, DIM, Min, Mout, smooth, Niter, Ps, Psout, converg, free_energy, allpoints, closeP, closeI, closePP, closeII = loadpcurve(args.file1)
P1=load_data(args.file2,N,DIM)
print '#######################################################'
print '############# pcurve run parameters  ##################'
print '#######################################################'
print '                                               '
print 'N      (Size of the CV vectors)            = ',N
print 'DIM    (# of CVs)                          = ',DIM
print 'Min    (#of points on the string)          = ',Min
print 'Mout   (# of final reparametrized points)  = ',Mout
print 'smooth (strength of the smoothing [0-1])   = ',smooth
print 'Niter  (# of iterations)                   = ',Niter
print '                                             '
print '#######################################################'

if args.choice==0:
	PPP=Ps
if args.choice==1:
	PPP=Psout
# calculate the minimal/maximal distance in the set of points
maxi=0.0
alldist=[]
for i in range(N-1):
	for j in range(i+1,N):
		dist=np.linalg.norm(P1[i]-P1[j])
		alldist.append(dist)
		if dist>maxi:
			maxi=dist

mini=maxi
meanD=sum(alldist)/len(alldist)
for i in range(N-1):
        for j in range(i+1,N):
                dist=np.linalg.norm(P1[i]-P1[j])
                if dist<mini:
                        mini=dist

deltaR=mini
max_gap=mini*5.0

#find the maximal density for point A
NR=int(meanD/mini)
##########################################################
##########################################################
##########################################################
allDmax=[]
for ww in range(len(PPP)):
        Point=PPP[ww]
	print 'Point with index', ww, 'out of ',(len(PPP)-1)
	# calculate the NR densities
	densityI=[]
	r0=0.0
	rall=[]
	for i in range(NR):
		c=0
		r0=r0+deltaR
		rall.append(r0)
		for j in range(N):
			dist=np.linalg.norm(P1[j]-Point)
			if dist<r0:
				c=c+1
		vol=(4/3)*np.pi*((r0)**3)
		dens=(c/vol)	
		densityI.append(dens)

	#choose the points for the fitting
	for i in range(len(densityI)):
		if densityI[i] != 0:
			I1=i
			break
	I1=I1+1
	presentD=densityI[I1]
	nextD=densityI[I1+1]
	fitD=[]
	fitR=[]
	while presentD>nextD:
		presentD=densityI[I1]
		nextD=densityI[I1+1]
		I1=I1+1
		fitD.append(nextD)
		fitR.append(I1+1)

	fitDD=fitD[0:(len(fitD)-1)]
	fitRR=fitR[0:(len(fitR)-1)]

	#plot function ln(d)=ln(C)+ln(e(-lr))=ln(C)-lr
	#where d=density, r=radius, l=decay constant 
	#(1) get ln(d)
	logDD=[]
	for i in range(len(fitDD)):
		dl=math.log(fitDD[i])
		logDD.append(dl)

	#(2) calculate the slopes
	slopeT=[]
	for i in range(1,len(logDD)):
		slope=(logDD[i]-logDD[i-1])/(fitRR[i]-fitRR[i-1])
		slopeT.append(slope)
	m=(sum(slopeT))/len(slopeT)	
	FF=logDD[0]-(m)*fitRR[0]
	Dmax=math.exp(FF)
	allDmax.append(Dmax)


#save 
with open(args.file3, 'w') as file:
        for item in allDmax:
                file.write('{}\n'.format(item))


XX=[]
for i in range(1,len(allDmax)+1):
        XX.append(i)
pl.ylabel('Density on the curve')
pl.xlabel('Points on the Curve')
pl.xlim(1,len(allDmax))
pl.ylim(0,max(allDmax)+0.1)
pl.plot(XX,allDmax,c='k', linewidth=2)
pl.fill_between(XX, allDmax, facecolor='blue', alpha=0.1)
pl.show()

