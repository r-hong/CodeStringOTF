#! /usr/bin/env python
###################################################
# Rolando P. Hong Enriquez (rolando.hong@iit.it) ##
###################################################
import sys
import argparse
import csv
import math
import os
#import matplotlib.pyplot as plt
import pylab as pl
import numpy as np
#from stringfunc import normvect, diffstring
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import random
import pickle 

#################################################
def getX2near(X1,index1,dist_to_X2,Pin,DIM):
#################################################################################################
# Find point X2 lying in the line between the DIM-dimensional points X1 and the point Pin[index1]
# located at distance dist_to_X2 from X1.  
#################################################################################################
        m1=X1
        m2=Pin[index1]
        tmp=[]
        for i in range(0,DIM):
                tmp.append('')
                for i, value in enumerate(tmp):
                        exec "H%s=[]" % (i)

        PP = [[0]*DIM]
        P=zip(*PP)
        dm1m2=np.linalg.norm(np.asarray(m1)-np.asarray(m2))
        if (dm1m2<dist_to_X2):
                sys.exit('Error (getX2near): distance to X2 is larger than the distance to the next point in the input vector.')
        else:
                for i in range(DIM):
                        step=(((math.fabs(m1[i]-m2[i]))*dist_to_X2)/dm1m2)
                        if (m1[i]<m2[i]):
                                coord=(m1[i] + step)
                        else:
                                coord=(m1[i]-step)

                        exec 'H%s.append(coord)' % (i)
                GG=[]
                for i in range(DIM):
                        exec 'GG.append(H%s)' % i
                Xa=np.asarray(GG)
                Ga=np.reshape(Xa,(1,DIM))
		X2=Ga[0]
		index2=index1
                return X2, index2

def getX2far(X1,index1,dist_to_X2,Pin,DIM):
# add description here!!!
        # check if dist_to_X2 is too long
        fail=0
        finalSum=0
        dc=[]
        for i in range(index1,len(Pin)-1):
                d=np.linalg.norm(Pin[i]-Pin[i+1])
                dc.append(d)
                dt=sum(dc)
        if (dist_to_X2>dt):
                fail=1
        del d, dc, dt
        if fail==1:
                sys.exit('Error (getX2far): the distance to the new point should be smaller than the remaining length of the vector.')
        else:
                #find index2 and X2
                dc=[]
                for i in range(index1,len(Pin)):
                        d=np.linalg.norm(X1-Pin[i])
                        dc.append(d)
                        dt=sum(dc)
                        if dt<dist_to_X2:
                                X1=Pin[i]
                                index2=i+1
                        else:
                                finalSum=sum(dc[0:(len(dc)-1)])
                                break
                residualD=dist_to_X2-finalSum
		X2, indextmp = getX2near(Pin[index2-1],index2,residualD,Pin,DIM)
                #tmp1, indextmp = getX2near(Pin[index2-1],index2,residualD,Pin,DIM)
                #tmp2=np.asarray(tmp1)
                #X2=np.reshape(tmp2,(3,1))
                return X2, index2

####################################################################
def load_data(filename,N,DIM):
#load N points of dimension DIM. 
#data should be in a file with N rows and DIM columns
        tmp=[]
        for i in range(0,DIM):
                tmp.append('')
                for i, value in enumerate(tmp):
                        exec "coord%s=[]" % (i)

        #### reading file ########
        c=0
        with open(filename) as file:
                reader = csv.reader(file, delimiter=' ')
                for row in reader:
                        for i in range(0,DIM):
                                exec "coord%s.append(float(row[i]))" % (i)
        Plist = [[0]*DIM]*N
        for i in range(N):
                T=[]
                for j in range(DIM):
                        exec 'T.append(coord%s[%s])' % (j,i)

                Plist[i]=T
        Parray=np.asarray(Plist)
        return Parray
###############################################################
def load_data_range(filename,N,CVini,CVend):
#load N points from columns N1 to N2 in a data file
##############################################################
	# check total number of columns (or CVs) in the input file
	with open(filename) as f:
		reader = csv.reader(f, delimiter=' ', skipinitialspace=True)
		first_row = next(reader)
		num_cols = len(first_row)	

	#error handling
	if (CVini >= CVend):
		sys.exit('Error (load_data_range): The index of the first (CVini) should be smaller than the index of the last CV (CVend).')
	if (CVini>num_cols) or (CVend>num_cols):
		print "The CV_file only contains", num_cols, "CVs."
		sys.exit('Error (load_data_range): The indexes of the CVs cannot be larger the the total number of CVs.')

	#initialization
	DIM=CVend-CVini+1
        for i in range(CVini,CVend+1):
		exec "coord%s=[]" % (i)

        #### reading file ########
        c=0
        with open(filename) as file:
                reader = csv.reader(file, delimiter=' ')
                for row in reader:
                        for i in range(CVini,CVend+1):
                                exec "coord%s.append(float(row[i]))" % (i)
        Plist = [[0]*DIM]*N
        for i in range(N):
                T=[]
                for j in range(CVini,CVend+1):
                        exec 'T.append(coord%s[%s])' % (j,i)

                Plist[i]=T
        Parray=np.asarray(Plist)
        return Parray

##################################################################
def maxvec(P1,P2):
# calculates the maximal distance among all points of vector P1 
# to each point of vector P2. Return the vector of maximal distances
# of size len(P2).
        maxvec=[]
        for j in range(len(P2)):
                maximal=0
                dist=[]
                for i in range(len(P1)):
                        dist=np.linalg.norm(P1[i]-P2[j]) 
                        if dist > maximal:
                                maximal=dist
                                index_max=i
                maxvec.append(maximal)
        return maxvec

###########################################################################################
def voroid_mean(P1,P2,DIM,maxdistances):
###################################################
# Obtaining the mean values of each voronoid cell # 
###################################################
# input:
# DIM = dimension of the points
# maxdistances = vector of maximal distances calculated with the function maxvec
# P1  = original data
# P2  = string vector, for each point of P2 (P2[i])we calculate a voronoid cell:
# the points of P1 closer to P2[i] than to any other point of P2.
# output:
# P3= mean points of each voronoid cell point (eachpointwith dimension DIM)
# fail = binary number. fail=1 means that there were empty voronoid cells, and the procedure
# should be restarted with other conditions(e.g., more points on P1 or a different P2)
##############################################################################################	
        kB=0.0019872041            #Boltzmann constant (kcal/mol/K)
        T=300                      #temperature (K) 
	fail=0
	tmp=[]
	for i in range(0,len(P2)):
		tmp.append('')
		for i, value in enumerate(tmp):
			exec "vor%s=[]" % (i)
			exec "ff%s=[]" % (i)
			exec "points%s=[]" % (i)
	#find the points on each cell
	for i in range(len(P1)):
		for j in range(len(P2)):
			dist=np.linalg.norm(P1[i]-P2[j])
			if j==0:
				minid=maxdistances[j]
				index=0
			if dist < minid:
				minid=dist
				index=j
		exec "vor%s.append(P1[%s])" % (index,i)
		exec "points%s.append(%s)" % (index,i)
	#put the indexes of the points within each cell in a list
	#e.g., allpoints[0] will contain the indexes of the points
	#within the first voronoi cell
	allpoints=[]
	for j in range(len(P2)):
		exec 'vorT%s=zip(*vor%s)' % (j,j)
		exec 'allpoints.append(points%s)' % (j)

	c=0
	for j in range(len(P2)):
		dd=[]
		exec 'dd.append(vorT%s)' % (j)
		if dd==[[]]:
			c=c+1
			fail=1
		else:
			for i in range(DIM):
				exec 'a=np.asarray(vorT%s)' % (j)
				b=sum(a[i])/len(zip(*a))
				exec "ff%s.append(b)" % (j)	

	if c>0:
		print 'Error (voroid_mean): There are', c, 'empty voronoid cells'
		sys.exit('Error (voroid_mean): Change the # of data points or the initial guess of the string.')
	else:
	        #Free Energy Calc       
		free_energy=[]
		for j in range(len(P2)):
                	exec "free_energy.append((-1)*kB*T*math.log(float(len(vor%s))/(len(P1))))" % (j)		

	P3 = [[0]*DIM]*(len(P2))
	for i in range(len(P2)):
		exec 'a=ff%s' % (i)
		P3[i]=a

        #find the closest point to the center of each voronoi cell
	#hhh=np.linalg.norm(P1[1]-P3[1])
	closeP=[]
	closeI=[]
        for j in range(len(P2)):
                for i in range(len(P1)):
                        dist=np.linalg.norm(P1[i]-P3[j])
			if i==0:
				minid=maxdistances[j]
                                index=0
			if dist < minid:
                                minid=dist
                                index=i
		closeP.append(P1[index])
		closeI.append(index)
	return P3, fail, free_energy, allpoints, closeP, closeI

##################################################################
def string_ini(P1,DIM,M):
#####################################
# Estimating an initial string (P2) #
#####################################
# This function first divides the total number of points in the data set (len(P1)) 
# by the number of points in the string obtaining 'S'. 
# Next we calculate M equidistant points starting from m1 (mean point of the first
# S points of P1) to m2 (mean point of the last S points of P1) 
# ###################################
# input:
# P1  = original data points
# DIM = dimension of the points
# M   = # of points of the vector P2 to be estimated
# output:
# P2  = initial string
##################################################################
	S=int(len(P1)/M)
	PA=P1[0:S]
	PA1=zip(*PA)
	m1=[]
	for i in range(DIM):
		m1.append(sum(PA1[i])/S)	
	PB=P1[len(P1)-S:len(P1)]
	PA2=zip(*PB)
	m2=[]
	for i in range(DIM):
        	m2.append(sum(PA2[i])/S)

	tmp=[]
	for i in range(0,DIM):
		tmp.append('')
		for i, value in enumerate(tmp):
			exec "H%s=[]" % (i)

	for i in range(DIM):
		step=(math.fabs(m1[i]-m2[i]))/M
		for j in range(M):
			if (m1[i]<m2[i]):
				coord=(m1[i])+(step*j)
			else:
				coord=(m1[i])-(step*j)
			exec 'H%s.append(coord)' % (i)
	GG=[]
	for i in range(DIM):
		exec 'GG.append(H%s)' % i

	P2=zip(*GG)
	return P2

######################################################################
def dm_vec(P3,Mout):
######################################################################
# calculate dm (mean value of the distances) between each i and i+1 
# point of vector P3 
#######################################################################
	P3A=np.asarray(P3)
	c=[]
	for i in range(0,len(P3A)-1):
        	dist=np.linalg.norm(P3A[i]-P3A[i+1])
        	c.append(dist)
	if Mout<3:
		print 'Warning: too few points for the reparametrization'
	#else:
	#	print 'mean distance calc: OK'
	total_dist=sum(c)	
	mean_dist=sum(c)/(Mout-1)
	return mean_dist, total_dist
###################################################################
def midpoint(a,b,DIM):
# Returns the midpoint between points a and b of dimension DIM
##################################################################
	s=[]
	for i in range(DIM):
        	tmp=((a[i]+b[i])/2)
        	s.append(tmp)
	return s

######################################
def reparam(Pin,DIM,Mout,meanD):
##########################################################################
# Reparametrizes the input vector Pin to a vector Plist with Mout points
# with distances equal to dist_to_X2 between them. DIM is the dimension of
# both input and output vectors.
###########################################################################
        Pout = [[0]*DIM]*Mout
        for i in range(Mout):
                if i==0:
                        Pout[i]=Pin[0]
                        X1=Pin[0]
                        index1=1
                elif i==(Mout-1):
                        Pout[i]=Pin[len(Pin)-1]
                else:
                        dist=np.linalg.norm(X1-Pin[index1])
                        if (dist>meanD):
                                X2, index2 = getX2near(X1,index1,meanD,Pin,DIM)
                        else:
                                X2, index2 = getX2far(X1,index1,meanD,Pin,DIM)

                        Pout[i]=X2
                        X1=X2
                        index1=index2
        return Pout


##############################
def smooth_DIM(Plist,s,DIM):
##############################################################################
#Smooth of the parametrize curve. The method followed correspond to eq. 48
#of Maragliano et al. J. Chem. Phys. 125, 024106 (2006)
#Input: Plist = reparametrized vector of dimension DIM
#       s     = factor that determines the strength of the smoothing           
###############################################################################
        Ps = [[0]*DIM]*(len(Plist))
        for i in range(len(Plist)):
                if i==0:
                        Ps[i]=Plist[0]
                elif i==(len(Plist)-1):
                        Ps[i]=Plist[len(Plist)-1]
                else:
                        Ps[i]=(((1-s)*Plist[i])+((s/2)*(Plist[i-1]+Plist[i+1])))
        return Ps

#####################################################################################################################
def savepcurve(N,DIM,Min,Mout,smooth,Niter,Ps,Psout,converg,free_energy,allpoints,closeP,closeI,closePP,closeII):  ##
#####################################################################################################################
##### Save selected variables as one "state" of the calculations ########
#########################################################################
# selected variables#
#####################
        allvars=[]
        allvars.append(N)
        allvars.append(DIM)
        allvars.append(Min)
        allvars.append(Mout)
        allvars.append(smooth)
        allvars.append(Niter)
        allvars.append(Ps)
        allvars.append(Psout)
        allvars.append(converg)
        allvars.append(free_energy)
        allvars.append(allpoints)
	allvars.append(closeP)
	allvars.append(closeI)
	allvars.append(closePP)
	allvars.append(closeII)
        # save to file pcurve_state.dat
        # note: the file that will use pcurve_state.dat should retrive the 
        # independent variables following the same order used here.
        fileOUT='pcurve_state.dat'
        file1 = open(fileOUT, 'w')
        pickle.dump(allvars, file1)

#######################################################################################
def loadpcurve(fileIN):  ##
#######################################################################################
##### Save selected variables as one "state" of the calculations ########
#########################################################################
# selected variables#
#####################
	file2 = open(fileIN, 'r')
	allvars = pickle.load(file2)
	N=allvars[0]
	DIM=allvars[1]
	Min=allvars[2]
	Mout=allvars[3]
	smooth=allvars[4]
	Niter=allvars[5]
	Ps=allvars[6]
	Psout=allvars[7]
	converg=allvars[8]
	free_energy=allvars[9]
	allpoints=allvars[10]
	closeP=allvars[11]
	closeI=allvars[12]
	closePP=allvars[13]
	closeII=allvars[14]
	return N, DIM, Min, Mout, smooth, Niter, Ps, Psout, converg, free_energy, allpoints, closeP, closeI, closePP, closeII


