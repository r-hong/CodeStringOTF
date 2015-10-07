#! /usr/bin/env python
"""
Secondary functions:
This module contains a series of functions that are called from all the other important scripts in the package.
"""

import csv
import math
import numpy as np
import sys
import argparse
import os
import pylab as pl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import random
import pickle
import pandas as pd

def normvect(L, normalizeTo=1):
	"""
	This function normalize the values of a list to make its max=normalizeTo
	"""
	vMax = max(L)
	if vMax<0.01:
		vMax=0.01
	return [ x/(vMax*1.0)*normalizeTo for x in L]

def diffstring(fileA,fileB,Ncv,Np):
	"""
	uses: function normvect (normalize vector) 
	This function loads data from files fileA and fileB containing 2 "strings" or iterations 
	of the string method. In these files:
	(1) each column represents a collective variable (Ncv = <N. of CVs>)
	(2) each row represent the value of the colective variables in one point of the path (Np = <N. of points in the path>)
	(3) columns are separated by a space character: ' ' 
	"""
        #### Initializing temp arrays ###
        tmp=[]
        for i in range(0,Ncv):
                tmp.append('')
        for i, value in enumerate(tmp):
                exec "cvA%s=[]" % (i)
                exec "cvB%s=[]" % (i)
                exec "NcvA%s=[]" % (i)
                exec "NcvB%s=[]" % (i)
                exec "diff%s=[]" % (i)
                exec "suma%s=[]" % (i)

        #### reading and normalizing CVs for fileA ########
        with open(fileA) as file:
                reader = csv.reader(file, delimiter=' ',skipinitialspace=True)
                for row in reader:
                        for i in range(0,Ncv):
                                exec "cvA%s.append(float(row[i]))" % (i)
                                #exec "NcvA%s=normvect(cvA%s)" % (i,i)

        #### reading and normalizing CVs for fileB ########
        with open(fileB) as file:
                reader = csv.reader(file, delimiter=' ',skipinitialspace=True)
                for row in reader:
                        for i in range(0,Ncv):
                                exec "cvB%s.append(float(row[i]))" % (i)
                                #exec "NcvB%s=normvect(cvB%s)" % (i,i)

        # Calculating the differences between each CV in fileA and fileB 
        for i in range(0,Ncv):
                for j in range(0,Np):
                        exec "diff%s.append(cvA%s[j]-cvB%s[j])" % (i,i,i)
                        exec "suma%s=abs(float(sum(diff%s)/%s))" % (i,i,Np)

        # adding all the differences to obtain a total difference between fileA and fileB (i.e., between stringA and stringB)
        total1=[]
	total2=[]
        for i in range(0,Ncv):
                #exec "total.append(suma%s+suma%s)" % (i,i)
		exec "total1.append(suma%s)" % (i)
                #break
	total2=sum(total1)
        return total2

def parse_conf(conf_file):
	"""
	Parse the parameter of the configuration file. 
	Information on these parameters can be found in the example provided (input_conf.dat)   
	usage: 
	POINT_INI, POINT_END, ITER, NAMD_EXE, NAMD_CONF_CORE, COOR_FILE, PSF_FILE, PARAM_FILE, MDSTEPS, TIME_STEP_AKMA, GAMMA, KIN, TCL_FILE_CORE, NUM_COLVAR = parse_conf(conf_file)
	"""
	#### Initializations #####
	CONF_FILE=conf_file   #'input_conf.dat'
	Nparm=14 #number of initial parameters in the input_conf.dat file
	tmp=[]
	for i in range(0,Nparm):
		tmp.append('')
	for i, value in enumerate(tmp):
		exec "FF%s=[]" % (i)

	#### Open file for reading ####
	file = open(CONF_FILE)
	c=0
	for line in file:
		exec "FF%s = line.strip().split()" % (c)  
		c=c+1

	# error handling for the parameters in the file (not elegant solution)
	if FF0[0]=='POINT_INI':
		POINT_INI=int(FF0[1])
	else:
		print 'Error!! POINT_INI not set or in the wrong position within the file'

	if FF1[0]=='POINT_END':
        	POINT_END=int(FF1[1])
	else:
        	print 'Error!! POINT_END not set or in the wrong position within the file'

	if FF2[0]=='ITER':
        	ITER=int(FF2[1])
	else:
        	print 'Error!! ITER not set or in the wrong position within the file'

	if FF3[0]=='NAMD_EXE':
        	NAMD_EXE=FF3[1]
	else:
        	print 'Error!! NAMD_EXE not set or in the wrong position within the file'

	if FF4[0]=='NAMD_CONF_CORE':
        	NAMD_CONF_CORE=FF4[1]
	else:
        	print 'Error!! NAMD_CONF_CORE not set or in the wrong position within the file'

	if FF5[0]=='COOR_FILE':
        	COOR_FILE=FF5[1]
	else:
        	print 'Error!! COOR_FILE not set or in the wrong position within the file'

	if FF6[0]=='PSF_FILE':
        	PSF_FILE=FF6[1]
	else:
        	print 'Error!! PSF_FILE not set or in the wrong position within the file'

	if FF7[0]=='PARAM_FILE':
        	PARAM_FILE=FF7[1]
	else:
        	print 'Error!! PARAM_FILE not set or in the wrong position within the file'

	if FF8[0]=='MDSTEPS':
        	MDSTEPS=int(FF8[1])
	else:
        	print 'Error!! MDSTEPS not set or in the wrong position within the file'

	if FF9[0]=='TIME_STEP_AKMA':
        	TIME_STEP_AKMA=float(FF9[1])
	else:
        	print 'Error!! TIME_STEP_AKMA not set or in the wrong position within the file'

	if FF10[0]=='GAMMA':
        	GAMMA=float(FF10[1])
	else:
        	print 'Error!! GAMMA not set or in the wrong position within the file'

	if FF11[0]=='KIN':
        	KIN=int(FF11[1])
	else:
        	print 'Error!! KIN not set or in the wrong position within the file'

	if FF12[0]=='TCL_FILE_CORE':
        	TCL_FILE_CORE=FF12[1]
	else:
        	print 'Error!! TCL_FILE_CORE not set or in the wrong position within the file'
	if FF13[0]=='NUM_COLVAR':
        	NUM_COLVAR=int(FF13[1])
	else:
        	print 'Error!! NUM_COLVAR not set or in the wrong position within the file'
	return POINT_INI, POINT_END, ITER, NAMD_EXE, NAMD_CONF_CORE, COOR_FILE, PSF_FILE, PARAM_FILE, MDSTEPS, TIME_STEP_AKMA, GAMMA, KIN, TCL_FILE_CORE, NUM_COLVAR

def get_cv_info(POINT_END,POINT_INI,CV_FILE):
	"""
	Getting type of CV and initial and final lines for each CV in the input file 
	Format of the input file:
	CV_NAME FORCE_CONSTANT
	LINE_A
	LINE_B
	LINE_1
	LINE_2
	LINE_3
	*
	
	EXPLANATION
	----------- 
	(1) In the fist line there is the name of the CV and the force constant to be applied to this CV 
	Example: DIST 100
	(2) LINE_A, LINE_B, ...etc. are the lines defining the indexes of the atoms or group of atoms for the CV
	Example for the CV DIST, there should be only LINE A and LINE B, on each line there will one or several atoms.
	is several atoms are entered, then DIST will be calculated using the center of mass for that group.
	(3) LINE_1, LINE_2, LINE_3,...LINE_NUN_POINTS, are the initial values of the CV on each point along the
	path between the initial and final state of the system. This path represents the first 'string' that will be iterated
	in the method to find the minimum free energy path.
	(4) The end of the CV definition should be marked with a line with the character '*'. 
	The next CV should be defined immediately after, without leaving empty lines.
	"""
	NUM_POINTS=(POINT_END - POINT_INI)+1
	#######################################################
	### Note: add reference to a new CV in this list ######
	DEFINED_COLVARS=['DIH', 'DIST'] #######################
	#######################################################

	#### getting number of lines in the file CV_FILE #####
	NUM_LINES=0
	file = open(CV_FILE)
	for line in file:
		NUM_LINES=NUM_LINES+1

	#### Initializing empty lists #####
	tmp=[]
	for i in range(0,NUM_LINES):
        	tmp.append('')
        	for i, value in enumerate(tmp):
                	exec "FF%s=[]" % (i)
                	exec "dd=[]"
			exec "cvDIH=[]"
			exec "cvDIST=[]"
			###############################################
			### Note: initialize empty list foranew CV #### 
 			###############################################

	#### Getting type of CV and initial and final lines for each CV in the input file ####
	c=zero=cv_ini=cv_end=cv_current=0
	file = open(CV_FILE)
	for line in file:
		dd=[]
		exec "FF%s = line.strip().split()" % (c)
		exec "dd.append(FF{}[{}])".format(c,zero)
		c=c+1
		if dd[0] in DEFINED_COLVARS:
			cv_ini=c
			cv_current=cv_current+1
			#print 'CV No. {} --- 1rst line: {}'.format(cv_current,cv_ini)
			if dd[0]=='DIH':
				cv_end=cv_ini+NUM_POINTS+4
				#print 'CV No. {} --- last line: {}'.format(cv_current,cv_end)
				#print 'call dihedral parser'
				exec "FF%s = line.strip().split()" % (c)
				exec "cvDIH.append(FF{}[{}])".format(c,0)
				cvDIH.append(cv_ini)
				cvDIH.append(cv_end)
			elif dd[0]=='DIST':
				cv_end=cv_ini+NUM_POINTS+2
				#print 'CV No. {} --- last line: {}'.format(cv_current,cv_end)
				#print 'call distance parser'
                        	exec "FF%s = line.strip().split()" % (c)
                        	exec "cvDIST.append(FF{}[{}])".format(c,zero)
                        	cvDIST.append(cv_ini)
                        	cvDIST.append(cv_end)
		################################################################
		# Note: add reference to a new CV in the 'if elif' conditions ##
		################################################################
	return NUM_LINES, cvDIH, cvDIST 


def diffpoint(fileA,fileB,Ncv,Np,a):
	"""
	uses: function normvect (normalize vector) 
	This function loads data from files fileA and fileB containing 2 "strings" or iterations 
	of the string method. In these files:
	(1) each column represents a collective variable (Ncv = <N. of CVs>)
	(2) each row represent the value of the colective variables in one point of the path (Np = <N. of points in the path>)
	(3) columns are separated by a space character: ' ' 
	"""
        #### Initializing temp arrays ###
        tmp=[]
        for i in range(0,Ncv):
                tmp.append('')
        for i, value in enumerate(tmp):
                exec "cvA%s=[]" % (i)
                exec "cvB%s=[]" % (i)

        #### reading and normalizing CVs for fileA ########
        with open(fileA) as file:
                reader = csv.reader(file, delimiter=' ',skipinitialspace=True)
		c=0
                for row in reader:
			c=c+1
			if (c==a):
                        	for i in range(0,Ncv):
                                	exec "cvA%s.append(float(row[i]))" % (i)
        #### reading and normalizing CVs for fileB ########
        with open(fileB) as file:
                reader = csv.reader(file, delimiter=' ',skipinitialspace=True)
		c=0
                for row in reader:
			c=c+1
			if (c==a):	
                        	for i in range(0,Ncv):
                                	exec "cvB%s.append(float(row[i]))" % (i)
	# calculating the distance between the image(row) 'a' in CV space 
	# extracted from fileA and fileB
	suma=0
	for j in range(Ncv):
		exec "AA=np.asarray(cvA%s)" % (j)
		exec "BB=np.asarray(cvB%s)" % (j)
		t1=(AA-BB)**2
		suma=suma+t1
	dist_i=math.sqrt(suma)
	dd=[]
	dd.append(dist_i)
	return dd


def file_lines(fname):
	"""
	Returns the number of lines in a file
	"""
        with open(fname) as f:
                for i, l in enumerate(f):
                        pass
        return i + 1


def file_cols(fname):
	"""
	Returns the number of columns in a file
	"""
        with open(fname) as f:
                reader = csv.reader(f, delimiter=' ')
                first_row = next(reader)
                num_cols = len(first_row)
        return num_cols

def readFile(fileName,nCol):
        """
        This function reads a file 'fileName' with 'cCol' columns into a matrix. 
	This code is used for instance in the code of mforce.py calculate 
	the maximum likehood of the forces in an arbitrary number of coordinates.
        """
        #### Initializing temp arrays ###
        tmp=[]
        for i in range(0,nCol):
                tmp.append('')
        for i, value in enumerate(tmp):
                exec "cv%s=[]" % (i)

        #### reading the columns for file fileName ########
        with open(fileName) as file:
                reader = csv.reader(file, delimiter=' ')
                for row in reader:
                        for i in range(0,nCol):
                                exec "cv%s.append(float(row[i]))" % (i)
        vecT=[]
        for i in range(0,nCol):
                        exec "vecT.append(cv%s)" % (i)
        vecT=np.asarray(vecT)
        return vecT



def delrn(A,index):
	"""
	This function returns an array in which the collective variable with index 'index' has been eliminated from the array 'A'
	"""
        P1T=zip(*A)
        S1=np.delete(P1T,index,0)
        S2=zip(*S1)
        P2=np.asarray(S2)
        return P2

def getE(P1,P2):
	"""
	Function to estimate the error of using Sammon mapping (dimensionality reduction), 
	that is, to make a projection from a vector of higher dimension P1 
	to the one of lower dimmension P2 
	"""
        ST=[]
        SD=[]
        for i in range(len(P1)):
                for j in range(i+1,len(P1)):
                        d1=np.linalg.norm(P1[i]-P1[j])
                        d2=np.linalg.norm(P2[i]-P2[j])
                        S1=((d1-d2)**2)/d1
                        ST.append(S1)
                        SD.append(d1)
        E=sum(ST)/sum(SD)
        return E


def load_data_range(filename,N,CVini,CVend):
	"""
	This function loads a range of data. That is, it loads N points (rows) from columns CVini to CVend in a data file 'filename'
	"""
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

def getX2near(X1,index1,dist_to_X2,Pin,DIM):
	"""
	Find point X2 lying in the line between the DIM-dimensional points X1 and the point Pin[index1]
	located at distance dist_to_X2 from X1.  
	"""
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
        """
        Find point X2 lying in the line between the DIM-dimensional points X1 and the point Pin[index1]
        located at distance larger than dist_to_X2 from X1.  
        """
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
                return X2, index2


def load_data(filename,N,DIM):
	"""
	load N points of dimension DIM. 
	data should be in a file 'filename' with N rows and DIM columns
	"""
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

def maxvec(P1,P2):
	"""
	Calculates the maximal distance among all points of vector P1 
	to each point of vector P2. Return the vector of maximal distances
	of size len(P2).
	"""
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


def voroid_mean(P1,P2,DIM,maxdistances,fixed):
	"""
	Obtaining the mean values of each voronoid cell 
	-----------------------------------------------
	input:
	DIM = dimension of the points
	maxdistances = vector of maximal distances calculated with the function maxvec
	P1  = original data
	P2  = string vector, for each point of P2 (P2[i])we calculate a voronoid cell:
	the points of P1 closer to P2[i] than to any other point of P2.
	fixed = [0,1]. If fixed=1 then the first and last voronoi centers are kept the same as in the 
	original string P2. Any other value for 'fixed' will find all the voronoi centers
	output:
	P3= mean points of each voronoid cell point (each point with dimension DIM)
	fail = binary number. fail=1 means that there were empty voronoid cells, and the procedure
	should be restarted with other conditions(e.g., more points on P1 or a different P2)
	"""	
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
	if fixed==1:		
		for i in range(len(P2)):
			exec 'a=ff%s' % (i)
			if i==(len(P2)-1):
				P3[i]=P2[len(P2)-1]
			elif i==0:
				P3[i]=P2[0]
			else:
				P3[i]=a
	else:
                for i in range(len(P2)):
                        exec 'a=ff%s' % (i)
                        P3[i]=a

        #find the closest point to the center of each voronoi cell
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

def string_ini(P1,DIM,M):
	"""
	This function generates an initial string for the calculation of the principal curve 
	---------------------------------
	The function first divides the total number of points in the data set (len(P1)) 
	by the number of points in the string obtaining 'S'. 
	Next we calculate M equidistant points starting from m1 (mean point of the first
	S points of P1) to m2 (mean point of the last S points of P1) 
	------------------------------------------------------------
	input:
	P1  = original data points
	DIM = dimension of the points
	M   = # of points of the vector P2 to be estimated
	output:
	P2  = initial string
	"""
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


def string_centers(m1,m2,DIM,M):
	"""
	This function generates an initial string for the calculation of the principal curve
	------------------------------------------------------------------------------------
	This function uses the center of two clusters (c1 and c2) as the extreme values to find an initial string. 
	A final string of M points is estimated by finding equidistant points between c1 and c2. 
	----------------------------------------------------------------------------------------
	input:
	c1  = center of the first cluster
	c2  = center of the second cluster	
	DIM = dimension of the points
	M   = # of points of the vector P2 to be estimated
	output:
	P2  = initial string
	"""
        tmp=[]
        for i in range(0,DIM):
                tmp.append('')
                for i, value in enumerate(tmp):
                        exec "H%s=[]" % (i)

        for i in range(DIM):
                step=(math.fabs(m1[i]-m2[i]))/(M-1)
                for j in range(M):
			#if j==M:
			#	coord=m2[i]
                        if (m1[i]<m2[i]):
                                coord=(m1[i])+(step*j)
                        else:
                                coord=(m1[i])-(step*j)
                        exec 'H%s.append(coord)' % (i)
        GG=[]
        for i in range(DIM):
                exec 'GG.append(H%s)' % i

        P2=zip(*GG)
        return GG

def dm_vec(P3,Mout):
	"""
	This function calculate dm (mean value of the distances) between each i and i+1 
	point of vector P3 
	"""
	P3A=np.asarray(P3)
	c=[]
	for i in range(0,len(P3A)-1):
        	dist=np.linalg.norm(P3A[i]-P3A[i+1])
        	c.append(dist)
	if Mout<3:
		print 'Warning: too few points for the reparametrization'
	total_dist=sum(c)	
	mean_dist=sum(c)/(Mout-1)
	return mean_dist, total_dist

def midpoint(a,b,DIM):
	"""
	This function returns the midpoint between points a and b of dimension DIM
	"""
	s=[]
	for i in range(DIM):
        	tmp=((a[i]+b[i])/2)
        	s.append(tmp)
	return s

def reparam(Pin,DIM,Mout,meanD):
	"""
	Reparametrizes the input vector Pin to a vector Plist with Mout points
	with distances equal to dist_to_X2 between them. DIM is the dimension of
	both input and output vectors.
	"""
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

def smooth_DIM(Plist,s,DIM):
	"""
	This function performs a smoothing of the parametrize curve. The method followed correspond to eq. 48
	of Maragliano et al. J. Chem. Phys. 125, 024106 (2006)
	-----------------------------------------------------
	Input: Plist = reparametrized vector of dimension DIM
	       s     = factor that determines the strength of the smoothing
	       DIM   = dimension of the curve	           
	"""
        Ps = [[0]*DIM]*(len(Plist))
        for i in range(len(Plist)):
                if i==0:
                        Ps[i]=Plist[0]
                elif i==(len(Plist)-1):
                        Ps[i]=Plist[len(Plist)-1]
                else:
                        Ps[i]=(((1-s)*Plist[i])+((s/2)*(Plist[i-1]+Plist[i+1])))
        return Ps

def savepcurve(N,DIM,Min,Mout,smooth,Niter,Ps,Psout,converg,free_energy,allpoints,closeP,closeI,closePP,closeII):
	"""
	This function saves the selected variables as one "state" of the calculations of the principal curve. 
	The saved variables can be loaded into the workspace using the function 'loadpcurve'
	"""
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
        fileOUT='pcurve_state.dat'
        file1 = open(fileOUT, 'w')
        pickle.dump(allvars, file1)

def loadpcurve(fileIN):  
	"""
	This function loads the selected variables as one "state" of the calculations of the principal curve, 
	specifically the variables must be saved with the function 'savepcurve'
	"""
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

def mergevor(allpoints,P1,a,b):
	"""
	This function merges the voronoid cells 'a' and 'b'. 'allpoints' the indexes of all the points of the voronoi cells   
	"""
        p0=allpoints[a]
        p1=allpoints[b]
        index=[]
        for i in range(len(p0)):
                ind=p0[i]
                index.append(ind)
        for i in range(len(p1)):
                ind=p1[i]
                index.append(ind)

        Pj=[]
        for i in range(len(index)):
                point=P1[index[i]]
                Pj.append(point)
        return index,Pj

def binom(n, m):
	"""
	This function returns all the combinations of n in m
	"""
	b = [0] * (n + 1)
	b[0] = 1
	for i in xrange(1, n + 1):
		b[i] = 1
		j = i - 1
		while j > 0:
			b[j] += b[j - 1]
			j -= 1
	return b[m]

def pdb2df(input_pdb):
	"""
	Load file in pdb format into a data frame using 'pandas'
	"""
	l1=1
	l2=file_lines(input_pdb)
	col1=1
	col2=2
	col3=3
	col4=4
	col5=5
	col6=6
	col7=7
	col8=8
	col9=9
	col10=10
	file_tmp='tmp.pdb'
	a1='sed -n '
	a1+=`l1`
	a1+=','
	a1+=`l2`
	a1+='p '
	a1+=input_pdb
	a1+=' | awk -v OFS='
	a1+="'"
	a1+="  "
	a1+="' '{print $"
	a1+=`col1`
	a1+=', $'
	a1+=`col2`
	a1+=', $'
	a1+=`col3`
	a1+=', $'
	a1+=`col4`
	a1+=', $'
	a1+=`col5`
	a1+=', $'
	a1+=`col6`
	a1+=', $'
	a1+=`col7`
	a1+=', $'
	a1+=`col8`
	a1+=', $'
	a1+=`col9`
	a1+=', $'
	a1+=`col10`
	a1+='}'
	a1+="'"
	a1+=' > '
	a1+=file_tmp
	os.system(a1)

	#### reading the coordinates from file_tmp ########
	c1=[]
	c2=[]
	c3=[]
	c4=[]
	c5=[]
	c6=[]
	c7=[]
	c8=[]
	c9=[]
	c10=[]
	with open(file_tmp) as file:
        	reader = csv.reader(file, delimiter=' ')
        	for row in reader:
                	c1.append(row[0])
                	c2.append(row[2])
                	c3.append(row[4])
                	c4.append(row[6])
                	c5.append(row[8])
                	c6.append(row[10])
                	c7.append(row[12])
                	c8.append(row[14])
                	c9.append(row[16])
                	c10.append(row[18])

	#dictionary for the pdb dataframe
	d={'atom': c1,
	   'atomid': c2,
	   'atomname': c3,
	   'resname': c4,
	   'resid': c5,
	   'x': c6,
	   'y': c7,
	   'z': c8,
	   'occupancy': c9,
	   'beta': c10}

	df = pd.DataFrame(d)
	return df

def df2pdb(inDF,outPDB):
	"""
	Save a data frame with molecular information (loaded using pandas) to a pdb file format
	---------------------------------------------------------------------------------------
	input:  inDF   = pdf file in data frame format (using pandas)
	output: outPDB = a file in pdb format
	""" 
        c1=inDF['atom']
        c2=inDF['atomid']
        c3=inDF['atomname']
        c4=inDF['resname']
        c5=inDF['resid']
        c6=inDF['x']
        c7=inDF['y']
        c8=inDF['z']
        c9=inDF['occupancy']
        c10=inDF['beta']

        f = open(outPDB, "w")
        for i in xrange(len(c1)):
                f.write("{:<4}{:>7} {:^4} {:<3}{:>6}{:>12}{:>8}{:>8}{:>6}{:>6}\n".format(c1[i], c2[i], c3[i], c4[i], c5[i], c6[i], c7[i], c8[i], c9[i],c10[i]))
        f.close()


def dist_atom(dfA,dfB,index):
	"""
	Given the pdb dataframes dfA and dfB calculate the distance between 
	the equivalent atoms of index equal to 'index'
	"""
        dist=np.sqrt((float(dfA['x'].iloc[index])-float(dfB['x'].iloc[index]))**2 + (float(dfA['y'].iloc[index])-float(dfB['y'].iloc[index]))**2 + (float(dfA['z'].iloc[index])-float(dfB['z'].iloc[index]))**2)
        return dist



def get_ligand_coord(filename,l1,l2,x,y,z):
	"""
	This function get the coordinates of a ligand.
	----------------------------------------------
	Input:
	filename = name of the pdb file
	l1,l2    = first and last lines on the pdb file containing the ligand
	x,y,z    = columns on the pdb file containing the coordinates
	Output:
	x,y,z,   = vectors containing the coordinates of every atom in the ligand 
	"""
        # creating the string with the external command. The command uses sed and awk to 
        # pre-process the initial pdb file and creating  a temporal file 
        # with the coordinates of the ligand.   
        file_tmp='tmp.pdb'
        a1='sed -n '
        a1+=`l1`
        a1+=','
        a1+=`l2`
        a1+='p '
        a1+=filename
        a1+=' | awk -v OFS='
        a1+="'"
        a1+="  "
        a1+="' '{print $"
        a1+=`x`
        a1+=', $'
        a1+=`y`
        a1+=', $'
        a1+=`z`
        a1+='}'
        a1+="'"
        a1+=' > '
        a1+=file_tmp
        os.system(a1)
        #### reading the coordinates from file_tmp ########
        x1=[]
        y1=[]
        z1=[]
        with open(file_tmp) as file:
                reader = csv.reader(file, delimiter=' ')
                for row in reader:
                        x1.append(row[0])
                        y1.append(row[2])
                        z1.append(row[4])

        x2=np.asarray(x1)
        y2=np.asarray(y1)
        z2=np.asarray(z1)
        xx=[]
        yy=[]
        zz=[]
        for i in range(len(x1)):
                xt=float(x2[i])
                yt=float(y2[i])
                zt=float(z2[i])
                xx.append(xt)
                yy.append(yt)
                zz.append(zt)
        return xx,yy,zz

def rotation_translation_3D(S1Ref, S1Target, S2Ref):   ###
	"""
	########   Application 1 ---> Ligand/Receptor system   ###################
	Implement ridig body transformations (rotations/translations) over a set of points following a reference
	Although this is a gereral problem in this context the application is inspired by the following problem:
	giving a receptor/ligand system, if we know the coordinates of the system and time1(reference) and the coordinates
	of the receptor at time2(target), estimate the coordinates of the ligand at time2(target).   
	#######   Application 2 ---> measuring conformational change #############
	We can have a ligand in:
	(1) ligand Reference position&conformation 
	(2) ligand Target    position&conformation, 
	if we simply set as the new reference, the same ligand reference (S2Ref=S1Ref), then the resulting output target 
	(S2Target) will contain the conformation of the ligand in the target and the position of the ligand in the reference,
	this is: 
	(3)  ligand ReferencePosition & TargetConformation
	 Therefore, the difference between (1) and (3) (i.e., between S2Target and S1Ref) if the ligand's conformational change
	--------------------------------------------------------------------------------------------------------
	 Input:
		S1Ref    = set of points 1 (in the reference system)
		S1Target = set of points 1 (in the target system)
		S2Ref    = set of points 2 (in the reference system)
	 Output:
		S2Target = set of points 2 (in the target system)
	"""
	assert ((len(S1Ref) == len(S1Target)) & (len(S1Target)==len(S2Ref))), "the input data sets should have equal length!!" 
	N = S1Ref.shape[0]; # getting the number of points
	centerS1Ref    = mean(S1Ref, axis=0)
	centerS1Target = mean(S1Target, axis=0)
    
	# Put the data sets at the origin of the coordinate system
	A = S1Ref    - tile(centerS1Ref,    (N, 1))
	B = S1Target - tile(centerS1Target, (N, 1))

	# get the covariance matrix covarM
	covarM = transpose(A) * B
	M1, M2, M3 = linalg.svd(covarM)
	Rotation = M3.T * M1.T
	# catch possible reflextion
	if linalg.det(Rotation) < 0:
		print "Special case: Reflection detected"
		M3[2,:] *= -1
		Rotation = M3.T * M1.T
	Translation = -Rotation*centerS1Ref.T + centerS1Target.T
	S2Target = (Rotation*S2Ref.T) + tile(Translation, (1, N))
	S2Target = S2Target.T
	return S2Target


def loadrmsd(fileIN):  
	"""
	Load rmsd data from a 'state' file
	"""
        file2 = open(fileIN, 'r')
        allvars = pickle.load(file2)
        NP=allvars[0]
        rmsd=allvars[1]
        X=allvars[2]
        return NP, rmsd, X

def savermsd(NP,rmsd,X,fileOUT):  
	"""
	Save rmsd data to a 'state' file
	"""
        allvars=[]
        allvars.append(NP)
        allvars.append(rmsd)
	allvars.append(X)
        file1 = open(fileOUT, 'w')
        pickle.dump(allvars, file1)

def bootstrapN(X,N):
	"""
	Bootstrap resampling.
	Uses as input a vector or arbitrary size and returns another vector of size N 
	composed by elements that are sampled (with replacement) from the input vector. 
	"""
        XresampledN=[]
        for i in range(0,N):
                 XresampledN.append(X[random.randrange(0, len(X))])
        return XresampledN

def readFile(fileName,nCol):
        '''
        Read data from file 'fileName' with 'nCol' columns into a matrix
        '''
        #### Initializing temp arrays ###
        tmp=[]
        for i in range(0,nCol):
                tmp.append('')
        for i, value in enumerate(tmp):
                exec "cv%s=[]" % (i)

        #### reading and normalizing CVs for fileA ########
        with open(fileName) as file:
                reader = csv.reader(file, delimiter=' ')
                for row in reader:
                        for i in range(0,nCol):
                                exec "cv%s.append(float(row[i]))" % (i)
        vecT=[]
        for i in range(0,nCol):
                        exec "vecT.append(cv%s)" % (i)
        vecT=np.asarray(vecT)
        return vecT

def estValB(data,K,N):
        '''
        This function finds the estimated value of each row in a matrix 'data' using bootstrap.
        Here:
        (1) data -- is the data matrix. len(data) should represent the number of columns of the matrix.
            finding the max likelihood value for each column is the goal of this function.
        (2) N    -- is the size of the resample bootstraped subvector taken from the values of each column.
        (3) K    -- is the number of times we will repeat the resampling of size N for each column.

        the function returns a list of size equal to the number of columns in 'data' in which each element
        is the maximum likelihood value for the specified column.       
        '''
        nCol=len(data)
        #### Initializing temp arrays ###
        tmp=[]
        for i in range(0,nCol):
                tmp.append('')
        for i, value in enumerate(tmp):
                exec "B%s=[]" % (i)
                exec "c%s=data[%s]" % (i,i)
        for j in range(1,K):
                tmp=[]
                for i in range(0,nCol):
                        tmp.append('')
                for i, value in enumerate(tmp):
                        exec "A%s=[]" % (i)
                for i in range(1,N):
                        for g in range(0,nCol):
                                exec "tmp%s=bootstrapN(c%s,N)" % (g,g)
                                exec "A%s.append(sum(tmp%s[0:%s])/len(tmp%s[0:%s])   )" % (g,g,i,g,i)
                for h in range(0,nCol):
                        exec "B%s.append(A%s[N-2])" % (h,h)
                allCalc=[]
                for p in range(0,nCol):
                        exec "allCalc.append(B%s)"% (p)
        return allCalc

def formatRes(myList,n,noFloat=False):
        '''
        This function formats the numeric elements of 'myList' to number with 'n' decimal values
        '''
        nCol=len(myList)
        tmp=[]
        formated=[]
        for i in range(0,nCol):
                tmp.append('')
        for i, value in enumerate(tmp):
                exec 'f%s=("{0:.%sf}".format(np.mean(myList[%s])))' % (i,n,i)
		if noFloat==False:
                	exec "formated.append(float(f%s))" % (i)
		else:
			exec "formated.append(int(f%s))" % (i)
        return formated


def cleanpdb(input_pdb):
	'''
	This is a custom function to clean a file pdb written by namd, that is, 
	eliminate the ions, remarks end end labels.
	This is a preprocess step before loading the pdb as a dataframe with using pandas.
	'''
        phrase1='REMARK'
        phrase2='END'
        phrase3='WAT'
        phrase4='Na+'
        f = open(input_pdb,"r")
        lines = f.readlines()
        f.close()

        f = open(input_pdb,"w")
        for line in lines:
                if (phrase1 in line) or (phrase2 in line) or (phrase3 in line) or (phrase4 in line):
                        a=1
                else:
                        f.write(line)


def ligCOM(df,ligand):   
	'''
	This function returns the center of mass (COM) of a ligand in a coordinate pdb file. 
	The pdb file has to be previously loaded into a data frame using pandas.
	input: (1) dataframe of the prot-ligand pdb
	       (2) name of the ligand (as in resname field of the pdb)
	output:(3) coordinates of the center of mass of the ligand  
	'''
        c=0
        x=0.0
        y=0.0
        z=0.0
        for i in range(len(df['beta'])):
                if (df['resname'].iloc[i]==ligand):
                        x=x+float(df['x'].loc[i])
                        y=y+float(df['y'].loc[i])
                        z=z+float(df['z'].loc[i])
                        c=c+1
        xc=x/c
        yc=y/c
        zc=z/c
        return xc, yc, zc


def getPdbFilename(a,b):  
	'''
	This is a custom function with very restricted use for the filenames we use in the original simulations.
	The function simply creates a file named namdout{a}-{b}.coor
	input: 
	   (1) a = index of the image in the string
	   (2) b = index of the iterationin the string method procedure
	output:
	   (3)  name of the pdb file
	'''
        ffname='namdout'
        ffname+=`a`
        ffname+='-'
        ffname+=`b`
        ffname+='.coor'
        return ffname

def getMeanString(allstrings): 
	'''
	This function calculates the mean string 
	input: 
	      allstrings = containts all the strings in the group
	output:
	      xm,ym,zm = vectors with the mean values (per image) of the coordinates for the mean string
	'''
        #initialize vector to save the sum of the coordinates for each string
        st=zip(*allstrings[0])
        x=st[0]
        Nimages=len(x)
        xt=[]
        yt=[]
        zt=[]
        for i in range(Nimages):
                xt.append(0.0)
                yt.append(0.0)
                zt.append(0.0)
        #calculate and save the sum of the coordinates for each string
        for i in range(len(allstrings)):
                st=zip(*allstrings[i])
                x=st[0]
                y=st[1]
                z=st[2]
                for j in range(Nimages):
                        xt[j]=xt[j]+x[j]
                        yt[j]=yt[j]+y[j]
                        zt[j]=zt[j]+z[j]
        #calculate the mean values of the coordinates for each string (the mean string)
        xm=[]
        ym=[]
        zm=[]
        for i in range(Nimages):
                xm.append(xt[i]/len(allstrings))
                ym.append(yt[i]/len(allstrings))
                zm.append(zt[i]/len(allstrings))
        return xm,ym,zm

