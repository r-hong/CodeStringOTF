#! /usr/bin/env python
###################################################
# Rolando P. Hong Enriquez (rolando.hong@iit.it) ##
###################################################
import argparse
import pickle
import math
import random
from functions import load_data_range, delrn, getE
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

#### starting conditions
CV_file='rmsd.dat'
N=50 #number of points
CVini=1
CVend=100
L=CVend-CVini #initial dimensions
#########################################

P1=load_data_range(CV_file,N,CVini,CVend)
print P1.shape
print L
rn=(random.randint(0, L))
print rn
print P1

print '#####################'
#call delrn
P2=delrn(P1,rn)
print P2
#get the error
E=getE(P1,P2)
print E

EE=[]
gradient=1000
tolerance=10**(-25)
E1=1.00
c=0
while gradient>tolerance:
	rn=(random.randint(0, L))
	P2=delrn(P1,rn)
	E2=getE(P1,P2)
	if E2<E1:
		EE.append(E2)
		tmp=E1
		E1=E2
		P1=P2
		L=L-1
		gradient=(E-E2)**2
		print 'gradient', gradient
		c=0
	c=c+1
	if c==1000:
		print 'exit, we get traped'
		break	

print 'L', L
print 'EE', EE

#plot
EF=EE[2:len(EE)]
XF=[]
for i in range(1,len(EF)+1):
        XF.append(i)

pl.plot(XF,EF,c='r', linewidth=2)
pl.ylabel('Error')
pl.xlabel('Accepted Iterations')
pl.title('Monte Carlo')
pl.xlim(1,len(EF))
pl.ylim(0,max(EF)+0.1)
pl.fill_between(XF, EF, facecolor='blue', alpha=0.1)
pl.show()











