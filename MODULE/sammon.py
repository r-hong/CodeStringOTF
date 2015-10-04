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
#from collections import defaultdict, Counter

e=2.718281828459
X=5
Y=e**X

print Y

##
D1=-1.100*(10**6)
D2=-1.000*(10**6)

result=(D2-D1)
print result
#ff=e**result
#print 'ff', result

ww1=(math.exp((D1-D2)/0.1))
ww2=e**((D1-D2)/0.1)

print ww1,ww2

#### starting conditions
CV_file='rmsd.dat'
N=6 #number of points
CVini=1
CVend=200
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

###############################
#### Monte Carlo
EE=[]
c=[]
E1=1000000.0000
KT=0.01
L
MC=100 #Monte Carlo iterations
for i in range (MC):
	rn=(random.randint(0, L))
	P2=delrn(P1,rn)
	E2=getE(P1,P2)
	E2=E2*(10**6)
	rnMC=(random.randint(0, 100000))/100000.00000
	#print rnMC
	print math.exp((E1-E2)/KT) 
	if (math.exp((E1-E2)/KT) > rnMC):
		#c.append(i)
		P1=P2
		E1=E2
		L=L-1
		EE.append(E2)


print 'L', L
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











