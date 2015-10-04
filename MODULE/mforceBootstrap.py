#! /usr/bin/env python
###################################################
# Rolando P. Hong Enriquez (rolando.hong@iit.it) ##
###################################################
import argparse
import os
import csv
import pylab as pl
import numpy as np
import math
import random
import scipy as sp
import scipy.stats
#######################################################
###### Parsing arguments from the command line ########
#######################################################
parser = argparse.ArgumentParser(description="Calculate the mean forces in time.",epilog="Thanks for using mforces.py!!")
parser.add_argument('IN_file', help="Input file")
parser.add_argument('N', type=int, help="Stride. Use only every Nth value from the input file.")
args = parser.parse_args()
###############################

#preprocessing the input file:
# saving every args.stride value of args.IN_file and saving to tmp.dat
ss="'"
ss+='0~'
ss+=`args.N`
ss+='p'
ss+="'"

cmd='sed -n '
cmd+=ss
cmd+=' '
cmd+=args.IN_file
cmd+=' > tmp.dat'
os.system(cmd)
##########################################

#### reading the file with forces (Fx,Fy,Fz) ########
c1=[]
c2=[]
c3=[]
with open('tmp.dat') as file:
	reader = csv.reader(file, delimiter=' ')
	for row in reader:
		c1.append(float(row[0]))
		c2.append(float(row[1]))
		c3.append(float(row[2]))



#reagange vectors x, y z for the calculation of the cosine similarity
vecT=[]
vecT.append(c1)
vecT.append(c2)
vecT.append(c3)
vecT=np.asarray(vecT)
vecF=vecT.T

#function for cosine similarity between 2 vectors
def cosine_similarity(v1,v2):
    "compute cosine similarity of v1 to v2: (v1 dot v1)/{||v1||*||v2||)"
    sumxx, sumxy, sumyy = 0, 0, 0
    for i in range(len(v1)):
        x = v1[i]; y = v2[i]
        sumxx += x*x
        sumyy += y*y
        sumxy += x*y
    return sumxy/math.sqrt(sumxx*sumyy)
#########################################
### bootstrap resampling function
########################################
def bootstrapN(X,N):
# returned a vector of size n that is a sample with replacement (boostrap)
	XresampledN=[]
	for i in range(0,N):
		 XresampledN.append(X[random.randrange(0, len(X))])	
	return XresampledN

#calculating the cummulative mean values K times (one for each bootstraped resampled vector of size N)
K=1000
N=100 #size of the resampled bootstraped vectors
Bx=[]
By=[]
Bz=[]
for j in range(1,K):
        aveBx=[]
        aveBy=[]
        aveBz=[]
        for i in range(1,N):
                tmpX=bootstrapN(c1, N)
		tmpY=bootstrapN(c2, N)
		tmpZ=bootstrapN(c3, N)
                aveBx.append(sum(tmpX[0:i])/len(tmpX[0:i]))
		aveBy.append(sum(tmpY[0:i])/len(tmpY[0:i]))
		aveBz.append(sum(tmpZ[0:i])/len(tmpZ[0:i]))
        Bx.append(aveBx[N-2])
	By.append(aveBy[N-2])
	Bz.append(aveBz[N-2])

#Calculatin confidence intervals for the bootstrap mean value stimation
confidence=0.95
#mean, standard error of the mean, standard deviations
mx, s1x, s2x = np.mean(np.array(Bx)), scipy.stats.sem(np.array(Bx)), np.std(np.array(Bx))
my, s1y, s2y = np.mean(np.array(By)), scipy.stats.sem(np.array(By)), np.std(np.array(By))
mz, s1z, s2z = np.mean(np.array(Bz)), scipy.stats.sem(np.array(Bz)), np.std(np.array(Bz))
hx = s1x * sp.stats.t._ppf((1+confidence)/2., len(Bx)-1)
hy = s1y * sp.stats.t._ppf((1+confidence)/2., len(By)-1)
hz = s1z * sp.stats.t._ppf((1+confidence)/2., len(Bz)-1)

#calculating
#print "Mean, CI of Bootstraped Xs",mx, 2*hx, mx-hx,mx+hx

##########################
##### ploting
#########################

xx=[]
for i in range(len(Bx)):
        xx.append(i)

lx = [mx for i in xrange(len(Bx))]
ly = [my for i in xrange(len(By))]
lz = [mz for i in xrange(len(Bz))]

lx1 = [mx-hx for i in xrange(len(Bx))]
lx2 = [mx+hx for i in xrange(len(Bx))]
ly1 = [my-hy for i in xrange(len(By))]
ly2 = [my+hy for i in xrange(len(By))]
lz1 = [mz-hz for i in xrange(len(Bz))]
lz2 = [mz+hz for i in xrange(len(Bz))]

'''
fig1=pl.figure()
pl.subplot(3, 1, 1)  # 2x1 grid, first plot
pl.scatter(xx, Bx,c='b',s=1.5,label='x')
pl.plot(xx,lx,c='b')
pl.plot(xx,lx1,alpha=0.8,c='k')
pl.plot(xx,lx2,alpha=0.8,c='k')
pl.legend(loc='upper right')
pl.title('ML estimators of the force.')

pl.subplot(3, 1, 2)  # 2x1 grid, second plot
pl.scatter(xx, By,c='g',s=1.5,alpha=0.5,label='y')
pl.plot(xx,ly,c='g')
pl.plot(xx,ly1,alpha=0.8,c='k')
pl.plot(xx,ly2,alpha=0.8,c='k')
pl.legend(loc='upper right')
pl.xlabel('x')

pl.subplot(3, 1, 3)  # 2x1 grid, second plot
pl.scatter(xx, Bz,c='r',s=1.5,alpha=0.5,label='z')
pl.plot(xx,lz,c='r')
pl.plot(xx,lz1,alpha=0.8,c='k')
pl.plot(xx,lz2,alpha=0.8,c='k')
pl.legend(loc='upper right')
pl.xlabel('Number of bootstrap samples')

pl.show()
pl.close(fig1)
'''

#################################
#formating and saving
#################################
fx= ("{0:.6f}".format(mx))
fy= ("{0:.6f}".format(my))
fz= ("{0:.6f}".format(mz))

ex= ("{0:.6f}".format(hx))
ey= ("{0:.6f}".format(hy))
ez= ("{0:.6f}".format(hz))

sx= ("{0:.6f}".format(s2x))
sy= ("{0:.6f}".format(s2y))
sz= ("{0:.6f}".format(s2z))


# mean(x), mean(y), mean(z), ci(x), ci(y), ci(z), std(x), std(y), std(z)  
print fx, fy, fz, ex, ey, ez, sx, sy, sz 

