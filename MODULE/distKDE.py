#! /usr/bin/env python

import csv
import argparse
import numpy as np
import pylab as pl
from functions import file_cols, commonArea
from scipy.stats import kde

############################################################################################
parser = argparse.ArgumentParser(description="This perform the following tasks: (1) load the data from columns r1 and r2 in two vectors, (2) makes an estimation of the probability density for the loaded vectors using Kernel Density Estimation (KDE) with a gaussian kernel, (3) plots the estimated probability densities for re and r1 as well as the common area between these distributions, and (4) estimates the total variation distance to quantify the difference between the KDE calculated probability densities",epilog="Thanks for using distKDE!!")
parser.add_argument('IN_file', help="input file.")
parser.add_argument('r1', type=int, help="Index of the first column.")
parser.add_argument('r2', type=int, help="Index of the second column.")
args = parser.parse_args()
############################################################################################

#Read only 2 columns the data file
nc = file_cols(args.IN_file)
x1=[]
x2=[]
with open(args.IN_file) as file:
	reader = csv.reader(file, delimiter=' ',skipinitialspace=True)
	for row in reader:
		for i in range(nc):
			if (i==args.r1):
				exec "x1.append(float(row[%s]))" % (i)
			elif (i==args.r2):
				exec "x2.append(float(row[%s]))" % (i)

xt=x1+x2
x1=np.asarray(x1)
x2=np.asarray(x2)
xt=np.asarray(xt)				

#kernel density stimation (non-parametric prob distrib)
NN=100
density1 = kde.gaussian_kde(x1)
xgrid1 = np.linspace(x1.min(), x1.max(), NN)
density2 = kde.gaussian_kde(x2)
xgrid2 = np.linspace(x2.min(), x2.max(), NN)
densityt = kde.gaussian_kde(xt)
xgridt = np.linspace(xt.min(), xt.max(), NN)

#Common area
ac, AC = commonArea(density1,density2,xgridt)
print 'AC',AC

#ploting

pl.plot(xgridt, density1(xgridt), 'b-',label='Column r1')
pl.plot(xgridt, density2(xgridt), 'r-', label='Column r2')
pl.plot(x1, np.zeros(x1.shape), 'b+', ms=12) #rug plot
pl.plot(x2, np.zeros(x2.shape), 'r+', ms=12) #rug plot
pl.plot(xgridt, ac, 'm-')
pl.fill_between(xgridt, ac, facecolor='magenta', alpha=0.3)
pl.xlabel("X")
pl.ylabel("Prob. Density Function")
pl.title("Kernel Density Estimation (KDE)")
pl.legend(loc="upper right")
pl.show()

