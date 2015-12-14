#! /usr/bin/env python
###################################################
# Rolando P. Hong Enriquez (rolando.hong@iit.it) ##
###################################################
import argparse
import pickle 
import math
from functions import * # load_data, load_data_range,loadpcurve 
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from collections import defaultdict, Counter
import sys
from scipy.stats.stats import pearsonr 
import pylab as pl
#######################################################
###### Parsing arguments from the command line ########
#######################################################
# example
#./map2img.py rmsdB_image20.dat rmsd_CA_image20.dat dist_dih_resid_image20.dat
parser = argparse.ArgumentParser(description="This tool bla bla bla.",epilog="Thanks for using bla bla...!!")
parser.add_argument('file1', help="Input file 1.")
parser.add_argument('file2', help="Input file 2.")
parser.add_argument('file3', help="Input file 3.")
#parser.add_argument('NN', type=int, help="NN.")
#parser.add_argument('CVend', type=int, help="Index (column) of the last CV to be mapped.")
args = parser.parse_args()


M1=file_cols(args.file1)
N1=file_lines(args.file1)

M2=file_cols(args.file2)
N2=file_lines(args.file2)

M3=file_cols(args.file3)
N3=file_lines(args.file3)



# loading CVs to be mapped on the pcurve 
P1=load_data(args.file1,N1,M1)
P2=load_data(args.file2,N2,M2)
P3=load_data(args.file3,N3,M3)

a=zip(*P1)
b=zip(*P2)
c=zip(*P3)
print len(a[0])
print len(b[1])
print len(c[0])

print pearsonr(a[0],b[0])
#sys.exit("stop here")

x=[]
corr_coef=[]
mean_dist=[]
for i in range(M2):
	corr_coef.append(pearsonr(a[0],b[i]))
	x.append(i)
	mean_dist.append(sum(c[i])/N1)

corr_coef1=zip(*corr_coef)
ss=corr_coef1[0]
mdist=[]
for j in range(len(mean_dist)):
	if mean_dist[j]<10:
		mdist.append(ss[j])
	else:
		mdist.append(0)



pl.title('Ligand RMSD vs. residue_id RMSD',alpha=0.6)
pl.xlabel('Residue #')
pl.ylabel('Pearson Corr. Coef')
pl.xlim(min(x), max(x))
pl.bar(x, corr_coef1[0],color='black',alpha=0.6)
for j in range(len(mean_dist)):
        if (mean_dist[j]>5) and (mean_dist[j]<10):
		pl.scatter(x[j],ss[j],color='red')
#pl.savefig('corr20.png')
pl.show()


#pl.title('mean distances',alpha=0.6)
#pl.xlabel('Residue #')
#pl.ylabel('mean dist')
#pl.plot(x, mean_dist)
#pl.show()




sys.exit("stop here")

