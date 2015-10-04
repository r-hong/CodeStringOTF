#! /usr/bin/env python
"""
This script implement a brute force approach to calculate the Sammon mapping [1] for a high dimansional vector that is given in an input file IN_file. A range of dimensions (columns) from CVini to CVend can be selected and a limited number of points (rows) equal to N can be selected too.
Sammon Mapping is a dimensionality reduction technique. Unlike principal component analysis (PCA) that performs a simple linear projection of the data while trying to keep as much variance from the original data as possible without any concern about the internal geometrical structure of the data, Sammon mapping tackles this problem. Simply put, a Sammon mapping minimizes the differences between the inter-point distances in the two spaces: the original spaced and the space resulting from the mapping. Furthermore, this type of mapping tries to preserved the topology of the original space.
-----------------------
[1] Sammon JW Jr. "A nonlinear mapping for data structure analysis". IEEE Transactions on Computers, vol C-18, No. 5: 401-409 (1969)  
"""

import argparse
import pickle
import math
import random
from functions import load_data_range, delrn, getE
import pylab as pl
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from operator import itemgetter


def _make_parser():
        p = argparse.ArgumentParser()
        p.add_argument('N', type=int, help="Number of points",default=50)
	p.add_argument('CVini', type=int, help="Initial column",default=1)
	p.add_argument('CVend', type=int, help="Number of points",default=100)
        p.add_argument('IN_file', help="Input file containing the collection of N dimensional points to be sammon projected in a low dimensioanl space")
        return p
#--------------------------------
_p    = _make_parser()
__doc__ += _p.format_help()


######################################################
if __name__ == '__main__':
        args = _p.parse_args()

	L=args.CVend-args.CVini+1 
	P1=load_data_range(args.IN_file,args.N,args.CVini,args.CVend)

	print L
	OI=[]
	dd=40
	for i in range(dd,L):
		print i

	EF=[]
	deletedCVs=[]
	for i in range(L-1):
		if i==0:
			Plast=P1

		#get maximal
        	for j in range (i,L):
                	Ptmp=delrn(Plast,j-i)
                	Etmp=getE(P1,Ptmp)
			if j==i:
				maxE=Etmp
				max_index=i
			else:
				if Etmp>maxE:
					maxE=Etmp
					max_index=j
	
		EE=[]
		minE=maxE
		
		for j in range (i,L):
			PP=delrn(Plast,j-i)
			E=getE(P1,PP)
			if E<minE:
				minE=E
				mini_index=j

		print 'mini_index', mini_index
		PPP=delrn(Plast,mini_index-i)
		Plast=PPP
		EF.append(minE)
		deletedCVs.append(mini_index-i)	


	print min(deletedCVs)
	print max(deletedCVs)
	print deletedCVs

	# Plot
	XF=[]
	for i in range(1,len(EF)+1):
        	XF.append(i)

	pl.plot(XF,EF,c='r', linewidth=2)
	pl.ylabel('Error')
	pl.xlabel('# of deleted CVs')
	pl.title('Brute force')
	pl.xlim(1,len(EF))
	pl.ylim(0,max(EF)+0.1)
	pl.fill_between(XF, EF, facecolor='blue', alpha=0.1)
	pl.show()

	####################################################
	#Plot of the eliminated indexes
	XF=[]
	for i in range(1,len(deletedCVs)+1):
        	XF.append(i)

	pl.plot(XF,deletedCVs,c='r', linewidth=2)
	pl.ylabel('Index of the delected CVs')
	pl.xlabel('# of deleted CVs')
	pl.title('Brute force')
	pl.xlim(1,len(deletedCVs))
	pl.ylim(0,max(deletedCVs)+0.1)
	pl.fill_between(XF, deletedCVs, facecolor='blue', alpha=0.1)
	pl.show()







