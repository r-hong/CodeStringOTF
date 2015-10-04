#! /usr/bin/env python
'''
'''
###################################################
# Rolando P. Hong Enriquez (rolando.hong@iit.it) ##
###################################################
import argparse
import pickle 
import math
from functions import load_data, load_data_range,loadpcurve 
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from collections import defaultdict, Counter

#######################################################
###### Parsing arguments from the command line ########
#######################################################
def _make_parser():
	parser = argparse.ArgumentParser(description="This tool maps collective variables (CVs) from an input file on a principal curve. Each CV in the data file should be place in a column. The rows in the CV data file should coincide in time with the CVs used to estimate the pricipal curve (as a particular case, the same CV data file can be used with both pcurve.py and map2curve.py). The initial and final indexes of the CVs to be mapped should be provided. The information of the principal curve is automatically loaded from the file 'pcurve_state.dat' (output from 'pcurve.py' ). To extract CVs from MD trajectories using VMD see the file getCVs.tcl.",epilog="Thanks for using map2curve.py!!")
	parser.add_argument('pcurveFile', help="Input file with the data from a calculated principal curve.")
	parser.add_argument('cvFile', help="Input file with CVs.")
	parser.add_argument('iniCV', type=int, help="Index (column) of the first CV to be mapped.")
	parser.add_argument('endCV', type=int, help="Index (column) of the last CV to be mapped.")
        return parser
#------------------------------
_p    = _make_parser()
__doc__ += _p.format_help()
######################################################
if __name__ == '__main__':
        args = _p.parse_args()

	#######################################
	############ MAP2CURVE ################
	####################################### 
	## loading pcurve state
	N, DIM, Min, Mout, smooth, Niter, Ps, Psout, converg, free_energy, allpoints, closeP, closeI, closePP, closeII = loadpcurve(args.pcurveFile) 
	# loading CVs to be mapped on the pcurve 
	P1=load_data_range(args.cvFile,N,args.iniCV,args.endCV)
	DIM_new=args.endCV-args.iniCV

	############################################################
	#Creating the matrix (Z values) of 'Min' X 'DIM_new' 
	#Calculation of mean and std values
	############################################################
	CVMT=[]
	CVDT=[]	
	CVDATAT=[]
	for q in range(Min):
		indexes=allpoints[q]
		CVM=[]
		CVD=[]
		CVDATA=[]
		for i in range(DIM_new):
			CVX= P1[:,i]
			vals_in_cell=[]	
			for j in range(N):
				if (j in indexes):
					vals_in_cell.append(CVX[j])
			media=np.mean(vals_in_cell)
			standev=np.std(vals_in_cell)
			CVM.append(media)
			CVD.append(standev)
			CVDATA.append(CVX[q])
		CVMT.append(CVM)
		CVDT.append(CVD)
		CVDATAT.append(CVDATA)

	ZZA=zip(*CVMT)
	ZZB=zip(*CVDT)
	ZZC=zip(*CVDATAT)
	zA=np.asarray(ZZA)
	zB=np.asarray(ZZB)
	zC=np.asarray(ZZC)

	###############################################
	#Creating 2 2D grids for the x and y bounds
	##############################################
	dx, dy = 1, 1
	y, x = np.mgrid[slice(1, DIM_new + dy, dy),
			slice(1, Min + dx, dx)]

	###################################
	### Contour  mean #################
	###################################
	zA = zA[:-1, :-1]
	levels = MaxNLocator(nbins=15).tick_values(zA.min(), zA.max())
	cmap = plt.get_cmap('PiYG')
	norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

	plt.contourf(x[:-1, :-1] + dx / 2.,
        	     y[:-1, :-1] + dy / 2., zA, levels=levels,
	             cmap=cmap)
	plt.colorbar()
	plt.title('Mean of the CVs along pcurve')
	plt.show()

	###################################
	### Contour  std #################
	###################################
	zB = zB[:-1, :-1]
	levels = MaxNLocator(nbins=15).tick_values(zB.min(), zB.max())
	cmap = plt.get_cmap('PiYG')
	norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

	plt.contourf(x[:-1, :-1] + dx / 2.,
        	     y[:-1, :-1] + dy / 2., zB, levels=levels,
	             cmap=cmap)
	plt.colorbar()
	plt.title('Std of the CVs along pcurve')
	plt.show()

	#########################################################################
	### Contour  closest data point to the principal curve #################
	########################################################################
	zC = zC[:-1, :-1]
	levels = MaxNLocator(nbins=15).tick_values(zC.min(), zC.max())
	cmap = plt.get_cmap('PiYG')
	norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

	plt.contourf(x[:-1, :-1] + dx / 2.,
	             y[:-1, :-1] + dy / 2., zC, levels=levels,
	             cmap=cmap)
	plt.colorbar()
	plt.title('Closest point to the Min pcurve')
	plt.show()

	####################################################################
	#A different calculation must be done for Mout
	############################################################
	CV_M01=[]
	for q in range(Mout):
	        CV_M02=[]
	        for i in range(DIM_new):
			CV_M02.append(P1[(closeII[q]),i])
	        CV_M01.append(CV_M02)

	ZZD=zip(*CV_M01)
	zD=np.asarray(ZZD)

	###############################################
	#Creating 2 2D grids for the x and y bounds
	##############################################
	dx, dy = 1, 1
	y, x = np.mgrid[slice(1, DIM_new + dy, dy),
	                slice(1, Mout + dx, dx)]

	################################################################
	### Contour  closest points to the Mout pcurve #################
	################################################################
	zD = zD[:-1, :-1]
	levels = MaxNLocator(nbins=15).tick_values(zD.min(), zD.max())
	cmap = plt.get_cmap('PiYG')
	norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

	plt.contourf(x[:-1, :-1] + dx / 2.,
	             y[:-1, :-1] + dy / 2., zD, levels=levels,
	             cmap=cmap)
	plt.colorbar()
	plt.title('Closest points to the Mout pcurve')
	plt.show()

	### Note: this list of indexes (closeII) is not built from the voronoi cells,
	### therefore the same point could be the closest one to different points on the Mout
	### principal curve. We here print a warning in case of duplicate points in the previous 
	### visualization
	ff=[]
	ff=[i for i, x in enumerate(closeII) if closeII.count(x) > 1]

	if ff != []:
		print "Warning. The existence of a unique neighbor point can be only garanteed for the points on the 'Min' principal curve calculated with pcurve.py (as this curve was calculated from voronoi cells). This warning raises because there are points on the 'Mout' curve sharing the same closest data point."
		print "The affected indexes on the Mout principal curve are:"
		print ff







