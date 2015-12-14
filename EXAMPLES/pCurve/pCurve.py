#! /usr/bin/env python
'''
'''
###################################################
# Rolando P. Hong Enriquez (rolando.hong@iit.it) ##
###################################################
import argparse
import sys
import csv
import math
import os
#import matplotlib.pyplot as plt
import pylab as pl
import numpy as np
#from stringfunc import normvect, diffstring
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from functions import load_data, file_lines, file_cols, dm_vec, midpoint, getX2near, getX2far, maxvec, voroid_mean, string_ini, reparam, smooth_DIM, savepcurve 
import random
import pickle 

#######################################################
###### Parsing arguments from the command line ########
#######################################################
def _make_parser():
	parser = argparse.ArgumentParser(description="Calculates the principal curve (pcurve) of a group of 'M' collective variables (columns) measured 'N' times (e.g., extracted from a MD trajectory, rows). The method evolves an initial guessed  string of 'Min' points to find the pcurve. At each new iteration, the new points on the string are estimated as the mean value of the DIM-dimensional voronoi cells. The new points are subsequently reparametrized to a string  with equidistant points which is also smoothed before the new iteration. The final converged string (or pcurve) can be optionally reparaterized to a higher number of points 'Mout'. Several plots are produced including the convergence of the iterations and an estimated free energy along the pcurve. The final state of the calculations is saved to the file 'pcurve_state.dat', which serves as input to other tools (e.g.,, map2curve.py).",epilog="Thanks for using pCurve.py!!")
	parser.add_argument('INPUT', help="Input file")
	parser.add_argument('Min', type=int, help="Number of points on the path (or string) located on the free energy landscape and joining the initial anf final states of the system.")
	parser.add_argument('Mout', type=int, help="The points in Min can be reparametrize to Mout points after convergence")
	parser.add_argument('smooth', type=float, help="A factor between [0,1] defining the strength of the smoothing performed after each reparametrization. ")
	parser.add_argument('Niter', type=int, help="Number of iterations updating the strings.")
	parser.add_argument('-CV', action='append', dest='CV2plot', default=[], help="Indexes between [0, (DIM-1)] of the CVs to be visualized (2D and 3D options are possible).")
	parser.add_argument('-plotdata', action='append', dest='VAL', default=[], help="Default is 'yes'. If the value 'no' is specified then the intermediate plots will be produced with the curves only.")
	return parser
#--------------------------------
_p    = _make_parser()
__doc__ += _p.format_help()


##################################
if __name__ == '__main__':
        args = _p.parse_args()

	########################################
	# Error handling for the visualization
	#######################################

	if args.CV2plot==[]:
		sys.exit('Error (-CV): Please enter the CVs to be analysed.')
	if (len(args.CV2plot)==1) or (len(args.CV2plot)>3):
		sys.exit('Error (-CV): only 2D or 3D visualizations are possible.')	

	#Getting axis labels for visualization
	######################################
	lab=[]
	if (len(args.CV2plot)==2):
		exec "lab.append('CV %s')" % (args.CV2plot[0])
		exec "lab.append('CV %s')" % (args.CV2plot[1])
		CVX=int(args.CV2plot[0])
		CVY=int(args.CV2plot[1])
		labx=lab[0]
		laby=lab[1]
	else:
		exec "lab.append('CV %s')" % (args.CV2plot[0])
		exec "lab.append('CV %s')" % (args.CV2plot[1])
		exec "lab.append('CV %s')" % (args.CV2plot[2])
		CVX=int(args.CV2plot[0])
		CVY=int(args.CV2plot[1])
		CVZ=int(args.CV2plot[2])
		labx=lab[0]
		laby=lab[1]
		labz=lab[2]
	# error handling for the optional -plotdata argument 
	plotdata='yes'
	if (args.VAL==[]):
		plotdata='yes'
	else:
		dd=args.VAL[0]
		if dd=='no':
			#print 'no data willbe plotted...'
			plotdata=dd
		elif (dd=='yes') or (args.VAL==[]):
			#print 'plot everything...'
			plotdata='yes'
		else:
			sys.exit("Error (-plotdata): Only 'yes' or 'no' values are admisible.")
	

	##########################
	####### MAIN PROGRAM #####
	##########################
	sys.stdout.write("\a")
	#loading the input file ('N' rows 'M' columns or colective variables...CVs)
	N=file_lines(args.INPUT)
	M=file_cols(args.INPUT)
	P1=load_data(args.INPUT,N,M)

	# estimate the initial string (P2) of size Min
	P2 = string_ini(P1,M,args.Min)
	Pini=string_ini(P1,M,args.Min)
	PiniT=zip(*Pini)

	# maximun distances of the data to the images of the string
	maxdistances=maxvec(P1,P2)

	# Estimate the mean of the voronoid cells
	P3, fail, free_energy, allpoints, closeP, closeI = voroid_mean(P1,P2,M,maxdistances, 0)

	# Reparametrization and smoothing of the centers of the voronoid cells
	dist_to_X2, totalD = dm_vec(P3,args.Min)
	Pin=np.asarray(P3)
	Plist=reparam(Pin,M,args.Min,dist_to_X2)

	#smooth
	Ps = smooth_DIM(Plist,args.smooth,M)

	Prep=zip(*Plist)
	PinT=zip(*Pin)
	PsT=zip(*Ps)

	# Iterations for N-dimensional case (3D visualization)
	#######################################################
	if (len(args.CV2plot)==3):
		fig = plt.figure()
		ax = fig.add_subplot(111,projection='3d')
		ax.set_title('Curve evolution', alpha=0.6)
		converg=[]
		for i in range(args.Niter):
			print 'Iteration ', (i+1), 'out of ', args.Niter
			# get voronoid centers
			maxdistances=maxvec(P1,Ps)
			Pm, fail, free_energy, allpoints, closeP, closeI = voroid_mean(P1,Ps,M,maxdistances,0)

			# Reparametrization and smoothing of the centers of the voronoid cells
			dist_to_X2, totalD = dm_vec(Pm,args.Min)
			Pin=np.asarray(Pm)
			Plist=reparam(Pin,M,args.Min,dist_to_X2)

			#smooth
			Ps1 = smooth_DIM(Plist,args.smooth,M)
			d=np.linalg.norm(np.asarray(Ps)-np.asarray(Ps1))
			converg.append(d)
			Ps=Ps1
			PsT=zip(*Ps)
			ax.plot(PsT[CVX],PsT[CVY],PsT[CVZ],c='g',linewidth=2,linestyle=':')

		if plotdata=='no':
			ax.plot(PiniT[CVX],PiniT[CVY],PiniT[CVZ],c='b',linewidth=2,label='curve 1')
			ax.plot(PsT[CVX],PsT[CVY],PsT[CVZ],c='r',linewidth=2,label='curve N')
			ax.scatter(PsT[CVX],PsT[CVY],PsT[CVZ],c='r',marker='o',s=20)
			ax.set_xlabel(labx)
			ax.set_ylabel(laby)
			ax.set_zlabel(labz)
			ax.legend(loc=3,fontsize=12,fancybox=True).get_frame().set_alpha(0.5)    
			plt.show()
		else:
	       	 	PP1=zip(*P1)
	        	ax.scatter(PP1[CVX],PP1[CVY],PP1[CVZ],c='m',marker='o',s=0.5)
	        	ax.plot(PiniT[CVX],PiniT[CVY],PiniT[CVZ],c='b',linewidth=2,label='curve 1')
	        	ax.plot(PsT[CVX],PsT[CVY],PsT[CVZ],c='r',linewidth=2,label='curve N')
	        	ax.scatter(PsT[CVX],PsT[CVY],PsT[CVZ],c='r',marker='o',s=20)
	        	ax.set_xlabel(labx)
	        	ax.set_ylabel(laby)
	        	ax.set_zlabel(labz)
			ax.legend(loc=3,fontsize=12,fancybox=True).get_frame().set_alpha(0.5)
	        	plt.show()

	else:
        	converg=[]
	        for i in range(args.Niter):
	                print 'Iteration ', (i+1), 'out of ', args.Niter
	                # get voronoid centers
	                maxdistances=maxvec(P1,Ps)
	                Pm, fail, free_energy, allpoints, closeP, closeI = voroid_mean(P1,Ps,M,maxdistances,0)

	                # Reparametrization and smoothing of the centers of the voronoid cells
	                dist_to_X2, totalD = dm_vec(Pm,args.Min)
	                Pin=np.asarray(Pm)
	                Plist=reparam(Pin,M,args.Min,dist_to_X2)

	                #smooth
	                Ps1 = smooth_DIM(Plist,args.smooth,M)
	                d=np.linalg.norm(np.asarray(Ps)-np.asarray(Ps1))
	                converg.append(d)
	                Ps=Ps1
	                PsT=zip(*Ps)
			pl.plot(PsT[CVX],PsT[CVY],c='g',linewidth=2,linestyle=':')

		if plotdata=='no':
	        	pl.plot(PiniT[CVX],PiniT[CVY],c='b',linewidth=2,label='curve 1')
	        	pl.plot(PsT[CVX],PsT[CVY],c='r',linewidth=2,label='curve N')
        		pl.scatter(PsT[CVX],PsT[CVY],c='r',marker='o',s=20)
        		pl.xlabel(labx)
        		pl.ylabel(laby)
			pl.legend(loc=3,fontsize=10,fancybox=True)
        		pl.show()
		else:
			PP1=zip(*P1)
	                pl.scatter(PP1[CVX],PP1[CVY],c='m',marker='o',s=0.5)
	                pl.plot(PiniT[CVX],PiniT[CVY],c='b',linewidth=2,label='curve 1')
	                pl.plot(PsT[CVX],PsT[CVY],c='r',linewidth=2,label='curve N')
	                pl.scatter(PsT[CVX],PsT[CVY],c='r',marker='o',s=20)
	                pl.xlabel(labx)
	                pl.ylabel(laby)
			pl.legend(loc=3,fontsize=10,fancybox=True)
	                pl.show()
		
	#Free energyplot
	XF=[]
	for i in range(1,len(Ps)+1):
	        XF.append(i)
	pl.plot(XF,free_energy,c='r', linewidth=2)
	pl.ylabel('Free Energy(kcal/mol)')
	pl.xlabel('Points on the Curve')
	pl.title('Free Energy on the Principal Curve')
	pl.xlim(1,len(Ps))
	pl.ylim(0,max(free_energy)+0.1)
	pl.fill_between(XF, free_energy, facecolor='blue', alpha=0.1)
	pl.show()

	#Convergence plot
	XX=[]
	for i in range(1,args.Niter+1):
	        XX.append(i)
	pl.ylabel('Curve[i]-Curve[i-1]')
	pl.xlabel('Iterations')
	pl.xlim(1,args.Niter)
	pl.ylim(0,max(converg)+0.1)
	pl.title('Convergence')
	pl.scatter(XX,converg,c='k', marker='o')
	pl.plot(XX,converg,c='k', linewidth=3)
	pl.fill_between(XX, converg, facecolor='blue', alpha=0.1)
	pl.show()

	#Plot to compare the last pcurve (Min points) and the corresponding closest 'Min' points.	
	PC=zip(*closeP)
	fig = plt.figure()
	ax = fig.add_subplot(111,projection='3d')
	ax.set_title('Actual data of the closest string to the final Min pcurve', alpha=0.6)
	ax.scatter(PC[CVX],PC[CVY],PC[CVZ],c='b',marker='o',s=20)
	ax.plot(PC[CVX],PC[CVY],PC[CVZ],c='b',linewidth=2,label='Closest String')
	ax.plot(PsT[CVX],PsT[CVY],PsT[CVZ],c='r',linewidth=2,label='Principal Curve')
	ax.scatter(PsT[CVX],PsT[CVY],PsT[CVZ],c='r',marker='o',s=20)
	ax.set_xlabel(labx)
	ax.set_ylabel(laby)
	ax.set_zlabel(labz)
	ax.legend(loc=3,fontsize=12,fancybox=True).get_frame().set_alpha(0.5)
	plt.show()

	# Reparametrization to Mout points after the iterations
	dist_to_X2, totalD = dm_vec(Ps,args.Mout)
	Pin=np.asarray(Ps)
	Plist=reparam(Pin,M,args.Mout,dist_to_X2)
	Psout = smooth_DIM(Plist,args.smooth,M)
	PsoutT=zip(*Psout)

	###########################################################
	#### plot of last reparametrization to Mout (2D and 3D)
	#legend
	lin='Min = '
	lin += `args.Min`
	lout='Mout = '
	lout += `args.Mout`

	if (len(args.CV2plot)==3):
		#### plot of last reparametrization to Mout (3D)
		fig = plt.figure()
		ax = fig.add_subplot(111,projection='3d')
		ax.plot(PsT[CVX],PsT[CVY],PsT[CVZ],c='b',marker='o',label=lin)
		ax.plot(PsoutT[CVX],PsoutT[CVY],PsoutT[CVZ],c='r',marker='o',label=lout)
		ax.set_xlabel(labx)
		ax.set_ylabel(laby)
		ax.set_zlabel(labz)
		ax.set_title('Reparametrization to Mout points after the iterations', alpha=0.6)
		ax.legend(loc=3,fontsize=10,fancybox=True).get_frame().set_alpha(0.5)
		plt.show()
	else:
		#### plot of last reparametrization to Mout (2D)
		pl.plot(PsT[CVX],PsT[CVY],c='b',marker='o',label=lin)
		pl.title('Reparametrization to Mout points after the iterations',alpha=0.5)
		pl.xlabel(labx)
		pl.ylabel(laby)
		pl.plot(PsoutT[CVX],PsoutT[CVY],c='r',marker='o',label=lout)
		pl.legend(loc=3,fontsize=10,fancybox=True)
		pl.show()

	#######################################################################################
	# Plot to compare the las reparametrization to Mout and the closest corresponding Mout
	#points in the data
	#(1) find the points the their indexes
	maxdistances=maxvec(P1,Psout)
	closePP=[]
	closeII=[]
	for j in range(len(Psout)):
		for i in range(len(P1)):
			dist=np.linalg.norm(P1[i]-Psout[j])
			if i==0:
				minid=maxdistances[j]
				index=0
			if dist < minid:
				minid=dist
				index=i
		closePP.append(P1[index])
		closeII.append(index)

	#(2) Plot
	PCC=zip(*closePP)
	fig = plt.figure()
	ax = fig.add_subplot(111,projection='3d')
	ax.set_title('Actual data of the closest string to the final Mout pcurve', alpha=0.6)
	ax.scatter(PCC[CVX],PCC[CVY],PCC[CVZ],c='b',marker='o',s=20)
	ax.plot(PCC[CVX],PCC[CVY],PCC[CVZ],c='b',linewidth=2,label='Closest String')
	ax.plot(PsoutT[CVX],PsoutT[CVY],PsoutT[CVZ],c='r',linewidth=2,label='Principal Curve')
	ax.scatter(PsoutT[CVX],PsoutT[CVY],PsoutT[CVZ],c='r',marker='o',s=20)
	ax.set_xlabel(labx)
	ax.set_ylabel(laby)
	ax.set_zlabel(labz)
	ax.legend(loc=3,fontsize=12,fancybox=True).get_frame().set_alpha(0.5)
	plt.show()
	
	################################################
	##### PLOTING of the final curve with data #####
	if fail==1:
		print 'Failed run... (no output)'
	else:
		PP1=zip(*P1)
		PP2=zip(*P2)
		PP3=zip(*P3)
		#ploting
	if (len(args.CV2plot)==3):
		#if args.DIM==3:
		fig = pl.figure()
		ax = fig.add_subplot(111,projection='3d')
	       	ax.set_xlabel(labx)
	       	ax.set_ylabel(laby)
	       	ax.set_zlabel(labz)
	       	ax.set_title('Final curve of Mout points with data', alpha=0.6)
		ax.scatter(PP1[CVX],PP1[CVY],PP1[CVZ],c='m',marker='o',s=0.5,alpha=0.5)
		ax.scatter(PsoutT[CVX],PsoutT[CVY],PsoutT[CVZ],c='r',marker='o',s=25)
		ax.plot(PsoutT[CVX],PsoutT[CVY],PsoutT[CVZ],c='r',linewidth=3)
		pl.show()
	else:
		#if args.DIM==2:
		pl.title('Final curve of Mout points with data',alpha=0.6)
		pl.scatter(PP1[CVX],PP1[CVY],c='m',marker='o',s=1.5,alpha=0.5)
		pl.scatter(PsoutT[CVX],PsoutT[CVY],c='r',marker='o',s=25)
		pl.plot(PsoutT[CVX],PsoutT[CVY],c='r',linewidth=3)
	       	pl.xlabel(labx)
	       	pl.ylabel(laby)
		pl.show()

	#########################################################################
	##### Save selected variables as one "state" of the calculations ########
	#########################################################################
	savepcurve(N,M,args.Min,args.Mout,args.smooth,args.Niter,Ps,Psout,converg,free_energy,allpoints,closeP,closeI,closePP,closeII)

