#! /usr/bin/env python
'''
'''
###################################################
# Rolando P. Hong Enriquez (rolando.hong@iit.it) ##
###################################################
import argparse
import pickle 
import sys
import math
from functions import load_data, load_data_range,loadpcurve 
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from collections import defaultdict, Counter

#######################################################
###### Parsing arguments from the command line ########
#######################################################
def _make_parser():
	parser = argparse.ArgumentParser(description="This tool displays (and/or saves) the content of a pcurve state file (output from pcurve.py).",epilog="Thanks for using check_pcurve.py!!")
	parser.add_argument('file1', help="Input pcurve_state file.")
	parser.add_argument('-plot', action='append', dest='PL', default=[], help="Plot data? [yes/no].")
	parser.add_argument('-save', action='append', dest='SV', default=[], help="Save curves? [yes/no].")
	return parser
#--------------------------------
_p    = _make_parser()
__doc__ += _p.format_help()


######################################################
if __name__ == '__main__':
        args = _p.parse_args()

	##############################
	# error handling
	#############################
	if args.PL==[]:
        	sys.exit('Error (-plot): Please specify if you want to plot the content of the state file.')
	if (len(args.PL)>1):
        	sys.exit('Error (-plot): only a single value [yes/no] is allowed.')
	if args.SV==[]:
        	sys.exit('Error (-save): Please specify if you want to save the curves of the state file.')
	if (len(args.SV)>1):
        	sys.exit('Error (-save): only a single value [yes/no] is allowed.')

	pl_val=args.PL[0]
	sv_val=args.SV[0]

	if (pl_val != 'yes') and (pl_val != 'no'):
		sys.exit('Error (-plot): Wrong value for -plot. Correct values: [yes/no].')
	if (sv_val != 'yes') and (sv_val != 'no'):
        	sys.exit('Error (-save): Wrong value for -save. Correct values: [yes/no].')

	#######################################
	############ CHECK_PCURVE ################
	####################################### 
	## loading pcurve state 1 and 2
	N, DIM, Min, Mout, smooth, Niter, Ps, Psout, converg, free_energy, allpoints, closeP, closeI, closePP, closeII = loadpcurve(args.file1)
	print '#######################################################' 
	print '############# pcurve run parameters  ##################'
	print '#######################################################'
	print '                                               '
	print 'N      (Size of the CV vectors)            = ',N
	print 'DIM    (# of CVs)                          = ',DIM
	print 'Min    (#of points on the string)          = ',Min
	print 'Mout   (# of final reparametrized points)  = ',Mout
	print 'smooth (strength of the smoothing [0-1])   = ',smooth
	print 'Niter  (# of iterations)                   = ',Niter
	print '                                             '
	print '#######################################################'

	# loading original data points
	#loading data points ('N' points for each one of the 'DIM' CVs)
	#P1=load_data(args.file3,N1,DIM1)

	if pl_val=='yes':
		print 'we are ploting...'
	if sv_val=='yes':
		print 'we are saving...'


	if pl_val=='yes':
		#plot smoothed Min and Mout curves
		PsoutT=zip(*Psout)
		PsT=zip(*Ps)
		fig = plt.figure()
		ax = fig.add_subplot(111,projection='3d')
		ax.plot(PsT[0],PsT[1],PsT[2],c='b',marker='o',label='Min curve')
		ax.plot(PsoutT[0],PsoutT[1],PsoutT[2],c='r',marker='o',label='Mout curve')
		ax.set_xlabel('X')
		ax.set_ylabel('Y')
		ax.set_zlabel('Z')
		ax.set_title('Min vs.Mout principal curves', alpha=0.6)
		ax.legend(loc=3,fontsize=10,fancybox=True).get_frame().set_alpha(0.5)
		plt.show()

		#Free energy plot
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
		for i in range(1,Niter+1):
        		XX.append(i)
		pl.ylabel('Curve[i]-Curve[i-1]')
		pl.xlabel('Iterations')
		pl.xlim(1,Niter)
		pl.ylim(0,max(converg)+0.1)
		pl.title('Convergence')
		pl.scatter(XX,converg,c='k', marker='o')
		pl.plot(XX,converg,c='k', linewidth=3)
		pl.fill_between(XX, converg, facecolor='blue', alpha=0.1)
		pl.show()


	if sv_val=='yes':
		#saving coordinates for pcurve (Mout) 1
		file1='x.dat'
		file2='y.dat'
		file3='z.dat'
		x=(PsoutT[0])
		y=(PsoutT[1])
		z=(PsoutT[2])
		with open(file1, 'w') as file:
        		for item in x:
                		file.write('{}\n'.format(item))

		with open(file2, 'w') as file:
        		for item in y:
                		file.write('{}\n'.format(item))
		with open(file3, 'w') as file:
        		for item in z:
                		file.write('{}\n'.format(item))

		#save the indexes of the points in the data thatare closest to the Mout pcurve
		file1='closeII.dat'
		with open(file1, 'w') as file:
        		for item in closeII:
                		file.write('{}\n'.format(item))


