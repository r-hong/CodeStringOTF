#! /usr/bin/env python
"""
"""
import argparse
import csv
import math
import pylab as pl
import numpy as np
from functions import *
 
def _make_parser():
	p = argparse.ArgumentParser(description="Tool to estimate the convergence of the string method for the calculation of the minimum free energy path. The calculation is based of the succesive differences between strings in the space of the normalized CVs. In the input files: (1) each column represents a CV, (2) each row represents a point along the path (string), and (3) the columns should be separated by a space character.", epilog="Thanks for using convergence.py!!")
	p.add_argument('Niter', type=int, help="number of iterations files updating the strings",default=100)
	p.add_argument('IN_file', help="Prefix of the files containing the string iterations. Example: if you have the files  file1 file2, etc., then IN_file='file'")
	return p
#--------------------------------
_p    = _make_parser()
__doc__ += _p.format_help()
######################################################
if __name__ == '__main__':
	args = _p.parse_args()
 
	#get Number of rows and number of columns from the first file
	fileA='%s''%s''.dat' % (args.IN_file,0)
	Ncv=file_cols(fileA)
	Npoints=file_lines(fileA)

	####################################################################
	# Initializing a temp array ###
	tf=[]
	for i in range(0,args.Niter):
		tf.append('')
	for i, value in enumerate(tf):
		exec "numS%s=[]" % (i)
	converg1=[0]*args.Niter
	#converg2=[0]*args.Niter
	X=[0]*args.Niter

	#### Calculating the convergence for all string files
	converg1=[]
	for j in range(0,args.Niter):
		fileB='%s''%s''.dat' % (args.IN_file,j+1)
		dd1=diffstring(fileA,fileB,Ncv,Npoints)
		converg1.append(dd1)

	converg2=normvect(np.asarray(converg1))
	final=[0]*len(converg1)
	for j in range(0,args.Niter):
		final[j]=1-converg2[j]
		X[j]=j

	pl.title('Normalized Frechet distance',alpha=0.6)
	pl.xlabel('Iterations')
	pl.ylabel('Convergence')
	pl.plot(X, final)
	pl.show()
	#####################################################################
	#plot per image
	####################################################################
	ImgX=[0]*Npoints
	for j in range(Npoints):
        	ImgX[j]=j+1

	import colorsys
	def get_color(color):
		for hue in range(color):
			hue = 1. * hue / color
			col = [int(x) for x in colorsys.hsv_to_rgb(hue, 1.0, 230)]
			yield "#{0:02x}{1:02x}{2:02x}".format(*col)

	color = get_color(args.Niter)
	dist_all=[]
	for i in range(args.Niter):
		fileB='%s''%s''.dat' % (args.IN_file,i+1)
		dist=[]
		for j in range(Npoints):
			dist.append(diffpoint(fileA,fileB,Ncv,Npoints,j+1))
		dist_all.append(dist)
	
		acolor = next(color)
		if (i==0):
			pl.plot(ImgX,dist_all[i],color='blue',linewidth=3)
		elif (i==(args.Niter-2)):
			pl.plot(ImgX,dist_all[i],color='red',linewidth=3)
		else:
			pl.plot(ImgX,dist_all[i],color=acolor,ls=':')
		pl.title('Evolution per image')
		pl.xlabel('Images')
		pl.ylabel('rmsd')
	pl.show()

