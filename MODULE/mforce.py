#! /usr/bin/env python
"""
Given afile with 3 columns corresponding to the force components on the x, y, and z coordinates this scrip calculates (and plots) the expected values of the mean forces. 
"""
###################################################
# Rolando P. Hong Enriquez (rolando.hong@iit.it) ##
###################################################
import argparse
import os
import csv
import pylab as pl
import numpy as np
import math
from functions import file_cols, readFile

def _make_parser():
	parser = argparse.ArgumentParser(description="Calculate the mean forces in time.",epilog="Thanks for using mforces.py!!")
	parser.add_argument('IN_file', help="Input file")
	parser.add_argument('OUT_img1', help="Output image (forces)")
	parser.add_argument('OUT_img2', help="Output image (cosine similarity)")
	parser.add_argument('N', type=int, help="Stride. Use only every Nth value from the input file.")
	return parser
#--------------------------------
_p    = _make_parser()
__doc__ += _p.format_help()
######################################################
if __name__ == '__main__':
        args = _p.parse_args()

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
	#### reading the file with forces (in colmns) ########
        nCol=file_cols('tmp.dat')
        vecF = readFile('tmp.dat',nCol)

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

	#calculating the cummulative mean values
	avex=[]
	avey=[]
	avez=[]
	avet=[]
	cosCUM=[]
	sumCS=0
	for i in range(1,len(c1)):
		avex.append(sum(c1[0:i])/len(c1[0:i]))
		avey.append(sum(c2[0:i])/len(c2[0:i]))
		avez.append(sum(c3[0:i])/len(c3[0:i]))
		avet.append(abs(avex[i-1])+abs(avey[i-1])+abs(avez[i-1]))
		cosSIM=cosine_similarity(vecF[0],vecF[i])
		#print cosSIM
		sumCS=sumCS+cosine_similarity(vecF[0],vecF[i])
		meanCS=sumCS/len(c1[0:i])
		cosCUM.append(meanCS)


	# ploting
	xx=[]
	for i in range(len(avex)):
		xx.append(i)

	fig1=pl.figure()
	pl.title('Mean force',alpha=0.6)
	pl.scatter(xx,avex,c='g',marker='o',s=1.5,alpha=0.5)
	pl.scatter(xx,avey,c='b',marker='o',s=1.5,alpha=0.5)
	pl.scatter(xx,avez,c='r',marker='o',s=1.5,alpha=0.5)
	#pl.scatter(xx,avet,c='k',marker='o',s=1.5,alpha=0.5)
	pl.plot(xx,avex,c='g',label='Fx')
	pl.plot(xx,avey,c='b',label='Fy')
	pl.plot(xx,avez,c='r',label='Fz')
	#pl.plot(xx,avet,c='k',label='abs(Fx)+abs(Fy)+abs(Fz)')
	pl.xlabel('time')
	pl.ylabel('force')
	#pl.legend()
	pl.savefig(args.OUT_img1)
	#pl.show()
	pl.close(fig1)

	fig2=pl.figure()
	pl.title("Mean cosine similarity",alpha=0.6)
	pl.scatter(xx,cosCUM,c='k',marker='o',s=1.5,alpha=0.5)
	pl.xlabel('time')
	pl.ylabel('cosine similarity')
	#pl.ylim((0,1))
	pl.savefig(args.OUT_img2)
	#pl.show()
	pl.close(fig2)

	NN=len(avex)-1
	xf= ("{0:.6f}".format(avex[NN-1]))
	yf= ("{0:.6f}".format(avey[NN-1]))
	zf= ("{0:.6f}".format(avez[NN-1]))
	print xf, yf, zf
