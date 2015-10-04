#! /usr/bin/env python
"""
Finds the maximum likelihood estimator of each column in a file. This version processes automatically files with any number of columns and the estimator is calculated using a bootstrap procedure. 
"""

###################################################
# Rolando P. Hong Enriquez (rolando.hong@iit.it) ##
###################################################
import argparse
import os
import csv
import numpy as np
import math
import random
import scipy as sp
import scipy.stats
from functions import file_cols,readFile,bootstrapN,estValB,formatRes

def _make_parser():
	parser = argparse.ArgumentParser(description="Finds the maximum likelihood estimator of each column in a file. This version processes automatically files with any number of columns and the estimator is calculated using a bootstrap procedure.",epilog="Thanks for using maxLikeBoot.py!!")
	parser.add_argument('fileName', help="Input file")
	parser.add_argument('N', type=int, help="Size of the resampled vector (on each column).")
	parser.add_argument('K', type=int, help="Number of times the resampling of size N will take place.")
	parser.add_argument('S', type=int, help="Stride. Use only every Nth value from the input file.")
	#args = parser.parse_args()
	return parser
#------------------------------
_p    = _make_parser()
__doc__ += _p.format_help()
######################################################
if __name__ == '__main__':
        args = _p.parse_args()

	#preprocessing the input file:
	# saving every args.stride value of args.IN_file and saving to tmp.dat
	ss="'"
	ss+='0~'
	ss+=`args.S`
	ss+='p'
	ss+="'"

	cmd='sed -n '
	cmd+=ss
	cmd+=' '
	cmd+=args.fileName
	cmd+=' > tmp.dat'
	os.system(cmd)
	##########################################
	#get Number of rows and number of columns from the file and read from the file into a matrix 'vecT'
	nCol = file_cols('tmp.dat')
	vecT = readFile('tmp.dat',nCol)


	#print len(vecT[0])

	#calculate the max likelihood estimators
	results = estValB(vecT,args.K,args.N)
	#print results[0]
	#format the results to 6 decimal points before printing
	myList=formatRes(results,6)
	print ' '.join(str(p) for p in myList)
