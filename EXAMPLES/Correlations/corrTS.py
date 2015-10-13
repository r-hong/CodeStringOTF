#! /usr/bin/env python
'''
'''
###################################################
# Rolando P. Hong Enriquez (rolando.hong@iit.it) ##
###################################################
import argparse
import csv
import numpy as np
import statsmodels.tsa.stattools as ts
import pandas as pd
import matplotlib.pyplot as plt
from functions import file_cols,plotCorrMatrix
#######################################################
###### Parsing arguments from the command line ########
#######################################################
def _make_parser():
	parser = argparse.ArgumentParser(description="This tool explores the correlations between time series. The input file should contain a time series on each column. Initially the data is loaded in a panda dataframe and later the correlation matrix is calculated and visualized.",epilog="Thanks for using corrTS.py!!")
	parser.add_argument('IN_file', help="Input file")
        return parser
#------------------------------
_p    = _make_parser()
__doc__ += _p.format_help()
######################################################
if __name__ == '__main__':
        args = _p.parse_args()

	#Read filename
	nc = file_cols(args.IN_file)
	for i in range(nc):
        	exec "c%s=[]" % (i)
	with open(args.IN_file) as file:
        	reader = csv.reader(file, delimiter=' ',skipinitialspace=True)
        	for row in reader:
			for i in range(nc):
                		exec "c%s.append(float(row[%s]))" % (i,i)
	allC=[]
	for i in range(nc):
		exec "allC.append(c%s)" % (i)
	#put everything in a panda data frame
	matrix=np.asmatrix(allC)
	matrixT=matrix.T
	df = pd.DataFrame(data=matrixT)

	plotCorrMatrix(df,10)


