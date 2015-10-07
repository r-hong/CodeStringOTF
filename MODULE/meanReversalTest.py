#! /usr/bin/env python
'''
'''
###################################################
# Rolando P. Hong Enriquez (rolando.hong@iit.it) ##
###################################################
import argparse
import os
import csv
import pylab as pl
import numpy as np
import statsmodels.tsa.stattools as ts
from functions import file_cols

#######################################################
###### Parsing arguments from the command line ########
#######################################################
def _make_parser():
	parser = argparse.ArgumentParser(description="This tool performs a mean reversion test, that is, the tendency of a time series to drift towards its long-term mean. As input file this tool expects a time series for each column; any number of columns and rows can be processes. The use of the test here was inspired by applications in financial mathematics; in our implementation the Augmented Dickey-Fuller test (ADFT) was used. With ADFT we will optain for each columns a '1' ---> mean reversal (stationarity) or '0' ---> NO time reversal (no stationarity) of the time series. Further explanations and examples can be found at https://gist.github.com/r-hong/82515ca44fb35090b8c8 ",epilog="Thanks for using meanReversalTest.py!!")
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

	#perform test
	for i in range(nc):
		exec "t%s  = ts.adfuller(c%s,1)" % (i,i)

	final=[]
	a="'1%'"
	b="'5%'"
	c="'10%'"
	dd=[]
	for i in range(nc):
		exec "tmp=t%s[4]" % (i)
		exec "a%s=( (t%s[0] < tmp[%s]) and (t%s[0] < tmp[%s]) and (t%s[0] < tmp[%s]) ) " % (i,i,a,i,b,i,c)
		exec "dd.append(a%s)" % (i)	
		if dd[i]==True:
        		final.append(1)
		else:
        		final.append(0)
	#print '1 --- mean reversion'
	#print '0 --- NO mean reversion'
	print ' '.join(str(p) for p in final)
