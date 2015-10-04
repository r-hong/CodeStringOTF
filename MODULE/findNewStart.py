#! /usr/bin/env python
'''
This is my script
'''
###################################################
# Rolando P. Hong Enriquez (rolando.hong@iit.it) ##
###################################################
import pydoc
import os
import sys
import csv
import pandas as pd
import argparse
from functions import pdb2df, df2pdb, cleanpdb, ligCOM, getPdbFilename, getMeanString, dm_vec, reparam,getX2near,getX2far
import numpy as np
import shutil
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

#######################################################
###### Parsing arguments from the command line ########
#######################################################
def _make_parser():
	parser = argparse.ArgumentParser(description="This uses coordinate data produced by a OTF string method calculation to find the mean string from the strings produced between iterations Nit1 and Nit2, then reparametrizes the found mean string to Mout points and finally fonds the coordinate files that are closest to each point in the Mout reparametrized string.",epilog="Thanks for using findNewStart.py!!")
	parser.add_argument('ligand', help="name of the ligand residue in the pdb file.")
	parser.add_argument('Nimages', type=int, help="Number of images in the string.")
	parser.add_argument('Nit1', type=int, help="Fist iteration to be analyzed.")
	parser.add_argument('Nit2', type=int, help="Last iteration to be analized.")
	parser.add_argument('Mout', type=int, help="Number ot points of the final reparametrized string.")
	return parser
#--------------------------------
_p    = _make_parser()
__doc__ += _p.format_help()

######################################################
if __name__ == '__main__':
        args = _p.parse_args()
	print '#######################################################'
	print '####################  Help  ###########################'
	print 'This script uses data from a string method calculation to:'
	print '(1) find the mean string from a group of strings between iteration Nit1 and Nit2'
	print '(2) reparametrizes the found mean string to Mout equidistant points'
	print '(3) finds the pdb files (using allways the data from Nit1 to Nit2) that are closest'
	print '    to each point in the Mout reparametrized string'
	print '(4) saves: (a) the coordinates of the center of mass (COM) of the new Mout points (CVs), '
	print '           (b) the indices of the pdb files between Nit1 and Nit2 used to create the new starting pdb files' 
	print '           (c) starting pdb files (Pi.pdb, 1<=i<=Mout) '
	print 'a new run will move the ligand from their initial positions in (c) to the positions specified in (a)     ' 
	print '########################################################'

	#Initial conditions
	pdb_tmp='temp.pdb'
	DIM=3

	#Put all strings in a single variable
	allstrings=[]
	existingF=[]
	for i in range(args.Nit1,args.Nit2):
		filea=getPdbFilename(1,i)
		string=[]
		if os.path.isfile(filea):
			#existingF.append(i) #existing files
			for j in range(1,args.Nimages+1):
				pdb_A=getPdbFilename(j,i)
				shutil.copy2(pdb_A, pdb_tmp)
				df=pdb2df(pdb_tmp)
				xc,yc,zc=ligCOM(df,args.ligand)
				string.append([xc,yc,zc])
		else:
			a=1
		if (not string):
			a=1
		else:
			allstrings.append(string)

	#get the mean string from a group of strings
	v1,v2,v3=getMeanString(allstrings)
	vecn=[v1,v2,v3]
	vecz=zip(*vecn)
	vec_array=np.asarray(vecz)

	#reparametrize the mean string
	meanD,totD=dm_vec(vecz,args.Mout)
	vecr1=reparam(vec_array,DIM,args.Mout,meanD)
	vecr2=zip(*vecr1)
	#############################################################
	#Saving the reparametrized center of mass (COM) coordinates
	##formating the COMs###
	x=vecr2[0]
	y=vecr2[1]
	z=vecr2[2]
	x1=[]
	y1=[]
	z1=[]
	for i in range(args.Mout):
		a= ("{0:.4f}".format(float(x[i])))
		b= ("{0:.4f}".format(float(y[i])))
		c= ("{0:.4f}".format(float(z[i])))
		x1.append(a)
		y1.append(b)
		z1.append(c)
	#putting the COMs in a dataframe
	dcom={'x': x1,
	      'y': y1,
	      'z': z1}
	dfcom = pd.DataFrame(dcom)
	#saving the COMs from the dataframe
	outf='COMcoords.txt'
	f = open(outf, "w")
	for i in xrange(args.Mout):
		f.write("{:<7} {:<7} {:<7}\n".format(x1[i], y1[i], z1[i]))
	f.close()

	################################################################
	# getting the closest pdb files to each image in the reparametrized (Mout) string
	indices=[] #indices (j,i) of the selected pdb files (namdout[j]-[i].coor)
	COMs=[]    #COM of the ligand in the selected pdb files
	for q in range(args.Mout):
		dist_final=5
		for i in range(args.Nit1,args.Nit2):
        		filea=getPdbFilename(1,i)
        		string=[]
        		if os.path.isfile(filea):
                		for j in range(1,args.Nimages+1):
                        		pdb_A=getPdbFilename(j,i)
                        		shutil.copy2(pdb_A, pdb_tmp)
                        		df=pdb2df(pdb_tmp)
                        		xc,yc,zc=ligCOM(df,args.ligand)
					ctemp=np.asarray([xc,yc,zc])
					#measure distance here
					d=np.linalg.norm(vecr1[q]-ctemp)
					if (d<dist_final):
						dist_final=d
						ind_final=[j,i]
						COM_final=[xc,yc,zc]
		indices.append(ind_final)
		COMs.append(COM_final)

	#saving the Mout pdb files with the coordinates
	for q in range(args.Mout):
		d=indices[q]
		ss=int(q+1)
		#print 'ss',ss
		j=d[0]
		i=d[1]
		pdb_A=getPdbFilename(j,i)
		prefix='P'
		prefix+=`ss`
		prefix+='.pdb'
		shutil.copy2(pdb_A, prefix)
		prefix=''

	#Saving the indices of the pdb files that were used to build the Mout coordinate pdb files
	outc='pdb_indices.txt'
	g = open(outc, "w")
	g.write("Indices (a,b) of the pdb files used to build a string with (Mout={:<3}) images\n".format(args.Mout))
	g.write("namdout[a]-[b].coor\n")
	g.write("Iterations taken into account for this calculation:\n")
	g.write("Initial iteration = {:<3}, Final iteration = {:<3}\n".format(args.Nit1,args.Nit2))
	g.write("The COM coordinates of the reparametrized string calculated from the mean string between the iterations are saved in the file {:<20}\n".format(outf))
	for q in xrange(args.Mout):
		d=indices[q]
		j=d[0]
		i=d[1]
        	g.write("{:<3} {:<3}\n".format(j, i))
	g.close()

