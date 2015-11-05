#! /usr/bin/env python
'''
'''
import argparse
from functions import *
from scipy.stats import kde
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

def _make_parser():
	############################################################################################
	p = argparse.ArgumentParser(description="This tool has a limited scope. It was designed with the idea of analyzing the differences between the monomers of PNP along the binding path of the ligand but its use can be extended in the folowing way. Use as input a file with N columns that can be divided into N/3 equivalent parts or monomers A, B, and C. The divisions are such that column i in monomer A is equivalent to columns i+ N/3 and i+2N/3 in monomers B and C respectively. This tool (1) uses kernel density estimation (KDE) to obtain the probability density distribution of each column, (2) analyses the differences between the KDE probability distributions of each set of equivalent columns; that is, for equivalent columns Aj, Bj, and Cj we calculate the distances d(Aj,Bj), d(Aj,Cj), and d(Bj,Cj). (3) next we determine the 'excess' in which a distance dominates the others. In practice, to see how far is Aj from the others we calculate Q_Aj = mean([d(Aj,Bj), d(Aj,Cj)]) - mean([d(Aj,Bj), d(Aj,Cj), d(Bj,Cj) ]), if mean([d(Aj,Bj), d(Aj,Cj)]) > mean([d(Aj,Bj), d(Aj,Cj), d(Bj,Cj) ])  (4) we also sum all Q_Aj, Q_Bj, and Q_Cj to have an estimate of the total excess of monomers A, B and C, respectively. The results are shown by column and by monomer for all the provided data files or images.",epilog="Thanks for using trimerKDE.py!!")
	p.add_argument('prefix', help="prefix of the input files.")
	p.add_argument('monoSize', type=int, help="Number of columns ocupied by a single unit (monomer). In this specific case for trimeric proteins this number should be 1/3 of the total number of columns in the input files.")
	p.add_argument('Nimg', type=int, help="Number of images. Each image should be in an independent file. E.g., if the prefix is 'file' then the data for the first image is file1.")
        return p
#--------------------------------
_p    = _make_parser()
__doc__ += _p.format_help()
######################################################
if __name__ == '__main__':
        args = _p.parse_args()
######################################################
	gridSize=100
	pname=args.prefix
	CVA=[]
	CVB=[]
	CVC=[]
	sumA=[]
	sumB=[]
	sumC=[]
	for q in range(1,args.Nimg+1):
		args.prefix +=str(q)
		N=file_lines(args.prefix)
		M=file_cols(args.prefix)
		P1=load_data_range(args.prefix,N,0,M-1)
		P2=np.asarray(P1)
		PZ=zip(*P2)
		args.prefix=pname

		betaA=[]
		betaB=[]
		betaC=[]
		for i in range(args.monoSize):
			ia=i
			ib=(i+args.monoSize)
			ic=(i+(2*args.monoSize))
		
			xa,xb,xab=getVec(PZ,ia,ib)
			xa,xc,xac=getVec(PZ,ia,ic)
			xb,xc,xbc=getVec(PZ,ib,ic)

			#common area A,B
			density1 = kde.gaussian_kde(xa)
			xgrid1 = np.linspace(xa.min(), xa.max(), gridSize)
			density2 = kde.gaussian_kde(xb)
			xgrid2 = np.linspace(xb.min(), xb.max(), gridSize)
			densityt = kde.gaussian_kde(xab)
			xgridt = np.linspace(xab.min(), xab.max(), gridSize)
			ac, D_ab = commonArea(density1,density2,xgridt)
	
	        	#common area A,C
	        	density1 = kde.gaussian_kde(xa)
	        	xgrid1 = np.linspace(xa.min(), xa.max(), gridSize)
	        	density2 = kde.gaussian_kde(xc)
	        	xgrid2 = np.linspace(xc.min(), xc.max(), gridSize)
	        	densityt = kde.gaussian_kde(xac)
	        	xgridt = np.linspace(xac.min(), xac.max(), gridSize)
	        	ac, D_ac = commonArea(density1,density2,xgridt)
	
	        	#common area B,C
	        	density1 = kde.gaussian_kde(xb)
	        	xgrid1 = np.linspace(xb.min(), xb.max(), gridSize)
	        	density2 = kde.gaussian_kde(xc)
	        	xgrid2 = np.linspace(xc.min(), xc.max(), gridSize)
	        	densityt = kde.gaussian_kde(xbc)
	        	xgridt = np.linspace(xbc.min(), xbc.max(), gridSize)
	        	ac, D_bc = commonArea(density1,density2,xgridt)
			dmean=np.mean([D_ab, D_ac, D_bc])
			dstd=np.std([D_ab, D_ac, D_bc])

			#betaA: How much the average distance to monomer A (from B and C) exceeds the average distance detween A,B, and C.
			#betaB and betaC defined in the same way. This definition is based on each residue (column)
			if (np.mean([D_ab, D_ac]) > np.mean([D_ab, D_ac, D_bc])):
				betaA.append((np.mean([D_ab, D_ac])-np.mean([D_ab, D_ac, D_bc])))
			else:
				betaA.append(0)

	        	if (np.mean([D_ab, D_bc]) > np.mean([D_ab, D_ac, D_bc])):
	                	betaB.append((np.mean([D_ab, D_bc])-np.mean([D_ab, D_ac, D_bc])))
	        	else:
	                	betaB.append(0)

	        	if (np.mean([D_ac, D_bc]) > np.mean([D_ab, D_ac, D_bc])):
	                	betaC.append((np.mean([D_ac, D_bc])-np.mean([D_ab, D_ac, D_bc])))
	        	else:
	                	betaC.append(0)

		sa=np.std(betaA)
		sb=np.std(betaB)
		sc=np.std(betaC)
		print 'image: ',q
		print 'mean +/- std'
		print 'A:',np.mean(betaA),'+/-',np.std(betaA)
		print 'B:',np.mean(betaB),'+/-',np.std(betaB)
		print 'C:',np.mean(betaC),'+/-',np.std(betaC)

		#sumA: is the sum of all the instances in which the distance to a residue (column) in A, exceeded the equivalent residues in B and C.
		#therefore sumA is the total excess of distance to A. The higher sumA the greater the total difference of A with respect to Band C.  
		print 'sum of all the differences'
		print 'A:',sum(betaA)
		print 'B:',sum(betaB)
		print 'C:',sum(betaC)
		CVA.append(betaA)
		CVB.append(betaB)
		CVC.append(betaC)
		sumA.append(sum(betaA))
		sumB.append(sum(betaB))
		sumC.append(sum(betaC))

	ZZA=zip(*CVA)
	ZZB=zip(*CVB)
	ZZC=zip(*CVC)
	zA=np.asarray(ZZA)
	zB=np.asarray(ZZB)
	zC=np.asarray(ZZC)

	###############################################
	#Creating 2 2D grids for the x and y bounds
	##############################################
	dx, dy = 1, 1
	y, x = np.mgrid[slice(1, args.monoSize + dy, dy),
	                slice(1, args.Nimg + dx, dx)]


	###################################
	### Contour  mean A ################
	###################################
	zA = zA[:-1, :-1]
	levels = MaxNLocator(nbins=15).tick_values(zA.min(), zA.max())
	cmap = plt.get_cmap('PiYG')
	norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
	plt.contourf(x[:-1, :-1] + dx / 2.,
	             y[:-1, :-1] + dy / 2., zA, levels=levels,
	             cmap=cmap)
	plt.colorbar()
	plt.title('RMSD change per residue along the binding trajectory (monomer A)')
	plt.show()

	###################################
	### Contour  mean B ###############
	###################################
	zB = zB[:-1, :-1]
	levels = MaxNLocator(nbins=15).tick_values(zB.min(), zB.max())
	cmap = plt.get_cmap('PiYG')
	norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
	plt.contourf(x[:-1, :-1] + dx / 2.,
	             y[:-1, :-1] + dy / 2., zB, levels=levels,
	             cmap=cmap)
	plt.colorbar()
	plt.title('RMSD change per residue along the binding trajectory (monomer B)')
	plt.show()

	###################################
	### Contour  mean C ###############
	###################################
	zC = zC[:-1, :-1]
	levels = MaxNLocator(nbins=15).tick_values(zC.min(), zC.max())
	cmap = plt.get_cmap('PiYG')
	norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
	plt.contourf(x[:-1, :-1] + dx / 2.,
	             y[:-1, :-1] + dy / 2., zC, levels=levels,
	             cmap=cmap)
	plt.colorbar()
	plt.title('RMSD change per residue along the binding trajectory (monomer C)')
	plt.show()

	#####################
	###Final Bar plot ###
	#####################
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ind = np.arange(args.Nimg)        # x locations for bars
	width = 0.35                      # width of the bars
	#non-stacked bars
	rects1 = ax.bar(ind,           sumA, width/2.0, color='yellow')
	rects2 = ax.bar(ind+width/2.0, sumB, width/2.0, color='green' )
	rects3 = ax.bar(ind+width,     sumC, width/2.0, color='blue'  )
	# axes and labels
	ax.set_xlim(-2*width,len(ind)+width)
	ax.set_ylabel('Total Excess')
	ax.set_title('Total exceess over the other monomers per image')
	xTickMarks = ['Image'+str(i) for i in range(1,args.Nimg+1)]
	ax.set_xticks(ind+width)
	xtickNames = ax.set_xticklabels(xTickMarks)
	plt.setp(xtickNames, rotation=45, fontsize=10)
	ax.legend( (rects1[0], rects2[0], rects3[0]), ('monomer A', 'monomer B', 'monomer C') )
	plt.show()
