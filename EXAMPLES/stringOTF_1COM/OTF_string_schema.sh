#!/bin/bash
###########################################################################################
#####  Initializations  ###################################################################
###########################################################################################
p1=1                                          # first point on the string
pN=3                                          # last point on the string
Nit=1000                                      # No. of iterations
Nsave=10                                      # save coordinates every Nsave iterations
LD_LIBRARY_PATH="/path1/:$LD_LIBRARY_PATH";   # path to namd library (tobe used with CUDA) 
export LD_LIBRARY_PATH
namd='/path2/namd2'                           # path to namd executable 

####################################################################################
###  comment the following lines if this run is a continuation of a previous one ###
cp INPUT/* OUTPUT/                                                               ###
cp INPUT/initialPcurve.dat string0.dat        # values of the initial string     ###
####################################################################################
####################
###  START      ####
####################

# Loop on the OTF iterations
for ((i=1; i<=$N; i++))
do
	>Get the CVs from old string (previous iteration)
	# Loop on the points of the string
	for ((j=1; j<=$M; j++))
	do
		>Prepare namd configuration file
		>Get coordinates, velocities, box size from the previous iteration
		>run MD to compute and evolve forces on CVs
		>save final coordinates, velocities and box size
		>save final values of the CVs
	done
	>reparametrize new string
done

