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
# Loop on iterations
for ((it=1; it<=$Nit; it++))
do
	save_coord=`echo $it $Nsave | awk '{save_coord=0;if($1 % $2 == 0) save_coord=1;print save_coord}'`	
	echo "   Iteration $it .. "
	# prepare the namd input file
	img=string0.dat
	mdstep=100                                # No. of MD steps on each iteration
	COORDXX=($(awk '{print $1}' $img))
	COORDYY=($(awk '{print $2}' $img))
	COORDZZ=($(awk '{print $3}' $img))
	# Loop on the pooints of the string (to prepare namd configuration file with CV values of previous iteration)
	# in this version the first point on the string (unbounded ligand) is keep fixed
	for ((i=$p1; i<=$pN; i++))
	do
                if ([ $i -eq 1 ])
                then
	                cp TMP_COORD_FIXED.namd COORD$i.namd
        	        echo "set XX0 ${COORDXX[$i-1]}" >> COORD$i.namd
                	echo "set YY0 ${COORDYY[$i-1]}" >> COORD$i.namd
                	echo "set ZZ0 ${COORDZZ[$i-1]}" >> COORD$i.namd
         	        echo "" >> COORD$i.namd
                	echo "run     100" >> COORD$i.namd
		else

			cp TMP_COORD.namd COORD$i.namd
			echo "set XX0 ${COORDXX[$i-1]}" >> COORD$i.namd
			echo "set YY0 ${COORDYY[$i-1]}" >> COORD$i.namd
			echo "set ZZ0 ${COORDZZ[$i-1]}" >> COORD$i.namd
			echo "" >> COORD$i.namd
			echo "run     100" >> COORD$i.namd
		fi
	done

	# Loop on the points of the string (to compute mean forces and evolve them)
	for ((i=$p1; i<=$pN; i++))
	do
		# namd coords input (always indexed *-0.coor)
                cp OUTPUT/namdout$i-0.coor T.coor
                cp OUTPUT/namdout$i-0.xsc T.xsc
                cp OUTPUT/namdout$i-0.vel T.vel
                cp OUTPUT/namdout$i-0.top T.top

		# run restrained MD
		echo "running namd for image $i ..."
		$namd +idlepoll +p6 COORD$i.namd > OUTPUT/namdout-$i

		# saving output files
		if [[ $save_coord -eq 1 ]]
		then
        		echo "save"
                        mv restr_coord.out OUTPUT/COORD$i.out
                        cp namdout.coor OUTPUT/namdout$i-0.coor
                        mv namdout.coor OUTPUT/namdout$i-$it.coor
                        cp namdout.xsc OUTPUT/namdout$i-0.xsc
                        mv namdout.xsc OUTPUT/namdout$i-$it.xsc
                        cp namdout.vel OUTPUT/namdout$i-0.vel
                        mv namdout.vel OUTPUT/namdout$i-$it.vel
                        cp namdout.colvars.traj OUTPUT/namdout$i-$it
		else
        		echo "don't save"
                        mv restr_coord.out OUTPUT/COORD$i.out
                        cp namdout.coor OUTPUT/namdout$i-0.coor
                        cp namdout.xsc OUTPUT/namdout$i-0.xsc
                        cp namdout.vel OUTPUT/namdout$i-0.vel
                        cp namdout.colvars.traj OUTPUT/namdout$i-$it.cv
		fi

		# read latest CV values
		XX=`tail -1 OUTPUT/COORD$i.out | awk '{print $2}'`
		YY=`tail -1 OUTPUT/COORD$i.out | awk '{print $3}'`
		ZZ=`tail -1 OUTPUT/COORD$i.out | awk '{print $4}'`
		echo "$XX $YY $ZZ" >> string.dat
	done
	# reparametrization
	mv string.dat images_norep.ls
	./reparametrization.x
	cp images_rep.ls string0.dat
	mv images_rep.ls string$it.dat
done

