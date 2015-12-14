# The force depends on the distance between the COMs of 2 groups
# the force is applied on group COM1_end in the ditection of COM2_end
  set COM1_end [addgroup {1 2 3}]
  set COM2_end [addgroup {2092 2093 2094}]
  set freq 100	


# OTF parameter
set timestep 0.02045 ;# Timestep in internal AKMA units (Check this?)
                        # 20 AKMA time units is .978 picoseconds
set gamma 100.0
set alpha [expr $timestep/$gamma];

 proc calcforces {} {
	global COM1_end  COM2_end COM1_ini COM2_ini lastCOMDIST 
	global k diff diffT alpha gamma timestep ts freq 
	loadcoords p
	set ts [getstep]
 
	#load COMDIST
	#set fp [open "lastCOMDIST.dat" r]
	#set lastCOMDIST [read $fp]
	#close $fp

	#print "coordinates at the begining of the procedure $p($COM_end)"
	set diff [vecsub $p($COM1_end) $p($COM2_end)]
	set diffT [expr abs([lindex $diff 0]) + abs([lindex $diff 1]) + abs([lindex $diff 2])]
	set gradient [list [expr abs([lindex $diff 0])] [expr abs([lindex $diff 1])] [expr abs([lindex $diff 2])]]

	# Calculate the "force" along the COM distance according to the harmonic constraint
	set force [expr -$k*($diffT)]

	# The force to be applied on each atom is proportional to its
	# corresponding gradient
	addforce $COM1_end [vecscale $diff $force]
	#addforce $COM2_end [vecscale -$gradient $force]

	set f1 [expr {$k*$diffT}]

	# Evolve COM
	set lastCOMDIST [expr {$lastCOMDIST+$f1*$alpha}]

	#print "coordinates at the end of the procedure $p($COM_end)"
	if {[expr [getstep] % $freq] == 0} {

		#save last COM1 and coordinates 
          	set outfile [open "lastCOMDIST.dat" w]
        	puts $outfile [format "%f" $diffT]
        	close $outfile

        	set trackfile [open "trackCOMDIST.dat" a]
        	puts $trackfile [format "%f" $lastCOMDIST]
        	close $trackfile

	}
}
