proc run {} {
#bigdcd getCOM rep7_frontal_6315_12000.dcd
#bigdcd rmsdCAi rep7_frontal_6315_12000.dcd
bigdcd getCOM p5.dcd
bigdcd_wait
}

proc caca { frame } {
}

proc rmsdCAi { frame } {
	# fiting before rmsd calculation. 
	# totalfit=1 (fit the selection 'selprot', i.e, the fit is for the complete protein)
	# totalfit=0 (fit the selection 'selprot' and the range of residues between res1 and 'resN') 
	set totalfit 0  
        set res1 1
        set resN 284
	set selprot "backbone"

	if {$totalfit == 1} {
		set ff $selprot
		set ref [atomselect top $selprot frame 0]
	} else {
		set w3 " to "
		set w1 " and "
		set w2 "resid "
		set ss $selprot$w1$w2$res1$w3$resN 
		set ff $ss
                set ref [atomselect top $ss frame 0]
        }


        #Perform the fit
        set all [atomselect top all]
	set sel [atomselect top $ff]
	$all move [measure fit $sel $ref]
	
        set R1 "resid "
	set space " "
	set RT ""
        for {set i $res1} {$i<=$resN} {incr i} {
                set R2 $R1$i
		set dd $selprot$w1$R2
		set ref [atomselect top $dd frame 0]
		set sel [atomselect top $dd]
		set tmp [measure rmsd $sel $ref]
		if {$i == $res1} {
			append rmsdi [format "%.5f" $tmp]
		} else {
			append rmsdi " " [format "%.5f" $tmp]
		}
        }
	#puts $rmsdi

        set outfile [open rmsd_backbone_all.dat a]
	puts $outfile $rmsdi
        close $outfile
}

proc getCOM {frame} {
# get coordinates of the center of mass for the selection
# to use see README.txt
set molid [molinfo top]
        set com_coord [COM_COORD top "resname DIH and resid 863"]
	set outfile [open COM.txt a]
        #set outfile [open LIGAND_COM_ALL.txt a]
	puts $outfile $com_coord
	close $outfile


}

proc COM_COORD {molid1 stxt1} {
        set s1 [atomselect $molid1 "$stxt1"]
        set v1 [measure center $s1]
        $s1 delete
        return $v1
}

proc rmsdCA { frame } {
	#calculate total rmsd for CA atoms
	set all [atomselect top all]
	set ref [atomselect top "name CA" frame 0]
	set sel [atomselect top "name CA"]
	####

	$all move [measure fit $sel $ref]
	set rmsdCA [measure rmsd $sel $ref]
	set outfile [open rmsdCA.txt a]
	puts $outfile $rmsdCA
	close $outfile
}

#####################################
####### Functions the bigdcd.tcl ####
##########################################################################
# the previous functions presented in this collection make reference to 
# this code. We have include all the functions in a single file to 
# facilitate the use. The code from bigdcd was not modified.
###########################################################################
###########################################################################  
# BigDCD-v2
# Justin Gullingsrud
# updates by Axel Kohlmeyer
# vmd@ks.uiuc.edu
################################################################################
# Purpose: Use this script to analyze one or more trajectory files that 
# don't fit into memory.  The script will arrage for your analysis function 
# to be called each time a frame is loaded, then delete the frame from memory.
# The analysis script must accept one argument; BigDCD will keep track of how
# many timesteps have been loaded and call your script with that number.
#
# How to include this function: either source the script directory, or 
# (better) place the script in one of the directories in your auto_path 
# variable and include "package require bigdcd" in your script.
#
# New in version 2:
#
# The bigdcd command accepts an (optional) argument that defines
# the format type of the trajectory file. Its default is "auto", i.e.
# let VMD guess it from the file name.
#
# The function bigdcd_wait can be called to have a script wait until
# all frames have been processed. This is most useful for batch processing
# of analysis script. By construction bigdcd executes in the background 
# and will return the control to VMD after the last trajectory file has 
# been scheduled for reading. Thus analysis scripts will terminate prematurely.
################################################################################

# Example 1: 
# This computes the center of mass for each frame in the DCD file.
#
# proc mycenter { frame } {
#   global all
#   puts "$frame: [measure center $all weight mass]"
# }
# set $mol [mol new alanin.psf type psf waitfor all]
# set all [atomselect $mol all]
# $all global
# bigdcd mycenter auto alanin.dcd 
#

# Example 2:
# This computes the RMS distance between each frame in
# a sequence of xyz files and a reference pdb file.  
# this example even works in batch mode.
#
# proc myrmsd { frame } {
#   global ref sel all
#   $all move [measure fit $sel $ref]
#   puts "$frame: [measure rmsd $sel $ref]"
# }
# set mol [mol new protein.psf type psf waitfor all]
# set all [atomselect $mol all]
# set ref [atomselect $mol "name CA" frame 0]
# set sel [atomselect $mol "name CA"]
# mol addfile protein.pdb waitfor all
# bigdcd myrmsd xyz eq01.xyz eq02.xyz eq03.xyz
# bigdcd_wait
# quit

proc bigdcd { script type args } {
    global bigdcd_frame bigdcd_proc bigdcd_firstframe vmd_frame bigdcd_running
  
    set bigdcd_running 1
    set bigdcd_frame 0
    set bigdcd_firstframe [molinfo top get numframes]
    set bigdcd_proc $script

    # backwards "compatibility". type flag is omitted.
    if {[file exists $type]} { 
        set args [linsert $args 0 $type] 
        set type auto
    }
  
    uplevel #0 trace variable vmd_frame w bigdcd_callback
    foreach dcd $args {
        if { $type == "auto" } {
            mol addfile $dcd waitfor 0
        } else {
            mol addfile $dcd type $type waitfor 0
        }
    }
    after idle bigdcd_wait
}

proc bigdcd_callback { tracedvar mol op } {
    global bigdcd_frame bigdcd_proc bigdcd_firstframe vmd_frame
    set msg {}
 
    # If we're out of frames, we're also done 
    # AK: (can this happen at all these days???). XXX
    set thisframe $vmd_frame($mol)
    if { $thisframe < $bigdcd_firstframe } {
        puts "end of frames"
        bigdcd_done
        return
    }
 
    incr bigdcd_frame
    if { [catch {uplevel #0 $bigdcd_proc $bigdcd_frame} msg] } { 
        puts stderr "bigdcd aborting at frame $bigdcd_frame\n$msg"
        bigdcd_done
        return
    }
    animate delete beg $thisframe end $thisframe $mol
    return $msg
}

proc bigdcd_done { } {
    global bigdcd_running
    
    if {$bigdcd_running > 0} then {
        uplevel #0 trace vdelete vmd_frame w bigdcd_callback
        puts "bigdcd_done"
        set bigdcd_running 0
    }
}

proc bigdcd_wait { } {
    global bigdcd_running bigdcd_frame
    while {$bigdcd_running > 0} {
        global bigdcd_oldframe
        set bigdcd_oldframe $bigdcd_frame
        # run global processing hooks (including loading of scheduled frames)
        display update ui
        # if we have read a new frame during then the two should be different.
        if { $bigdcd_oldframe == $bigdcd_frame } {bigdcd_done}
    }
}
