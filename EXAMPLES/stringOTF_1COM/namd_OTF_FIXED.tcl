  # The IDs of the four atoms defining the dihedral
  # atom indexes of the ligand
  set COM [addgroup {1 2 3 4 5 6}]
  set k 100 

# Output parameters

  set ofil "restr_coord.out"
  set freq  100 ;             # frequency for writing
  
# OTF parameter
  set timestep 0.02045 ;# Timestep in internal AKMA units (Check this?)
                        # 20 AKMA time units is .978 picoseconds
  set gamma 500.0
  set alpha [expr $timestep/$gamma];


  proc calcforces {} {
  
    global COM XX0 YY0 ZZ0 k ofil fout
    global alpha freq 
    
    loadcoords p

    # Calculate the current COORDS X, Y, Z.
    set XX [lindex $p($COM) 0]
    set YY [lindex $p($COM) 1]
    set ZZ [lindex $p($COM) 2]	

    set delX [expr $XX-$XX0]
    set delY [expr $YY-$YY0]
    set delZ [expr $ZZ-$ZZ0]

    # (optional) Add this constraining energy to "MISC" in the energy output
    addenergy [expr $k*$delX*$delX/2.0]
    addenergy [expr $k*$delY*$delY/2.0]
    addenergy [expr $k*$delZ*$delZ/2.0]
    
    # Calculate the "force" along the coordinates according to the harmonic constraint
    set forceX [expr -$k*($delX)]
    set forceY [expr -$k*($delY)]
    set forceZ [expr -$k*($delZ)]
    lappend force $forceX $forceY $forceZ

    # The force to be applied on each atom is proportional to its
    # corresponding gradient
    addforce $COM $force

    #set fX [expr {$k*$delX}]
    #set fY [expr {$k*$delY}]
    #set fZ [expr {$k*$delZ}]

    # Evolve phi and psi
      #set XX0 [expr $XX0+$fX*$alpha]
      #set YY0 [expr $YY0+$fY*$alpha]
      #set ZZ0 [expr $ZZ0+$fZ*$alpha]

    #printout

    if {[expr [getstep] % $freq] == 0} {
        set fout [open $ofil a+]
        puts $fout [format "%f %f %f %f" [getstep] $XX0 $YY0 $ZZ0]
        close $fout
    }
   
   return

  }


