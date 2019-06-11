#if { $argc != 2 } {
#    puts "ERROR: PDBNAME and DATASET arguments required."
#} else {
#    set pdbname [lindex $argv 0]
#    set dataset [lindex $argv 1]
#}

set pdbname "10gs"
set dataset "refined"

mol delete all
display resetview

mol new ${dataset}/${pdbname}/dock.pdb
mol modstyle 0 0 Licorice 0.30 12.0 12.0

mol new ${dataset}/${pdbname}/flex.pdb
mol modstyle 0 1 Licorice 0.30 12.0 12.0
mol modcolor 0 1 ColorID 1

mol new ../PDBbind18/${dataset}/${pdbname}/${pdbname}_protein.pdb
mol modstyle 0 2 NewCartoon 0.30 10.0 4.1 0
mol modcolor 0 2 ColorID 6

axes location off
