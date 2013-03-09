mol load mol2 /home/mlawrenz//FKBP-dock/docking/agonist/apo-aln-state2665.cent/stereo-44209747-0_000.mol2
set sel [ atomselect top "all"]
set outname "stereo-44209747-0_000"
set name [format "%s.pdb" $outname ]
$sel writepdb $name
exit
