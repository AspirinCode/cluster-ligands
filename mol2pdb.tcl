mol load mol2 /home/mlawrenz//FKBP-dock/docking/agonist/apo-aln-stateSTATE.cent/FILE_000.mol2
set sel [ atomselect top "all"]
set outname "FILE_000"
set name [format "%s.pdb" $outname ]
$sel writepdb $name
exit
