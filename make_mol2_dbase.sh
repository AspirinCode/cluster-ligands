#!/bin/bash

#list='182176 2083'

duplicates='16738803 44211748 19826083 182176 2083 7028'	 
count=0
state=$1
rm -f eon-dbase.list
rm -f eon-dbase.mol2
for file in `cat reduce-* | awk '{print $1}'`
do
tmp=${file##*stereo-}
tmp2=${tmp%%-*}
if [[ "$duplicates" != *$tmp2* ]]
then
cat /home/mlawrenz//FKBP-dock/docking/agonist/apo-aln-state${state}.cent/${file}_000.mol2 >> eon-dbase.mol2
echo /home/mlawrenz//FKBP-dock/docking/agonist/apo-aln-state${state}.cent/${file}_000.mol2 >> eon-dbase.list
fi
done
