#!/bin/bash

for dir in `ls -d state*/`
do
echo $dir
awk '{print $1}' < ranked-${dir%%/}-docking-scores.txt > files.tmp
for file in `cat files.tmp` 
do
file=`ls $dir/molecule*_${file}_000.mol2`
if [ ! -z $file ]
then
ls /home/mlawrenz/wontkill/new-results/bi-bindb-pgeomx/$file >> ranked-${dir%%/}-dbase.list
fi
done
done

