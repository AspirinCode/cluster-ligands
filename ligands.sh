#!/bin/bash

count=0
state=$1
rm -f eon-dbase.list
rm -f eon-dbase.pdb
for file in `cat reduce-* | awk '{print $1}'`
do
sed "s/FILE/${file}/g" < mol2pdb.tcl | sed "s/STATE/${state}/g" > tmp.tcl
vmd -dispdev text -e tmp.tcl
sed "s/NAME/${file}/g" < header.txt | sed "s/XXX/${count}/g" > tmp
cat ${file}_000.pdb >> tmp
mv tmp eon-${file}.pdb
echo eon-$file >> eon-dbase.list
cat eon-$file.pdb >> eon-dbase.pdb
let count=($count + 1)
done
