#!/bin/bash

n=1
for file in `cat eon-dbase.list`
do
echo $file > tmp
newfile=`sed -f slash.sed < tmp`
sed "s/NNN/${n}/g" < eon-only-template.parm  | sed "s/MOLECULE/${newfile}/g" > eon-only-molecule${n}.parm
eon -param eon-only-molecule${n}.parm 
let n=($n + 1)
done
