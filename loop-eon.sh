#!/bin/bash

n=1
for file in `cat eon-dbase.list`
do
sed "s/NNN/${n}/g" < eon-only-template.parm  | sed "s/MOLECULE/${file}/g" > eon-only-molecule${n}.parm
eon -param eon-only-molecule${n}.parm 
let n=($n + 1)
done
