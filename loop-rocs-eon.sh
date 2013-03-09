#!/bin/bash

n=1
for file in `cat dbase.list`
do
sed "s/NNN/${n}/g" < rocs-template.parm | sed "s/MOLECULE/${file}/g" > rocs-molecule${n}.parm
sed "s/NNN/${n}/g" < eon-rocs-template.parm > eon-rocs-molecule${n}.parm
rocs -param rocs-molecule${n}.parm > rocs-molecule${n}.pdb
eon -param eon-rocs-molecule${n}.parm 
let n=($n + 1)
done
