#!/bin/bash

for x in `cat reference.list`
do 
test=`grep -m 1 ${x} ../apo-2665/eon*bind*rpt | head -1`
if [ -z "$test" ]
then 
echo $x
fi
done
