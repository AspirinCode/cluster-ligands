#!/bin/bash
name1=$1
name2=$2
echo $name1
echo $name2
for file in `ls mod*rpt`
do 
awk -v n=$name1 '{if ($2==n) {print $0}}' < $file > found-$file
done
rm `ls -l found* | awk '{if ($5==0) {print $9}}'`
awk -v n=$name2 '{if ($1==n) {print 2-($4+$7)}}' < found-*

