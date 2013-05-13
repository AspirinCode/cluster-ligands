#!/bin/bash

#for dir in `ls -d state*/`
for dir in `ls -d state1018/`
do
tmp=${dir##*state}
state=${tmp%%/*}
echo "on $dir"
rm -f $dir/new*eon*rpt
rm -f $dir/new*eon*mol2
for x in `seq 1 3200`
do
file=$dir/eon-only-${state}-${x}.rpt
base=`basename $file`
reference=`awk -v n=${x} '{if (NR==n) {print $0}}' < mapped-dbase.list`
awk -v s="$state" '{print s"_"$1}' < $file  > tmp
awk -v r="$reference" '{print r, $3, $4, $5, $6, $7}' < $file > tmp2
paste tmp tmp2 >> $dir/new-eon-only-$reference.rpt
###### format mol2 files
echo "on hits file"
file=$dir/eon-only-${state}-${x}_hits.mol2
cp $file tmp0
count=0
for x in `grep -A 1 MOLE $file | grep -v MOL | grep -v "\-\-"`;
do 
let newcount=($count + 1)
if [ $count -eq 160 ]
then 
sed "s/^${x}/${state}_${x}/g" < tmp${count} > $dir/new-eon-only-$reference_hits.mol2
break
else 
sed "s/^${x}/${state}_${x}/g" < tmp${count} > tmp${newcount}
let count=($count + 1)
fi
done
done
done
