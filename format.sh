#!/bin/bash
rm -f format*list
for dir in `ls -d state*/`
do
state=${dir%%/}
for file in `cat ranked-${state}-dbase.list | grep -v "ICI" | grep -v CL | grep -v Zil | grep -v rac_2 | grep -v rac_3`
do
dir=${file%%molec*}
tmp=${file##*molecule_}
name=${tmp%%_000.mol2}
newfile=format-molecule_${name}_000.mol2
sed "s/title\: molecule 1/$name/g" < $file > $dir/$newfile
echo /home/mlawrenz/wontkill/new-results/bi-bindb-pgeomx/${state}/$newfile >> format-${state}-dbase.list
done
for file in  `cat ranked-${state}-dbase.list | grep "_ICI118551"`
do
dir=${file%%molec*}
tmp=${file##*molecule__}
name=${tmp%%_000.mol2}
newfile=format-molecule_${name}_000.mol2
sed "s/title\: molecule 1/ICI118551/g" < $file  > $dir/$newfile
echo /home/mlawrenz/wontkill/new-results/bi-bindb-pgeomx/${state}/$newfile >> format-${state}-dbase.list
done
for file in  `cat ranked-${state}-dbase.list | grep "_ICI_"`
do
dir=${file%%molec*}
tmp=${file##*molecule_}
name=${tmp%%_000.mol2}
newfile=format-molecule_ICI89406_000.mol2
sed "s/title\: molecule 1/ICI89406/g" < $file  > $dir/$newfile
echo /home/mlawrenz/wontkill/new-results/bi-bindb-pgeomx/${state}/$newfile >> format-${state}-dbase.list
done
for file in  `cat ranked-${state}-dbase.list | grep CL`
do
dir=${file%%molec*}
tmp=${file##*molecule_}
name=${tmp%%_000.mol2}
newfile=format-molecule_${name}_000.mol2
sed "s/title\: molecule 1/CL316243/g" < $file  > $dir/$newfile
echo /home/mlawrenz/wontkill/new-results/bi-bindb-pgeomx/${state}/$newfile >> format-${state}-dbase.list
done
for file in  `cat ranked-${state}-dbase.list | grep rac_Zilpaterol`
do
dir=${file%%molec*}
tmp=${file##*molecule__}
name=${tmp%%_000.mol2}
newfile=format-molecule_rac-Zilpaterol_000.mol2
sed "s/title\: molecule 1/rac-Zilpaterol/g" < $file  > $dir/$newfile
echo /home/mlawrenz/wontkill/new-results/bi-bindb-pgeomx/${state}/$newfile >> format-${state}-dbase.list
done
for file in  `cat ranked-${state}-dbase.list | grep pos_Zilpaterol`
do
dir=${file%%molec*}
tmp=${file##*molecule__}
name=${tmp%%_000.mol2}
newfile=format-molecule_pos-Zilpaterol_000.mol2
sed "s/title\: molecule 1/pos-Zilpaterol/g" < $file  > $dir/$newfile
echo /home/mlawrenz/wontkill/new-results/bi-bindb-pgeomx/${state}/$newfile >> format-${state}-dbase.list
done
for file in  `cat ranked-${state}-dbase.list | grep neg_Zilpaterol`
do
dir=${file%%molec*}
tmp=${file##*molecule__}
name=${tmp%%_000.mol2}
newfile=format-molecule_neg-Zilpaterol_000.mol2
sed "s/title\: molecule 1/neg-Zilpaterol/g" < $file  > $dir/$newfile
echo /home/mlawrenz/wontkill/new-results/bi-bindb-pgeomx/${state}/$newfile >> format-${state}-dbase.list
done
for file in  `cat ranked-${state}-dbase.list | grep rac_2`
do
dir=${file%%molec*}
tmp=${file##*molecule__}
name=${tmp%%_000.mol2}
newfile=format-molecule_Clenbutarol_000.mol2
sed "s/title\: molecule 1/Clenbutarol/g" < $file  > $dir/$newfile
echo /home/mlawrenz/wontkill/new-results/bi-bindb-pgeomx/${state}/$newfile >> format-${state}-dbase.list
done
for file in  `cat ranked-${state}-dbase.list | grep rac_3`
do
dir=${file%%molec*}
tmp=${file##*molecule__}
name=${tmp%%_000.mol2}
newfile=format-molecule_rac-3_000.mol2
sed "s/title\: molecule 1/rac-3/g" < $file  > $dir/$newfile
echo /home/mlawrenz/wontkill/new-results/bi-bindb-pgeomx/${state}/$newfile >> format-${state}-dbase.list
done
echo "done with $state"
done
