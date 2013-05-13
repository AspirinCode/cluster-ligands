awk '{ if ($1~"-" || $1 =="Terbutaline") {print $0}}'

for file in `ls eon*rpt`
do
sed "s/Terbutaline/Terbutaline-1/g" < $file | sed "s/Xamoterol/Xamoterol-1/g" | sed "s/Timolol/Timolol-1/g" | awk '{ if ($1~"-") {print $0}}' | sed "s/(-)-Zilpaterol/neg-Zilpaterol/g" | sed "s/(+)-Zilpaterol/pos-Zilpaterol/g" | sed | "s/(rac)-Zilpaterol/rac-Zilpaterol/g" | sed "s/(rac)-2/rac-2/g" | sed "s/(rac)-3/rac-3/g" > mod-$file 
done
