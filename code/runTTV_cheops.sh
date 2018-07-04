#!/bin/bash
#Bash script to run TTVFast on the system clones in the CHEOPS range

#~ cd code
rm transAmplCheops.csv
rm ampErrorCheops.csv
rm AmplPeriod.csv

python CHEOPSClones.py
cd ..

cp ./code/CHEOPS_clones/*.in ./TTVFast/c_version/CHEOPS/
cp ./code/CHEOPS_numP.csv ./TTVFast/c_version
cp ./code/wtm-CHEOPS_RA_Dec.csv ./TTVFast/c_version


cd TTVFast/c_version

readarray -t numPlanet < <(cut -d, -f2 CHEOPS_numP.csv)
awk -F "\"*,\"*" '{$3=$3*27; if(NR>1)print $3}' wtm-CHEOPS_RA_Dec.csv > maxTime.csv
awk -F "\"*,\"*" '{$3=$3*27;if(NR>1)print $4}' wtm-CHEOPS_RA_Dec.csv > minTime.csv
awk -F "\"*,\"*" '{$3=$3*27;if(NR>1)print $5}' wtm-CHEOPS_RA_Dec.csv > avgTime.csv
readarray -t maxTime < <(cut -d, -f2 maxTime.csv)
readarray -t minTime < <(cut -d, -f2 minTime.csv)
readarray -t avgTime < <(cut -d, -f2 avgTime.csv)


number=0
counter=1
for file1 in `ls CHEOPS/*.in | sort --version-sort`;
do
	readarray -t data < <(cut -f2 $file1)
	period=$(echo ${data[3]} | awk '{$1=$1/20; print $1;}')
	echo -e "$file1\n0\n$period\n365\n${numPlanet[$number]}\n0" > setup_file.txt;  		#${maxTime[$number]}
	./run_TTVFast setup_file.txt Times RV_file RV_out;
	mv Times CHEOPS_out;
	mv ./CHEOPS_out/Times CHEOPS_out/$file1;
	if (( $counter % 100 == 0 )) 
	then
		number=$((number+1));
		echo $number;
	fi
	counter=$((counter+1));

	#~ read -p "Press key to continue.. " -n1 -s

done


cd ../..
mv TTVFast/c_version/CHEOPS_out/CHEOPS/*.in code/CHEOPS_in/

cd code
count=0
for file2 in `ls CHEOPS_in/*.in | sort --version-sort`;
do
	echo $file2
	python transits_cheops.py $file2 $count
	count=$((count+1))
	#~ read -p "Press key to continue.. " -n1 -s
done
