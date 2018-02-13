#!/bin/bash
#Bash script to run TTVFast

rm -r ./TTVFast/c_version/input

cd code
python inputTTV.py
cd ..

mv ./code/input ./TTVFast/c_version
mv ./code/numberPlanets.csv

cd TTVFast/c_version

readarray -t numPlanet < <(cut -d, -f2 numberPlanets.csv)
echo ${numPlanet[@]}
$number="0"
for file1 in input/*.in
do
	echo -e "$file1\n0\n0.54\n1460\n${numPlanet[$number]}\n0" > setup_file.txt
	./run_TTVFast setup_file.txt Times RV_file RV_out
	mv Times output
	mv ./output/Times output/$file1
	#~ echo $file1
	number=$((number+1))
done


cd ../..
mv TTVFast/c_version/output/input/*.in code/times/input/

cd code
for file2 in times/input/*.in
do
	#~ echo $file2
	python transits.py $file2
done


