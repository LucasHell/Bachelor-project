#!/bin/bash
#Bash script to generate files in the correct format for TTVFast, run TTVFast and plot results

#Remove input and output files from when program ran last time
rm -r ./TTVFast/c_version/input
rm ./TTVFast/c_version/numberPlanets.csv
rm ./code/timingErrors.csv
rm  -r ./code/times/*

#Run code to generate data in supported format for TTVFast
cd code
python inputTTV.py
cd ..

#Move input data to TTVFast directory
echo "TEST1"
mv ./code/input ./TTVFast/c_version
mv ./code/numberPlanets.csv ./TTVFast/c_version

echo "TEST2"

cd TTVFast/c_version

#Run TTVFast by creating a setup file for every system and rename the output file to the KOI name of the first planet in the system
readarray -t numPlanet < <(cut -d, -f2 numberPlanets.csv)
number=0
for file1 in input/*.in
do
	echo -e "$file1\n0\n0.54\n1460\n${numPlanet[$number]}\n0" > setup_file.txt
	./run_TTVFast setup_file.txt Times RV_file RV_out
	mv Times output
	mv ./output/Times output/$file1
	echo $file1
	number=$((number+1))
done

#Move output data from TTVFast to directory for plotting
cd ../..
mv TTVFast/c_version/output/input/*.in code/times/

#Run python script to read data generated by TTVFast and plot it
cd code
count=0
for file2 in times/*.in
do
	echo $file2
	python transits.py $file2 $count
	count=$((count+1))
done


