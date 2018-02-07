#!/bin/bash
#Bash script to run the code for my bachelor project

rm -r ./TTVFast/c_version/input

cd code
python inputTTV.py
cd ..

mv ./code/input ./TTVFast/c_version

cd TTVFast/c_version

for file1 in input/*.in
do
	echo -e "$file1\n-1045\n0.54\n1700\n3\n0" > setup_file.txt
	./run_TTVFast setup_file.txt Times RV_file RV_out
	mv Times output
	mv ./output/Times output/$file1
	echo $file1
done

cd ../..
mv TTVFast/c_version/output/input/*.in code/times/input/

cd code
for file2 in times/input/*.in
do
	echo $file2
	python transits.py $file2
done


