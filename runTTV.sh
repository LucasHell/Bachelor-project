#!/bin/bash
#Bash script to run the code for my bachelor project

cd code
python inputTTV.py
cd ..

mv ./code/input ./TTVFast/c_version

cd TTVFast/c_version

for file in input/*.in
do
	echo -e "$file \n-1045\n0.54\n1700\n2\0" > setup_file 
done

rm -r input

# ./run_TTVFast setup_file Times RV_file RV_out
