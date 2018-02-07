#!/bin/bash
#Bash script to run the code for my bachelor project

python inputTTV.py

mv /code/input TTVFast/c_version/input

cd /TTVFast/c_version

for file in input/*.in
do
	echo "test" > $file
done
#~ cd ../TTVFast/c_version
#~ ./run_TTVFast setup_file Times RV_file RV_out
