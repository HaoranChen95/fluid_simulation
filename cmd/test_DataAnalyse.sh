#!/usr/bin/bash

# DataDir="${HOME}/Data/s*"
DataDir="${HOME}/Archive/s_diff_density_phi_0.05_kT_1.0_gamma_0/"

cd ..

if [ -z `echo $PYTHONPATH | grep $(pwd)`]; then
	echo `pwd`
	export PYTHONPATH=$PYTHONPATH:`pwd`
fi

echo $PYTHONPATH
cp cmd/test.py $DataDir.

cd $DataDir
pwd

python test.py
