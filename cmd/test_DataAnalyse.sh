#!/usr/bin/bash

DataDir="${HOME}/Data/s*"

if [ -z `echo $PYTHONPATH | grep $(pwd)/fluid_analyse`]; then
	export PYTHONPATH=$PYTHONPATH:`pwd`/fluid_analyse
fi

cd $DataDir
pwd

# cd temp

# python all_in_one.py
