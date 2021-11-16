#!/usr/bin/bash

export OMP_NUM_THREADS=4

make clean
make

cp ./simple_fluid_simulation ./temp/.
cp ./config/config.txt ./temp/.

cd temp 
rm -rf cfg* read* energy*

./simple_fluid_simulation 1. 2. 3.