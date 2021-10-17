#!/usr/bin/bash

export OMP_NUM_THREADS=4

make clean
make

cp ./sf_simu ./temp/.
cp ./config/config.txt ./temp/.

cd temp 
rm -rf cfg* read*

./sf_simu