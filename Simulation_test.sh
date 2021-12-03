#!/usr/bin/bash

export OMP_NUM_THREADS=4
MD_time=100
time_step=1e-4
density=0.1
gamma=10 # 0 MD >0 BD Simulation

make clean
make

cp ./simple_fluid_simulation ./temp/.
cp ./config/config.txt ./temp/.

cd temp 
rm -rf cfg* read* energy*

./simple_fluid_simulation $MD_time $time_step $density $gamma