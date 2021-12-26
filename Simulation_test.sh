#!/usr/bin/bash

export OMP_NUM_THREADS=4
MD_time=1
time_step=1e-4
density=0.1
gamma=0 # 0 MD >0 BD Simulation

cd build

cmake ..

make clean
make

pwd
cp ./main ../temp/.
cp ../config/config.txt ../temp/.

cd ../temp
mv main fluid_simulation
rm -rf cfg* read* energy*

./fluid_simulation $MD_time $time_step $density $gamma