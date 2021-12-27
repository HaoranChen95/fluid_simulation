#!/usr/bin/bash

export OMP_NUM_THREADS=4
MD_time=100
time_step=1e-4
density=0.1
gamma=1 # 0 MD >0 BD Simulation

cd build

cmake ..

make clean
make

pwd
cp ./main ../temp_/.
cp ../config/config.txt ../temp_/.

cd ../temp_
mv main fluid_simulation
rm -rf cfg* read* energy*

./fluid_simulation $MD_time $time_step $density $gamma