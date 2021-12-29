#!/usr/bin/bash

#SBATCH --out=out_%j.txt
#SBATCH --cpus-per-task=2
#SBATCH --time=1-0

module purge
module load DEVELOP cmake/3.21.1 intel/19.1
module list
whereis cmake
which cmake
cmake --version


export OMP_NUM_THREADS=2
MD_time=1
time_step=1e-4
density=0.4
kT=1.0
gamma=0 # 0 MD >0 BD Simulation

cd $HOME/fluid_simulation

mkdir build
cd build
cmake ..
make

pwd
cp ./main ../temp/.
cp ../config/config.txt ../temp/.

cd ../temp
mv main fluid_simulation
rm -rf cfg* read* energy*

./fluid_simulation $MD_time $time_step $density $kT $gamma
