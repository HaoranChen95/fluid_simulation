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

source compiler-select intel-fi

export OMP_NUM_THREADS=2
MD_time=500
time_step=1e-3
density=0.4
kT=1
gamma=1 # 0 MD >0 BD Simulation

cd ..
pwd

mkdir build
cd build
cmake ..
# cmake -DCMAKE_C_COMPILER=/home/opt/intel/oneapi/compiler/2021.1.1/linux/bin/intel64/icc -DCMAKE_CXX_COMPILER=/home/opt/intel/oneapi/compiler/2021.1.1/linux/bin/intel64/icpc ..
make

pwd
cp ./main ../temp/.
cp ../config/config.txt ../temp/.

cd ../temp
mv main fluid_simulation
rm -rf cfg* read* energy*

./fluid_simulation $MD_time $time_step $density $kT $gamma
