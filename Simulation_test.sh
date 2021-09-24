#!/usr/bin/bash

Pe=500
dt=1.e-3
MDt=5000
phi=0.1
Nm1p=10
Kb=0.0
ic=0
sname=0

export OMP_NUM_THREADS=4

make clean
make

cp ./sf_simu ./temp/.
cp ./config/config.txt ./temp/.

cd temp 

./sf_simu