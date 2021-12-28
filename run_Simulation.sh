#!/bin/bash

bash
core=2
export OMP_NUM_THREADS=$core

CPU_NUM=$(cat /proc/cpuinfo | grep 'processor' | wc -l)
echo "The Partition have $CPU_NUM cores"

# setting the path for running code and storing data
RunDir=$HOME"/RunData/"
StoreDir=$HOME"/TempData/"
InitDir=$HOME"/InitCfg/"
# RunDir=$HOME"/Data/RunData/"
# StoreDir=$HOME"/Data/TempData/"
# InitDir=$HOME"/InitCfg/"
# set the email to send
email=ha.chen@fz-juelich.de

exe_suffix="th2_${core}_core"
# setting the system parameters of simulation

array_phi=(0.40)
dt=1e-4
MDt=1.1e3
array_kT=(1.0 2.0 3.0)
array_gamma=(0)

# setting the data series name
snum=0

# go to simulation dir and comile the code
cd $HOME/fluid_simulation/build/

source compiler-select intel-fi
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/tmp_mnt/el7/impi/include/
cmake ..
make clean
make
cp ./main ../fluid_simulation_$exe_suffix

# send the implement file to the running dir
for phi in "${array_phi[@]}"; do
	for kT in "${array_kT[@]}"; do
		for gamma in "${array_gamma[@]}"; do
			fname="s_${snum}_phi_${phi}_kT_${kT}_gamma_${gamma}"

			mkdir $RunDir$fname

			echo "start copy at"
			pwd

			cp ../fluid_simulation_$exe_suffix $RunDir$fname
			cp ../config/config.txt $RunDir$fname

			cd $RunDir$fname
			echo "[Running] the simulation in"
			pwd

			mv fluid_simulation_$exe_suffix $fname
			# srun --partition=th2 --account=chen --out=out_%j.txt \
			nohup ./$fname $MDt $dt $phi $kT $gamma >> out_simulation.txt &&
				mv $RunDir$fname ${StoreDir}/. &
			sleep 2s

			cd -

			while :; do
				ps_num=$(($(ps | grep ic_ | wc -l)))
				echo "Process Num: $ps_num"
				if (($ps_num < $CPU_NUM / ${core})); then
					sleep 2
					break
				else
					echo "sleeping"
					sleep 30m
				fi
			done
		done
	done
done

while :; do
	ps_num=$(($(ps | grep ic_ | wc -l)))
	echo "The last Process Num: $ps_num"
	if (($ps_num < 1)); then
		sleep 2
		break
	else
		echo "sleeping"
		sleep 30m
	fi
done
