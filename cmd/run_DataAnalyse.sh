#! /usr/bin/bash

AnalyseCode="v_corr.py"
AnalyseOut="out_v_corr.txt"
AnalyseData="cfg_*.xyz"
DataFile="$HOME/FS_Data_2/"
# DataFile="$HOME/TempData/"
OutputDir="analyse/"
OutputFile="v_corr*"

cd ..
if [ -z $(echo $PYTHONPATH | grep $(pwd))]; then
	echo $(pwd)
	export PYTHONPATH=$PYTHONPATH:$(pwd)
fi

for fn in $(find $DataFile -type d -name "s_*"); do
	cd ${fn}
	pwd
	wd=$(pwd)
	wd=${wd##*/s_}
	s=${wd%%_phi*}
	echo $s

	cp $HOME/fluid_simulation/cmd/script/$AnalyseCode ${fn}/$AnalyseCode
	rm -rf ./$AnalyseData $AnalyseOut
	nohup python ${fn}/$AnalyseCode >>$AnalyseOut &&
		cd $OutputDir &&
		for Ofn in $OutputFile; do cp $Ofn $HOME/Data/TempData/s_${s}_$Ofn; done &&
		cd .. &&
		# cp ./$OutputFile $HOME/Data/TempData/. &&
		mv ./$AnalyseFile/$AnalyseData $HOME/Data/TempData/. &

	while :; do
		ps_num=$(($(ps | grep python | wc -l)))
		echo ${ps_num}
		if ((${ps_num} < 20)); then
			sleep 2
			break
		else
			echo "sleeping"
			sleep 30m
		fi
	done
done

while :; do
	ps_num=$(($(ps | grep python | wc -l)))
	echo ${ps_num}
	if ((${ps_num} < 1)); then
		sleep 2
		break
	else
		echo "sleeping"
		sleep 30m
	fi
done
