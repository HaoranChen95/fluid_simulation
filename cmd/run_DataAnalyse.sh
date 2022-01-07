#! /usr/bin/bash

AnalyseCode="v_corr.py"
AnalyseOut="out_v_corr.txt"
AnalyseData="cfg_*.xyz"
DataFile="$HOME/FS_Data_3/"
# DataFile="$HOME/TempData/"
OutputFile="analyse/v_corr*"

cd ..
if [ -z $(echo $PYTHONPATH | grep $(pwd))]; then
	echo $(pwd)
	export PYTHONPATH=$PYTHONPATH:$(pwd)
fi

for fn in $(find $DataFile -type d -name "s_*"); do
	cd ${fn}
	pwd

	cp $HOME/fluid_simulation/cmd/script/$AnalyseCode ${fn}/$AnalyseCode
	rm -rf ./$AnalyseData $AnalyseOut
	nohup python ${fn}/$AnalyseCode >>$AnalyseOut &&
		cp ./$OutputFile $HOME/Data/TempData/. &&
		mv ./$AnalyseFile/cfg*.xyz $HOME/Data/TempData/. &

	while :; do
		ps_num=$(($(ps | grep python | wc -l)))
		echo ${ps_num}
		if ((${ps_num} < 30)); then
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
