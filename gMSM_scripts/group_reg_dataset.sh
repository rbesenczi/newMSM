#!/bin/bash -l

dataset=HCP
workdir=$HOME/groupwise/${dataset}
grouplist=$workdir/group_list.txt

while IFS=',' read -r group_id size
do
	ntasks=8
	hours=12
	mem=24G

	if [ $size -gt 15 ]
	then
		ntasks=16
		hours=24
		mem=48G
	fi
	if [ $size -gt 31 ]
	then
		ntasks=32
		hours=48
		mem=96G
	fi
	if [ $size -gt 63 ]
	then
		ntasks=64
		hours=48
		mem=192G
	fi
	if [ $size -gt 127 ]
	then
		ntasks=64
		hours=48
		mem=384G
	fi

	echo "sbatch --partition=cpu --nodes=1 --job-name=${group_id} --output=${workdir}/logs/${group_id}.txt --ntasks=${ntasks} --time=0-${hours}:00 --mem=${mem} --export=ALL,ntasks=${ntasks},group_id=${group_id} run_gMSM.sh"
	sbatch --partition=cpu --nodes=1 --job-name=${group_id} --output=${workdir}/logs/${group_id}.txt --ntasks=${ntasks} --time=0-${hours}:00 --mem=${mem} --export=ALL,ntasks=${ntasks},group_id=${group_id} run_gMSM.sh

done < $grouplist
