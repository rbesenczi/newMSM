#!/bin/bash

dataset=HCP
workdir=$HOME/groupwise/${dataset}
grouplist=$workdir/group_list.txt

while IFS=',' read -r group_id size
do
	ntasks=4
	hours=24
	mem=32G

	if [ $size -gt 15 ]
	then
		ntasks=8
		hours=48
		mem=64G
	fi
	if [ $size -gt 31 ]
	then
		ntasks=16
		hours=48
		mem=256G
	fi
	if [ $size -gt 63 ]
	then
		ntasks=16
		hours=48
		mem=256G
	fi

	echo "sbatch --partition=cpu --nodes=1 --job-name=${group_id} --output=${workdir}/logs/${group_id}.txt --ntasks=${ntasks} --time=0-${hours}:00 --mem=${mem} --export=ALL,ntasks=${ntasks},group_id=${group_id} run_gMSM.sh"
	sbatch --partition=cpu --nodes=1 --job-name=${group_id} --output=${workdir}/logs/${group_id}.txt --ntasks=${ntasks} --time=0-${hours}:00 --mem=${mem} --export=ALL,ntasks=${ntasks},group_id=${group_id} run_gMSM.sh

done < $grouplist
