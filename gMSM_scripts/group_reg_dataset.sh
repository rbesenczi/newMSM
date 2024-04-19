#!/bin/bash -l

dataset=HCP
#workdir=/scratch/users/k2258483/groupwise/${dataset}
workdir=$HOME/groupwise/${dataset}
grouplist=$workdir/group_list.txt

while IFS=',' read -r group_id size
do
	ntasks=8
	hours=6
	mem=16384

	if [ $size -gt 16 ]
	then
		ntasks=16
		hours=12
		mem=32768
	fi
	if [ $size -gt 32 ]
	then
		ntasks=32
		hours=24
		mem=65536
	fi
	if [ $size -gt 64 ]
	then
		ntasks=64
		hours=48
		mem=131072
	fi
	if [ $size -gt 128 ]
	then
		ntasks=64
		hours=48
		mem=262144
	fi

	echo "sbatch --partition=cpu --nodes=1 --job-name=${group_id} --output=${workdir}/logs/groupwise_reg_${group_id}.txt --ntasks=${ntasks} --time=0-${hours}:00 --mem=${mem} --export=ALL,ntasks=${ntasks},group_id=${group_id} run_gMSM.sh"
	sbatch --partition=cpu --nodes=1 --job-name=${group_id} --output=${workdir}/logs/groupwise_reg_${group_id}.txt --ntasks=${ntasks} --time=0-${hours}:00 --mem=${mem} --export=ALL,ntasks=${ntasks},group_id=${group_id} run_gMSM.sh

done < $grouplist
