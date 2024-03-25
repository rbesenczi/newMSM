#!/bin/bash -l

#SBATCH --partition=cpu
#SBATCH --time=0-48:00
#SBATCH --nodes=1
#SBATCH --mem=1024
#SBATCH --ntasks=1
#SBATCH --job-name=cgMSM

dataset=HCP
workdir=/scratch/prj/cortical_imaging/Renato/groupwise/${dataset}
blocks_folder=$workdir/blocks
blocks_len=$(($(ls -l $blocks_folder | wc -l)-1))

for (( filenum=1; filenum<=blocks_len; filenum++ ))
do  
	filelen=$(cat $blocks_folder/block_$filenum.txt | wc -l)
    echo "sbatch --wait --array=1-${filelen} --output=./logs/groupwise_bunch_${filenum}_%a.txt --export=ALL,FILENUM=${filenum} cross_register.sh"
	time sbatch --wait --array=1-${filelen} --output=./logs/groupwise_bunch_${filenum}_%a.txt --export=ALL,FILENUM=${filenum} cross_register.sh
	
	if [ $(head -c 1 $blocks_folder/block_$filenum.txt) -eq 1 ]
	then
		cat $workdir/new_clusters/frontal_subject_clusters_HCP_* >> $workdir/frontal_subject_clusters_HCP.csv
		rm $workdir/new_clusters/frontal_subject_clusters_HCP_*
	fi
done
