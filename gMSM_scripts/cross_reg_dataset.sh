#!/bin/bash -l

#SBATCH --partition=cpu
#SBATCH --time=0-48:00
#SBATCH --nodes=1
#SBATCH --mem=1024
#SBATCH --ntasks=1
#SBATCH --job-name=cgMSMctr

dataset=HCP
workdir=/scratch/users/k2258483/groupwise/${dataset}
blocks_folder=$workdir/blocks
blocks_len=$(($(ls -l $blocks_folder | wc -l)-1))

for (( filenum=1; filenum<=blocks_len; filenum++ ))
do  
	filelen=$(cat $blocks_folder/block_$filenum.txt | wc -l)
	
	if [ $(head -c 1 $blocks_folder/block_$filenum.txt) -eq 0 ]
	then
	    echo "sbatch --wait --time=0-1:00 --ntasks=2 --array=1-${filelen} --output=${workdir}/logs/groupwise_bunch_${filenum}_%a.txt --export=ALL,FILENUM=${filenum} cross_register.sh"
		time sbatch --wait --time=0-1:00 --ntasks=2 --array=1-${filelen} --output=${workdir}/logs/groupwise_bunch_${filenum}_%a.txt --export=ALL,FILENUM=${filenum} cross_register.sh
	fi

    retval=$?
    if [ $retval -ne 0 ]; then
        echo "Error at step $filenum"
        exit $retval
    fi

	if [ $(head -c 1 $blocks_folder/block_$filenum.txt) -eq 1 ]
	then
	    echo "sbatch --wait --time=0-1:00 --ntasks=1 --array=1-${filelen} --output=${workdir}/logs/groupwise_bunch_${filenum}_%a.txt --export=ALL,FILENUM=${filenum} cross_register.sh"
		time sbatch --wait --time=0-1:00 --ntasks=1 --array=1-${filelen} --output=${workdir}/logs/groupwise_bunch_${filenum}_%a.txt --export=ALL,FILENUM=${filenum} cross_register.sh
		cat $workdir/new_clusters/frontal_subject_clusters_HCP_* >> $workdir/frontal_subject_clusters_HCP.csv
		rm $workdir/new_clusters/frontal_subject_clusters_HCP_*
	fi

    retval=$?
    if [ $retval -ne 0 ]; then
        echo "Error at step $filenum"
        exit $retval
    fi
done
