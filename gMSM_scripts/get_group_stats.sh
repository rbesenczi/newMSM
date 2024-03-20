#!/bin/bash

dataset="HCP"
workdir="$HOME/${dataset}_to_template"
reg_folder=$workdir/output/
metrics="$workdir/results"
mesh=$workdir/sunet.ico-6.sphere.surf.gii
group_list=$HOME/groupwise/${dataset}/group_list.txt
clustering=$HOME/groupwise/data/frontal_subject_clusters_${dataset}.csv

mkdir $metrics

#groups=(NODE1860 NODE1885 NODE1997) #for testing
groups=( $(cat $group_list | cut -d ',' -f1) ) # for all groups

mkdir $metrics/distortion_files/

for group in "${groups[@]}"
do
    echo "Calculating metrics of group $group."

    mkdir $metrics/distortion_files/$group

    subjects=()
    while IFS="," read -r line subject group_id
    do
        if [ $group_id = $group ]; then
          subjects+=($subject)
        fi
    done < $clustering

    datamerge="wb_command -metric-merge $metrics/$group.merge.sulc.affine.dedrifted.ico6.shape.gii "
    arealmerge="wb_command -metric-merge $metrics/$group.areal.distortion.merge.sulc.affine.dedrifted.ico6.shape.gii "
    shapemerge="wb_command -metric-merge $metrics/$group.shape.distortion.merge.sulc.affine.dedrifted.ico6.shape.gii "

    for subject in "${subjects[@]}"
    do
        datamerge+="-metric $reg_folder/$subject.MSMSulc.ico6.transformed_and_reprojected.func.gii "
        wb_command -surface-distortion $mesh $reg_folder/$subject.MSMSulc.ico6.sphere.reg.surf.gii $metrics/distortion_files/$group/$subject.MSMSulc.ico6.sphere.distortion.func.gii -local-affine-method -log2
        arealmerge+="-metric $metrics/distortion_files/$group/$subject.MSMSulc.ico6.sphere.distortion.func.gii -column 1 "
        shapemerge+="-metric $metrics/distortion_files/$group/$subject.MSMSulc.ico6.sphere.distortion.func.gii -column 2 "
    done

    $datamerge
    $arealmerge
    $shapemerge

    wb_command -metric-reduce $metrics/$group.merge.sulc.affine.dedrifted.ico6.shape.gii MEAN $metrics/$group.mean.sulc.affine.dedrifted.ico6.shape.gii
    wb_command -metric-reduce $metrics/$group.merge.sulc.affine.dedrifted.ico6.shape.gii STDEV $metrics/$group.stdev.sulc.affine.dedrifted.ico6.shape.gii

    wb_command -set-structure $metrics/$group.merge.sulc.affine.dedrifted.ico6.shape.gii CORTEX_LEFT
    wb_command -set-structure $metrics/$group.mean.sulc.affine.dedrifted.ico6.shape.gii CORTEX_LEFT
    wb_command -set-structure $metrics/$group.stdev.sulc.affine.dedrifted.ico6.shape.gii CORTEX_LEFT
done
