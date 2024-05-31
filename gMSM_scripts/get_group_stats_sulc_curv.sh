#!/bin/bash

dataset="HCP"
workdir="$HOME/${dataset}Sulc_Curv_to_template"
reg_folder=$workdir/output/
metrics="$workdir/results"
mesh=$workdir/sunet.ico-6.sphere.surf.gii
group_list=$HOME/groupwise/${dataset}/group_list.txt
clustering=$HOME/groupwise/${dataset}/frontal_subject_clusters_${dataset}.csv

mkdir $metrics

#groups=( $(cat $group_list | cut -d ',' -f1) ) #for all groups
groups=(NODE2159) #for testing

for group in "${groups[@]}"
do
    subjects=()
    while IFS="," read -r line subject group_id
    do
        if [ $group_id = $group ]; then
          subjects+=($subject)
        fi
    done < $clustering

    datasulcmerge="wb_command -metric-merge $metrics/$group.merge.sulc.affine.dedrifted.ico6.shape.gii "
    datacurvmerge="wb_command -metric-merge $metrics/$group.merge.curv.affine.dedrifted.ico6.shape.gii "
    arealmerge="wb_command -metric-merge $metrics/$group.areal.distortion.merge.curv.affine.dedrifted.ico6.shape.gii "
    shapemerge="wb_command -metric-merge $metrics/$group.shape.distortion.merge.curv.affine.dedrifted.ico6.shape.gii "

    for subject in "${subjects[@]}"
    do
        datasulcmerge+="-metric $reg_folder/$subject.MSMSulc.Curv.ico6.transformed_and_reprojected.func.gii -column 1 "
        datacurvmerge+="-metric $reg_folder/$subject.MSMSulc.Curv.ico6.transformed_and_reprojected.func.gii -column 2 "
        arealmerge+="-metric $reg_folder/$subject.MSMSulc.Curv.ico6.sphere.distortion.func.gii -column 1 "
        shapemerge+="-metric $reg_folder/$subject.MSMSulc.Curv.ico6.sphere.distortion.func.gii -column 2 "
    done

    echo "Calculating metrics of group $group with ${#subjects[@]} subjects."

    $datasulcmerge
    $datacurvmerge
    $arealmerge
    $shapemerge

    wb_command -metric-reduce $metrics/$group.merge.sulc.affine.dedrifted.ico6.shape.gii MEAN $metrics/$group.mean.sulc.affine.dedrifted.ico6.shape.gii
    wb_command -metric-reduce $metrics/$group.merge.sulc.affine.dedrifted.ico6.shape.gii STDEV $metrics/$group.stdev.sulc.affine.dedrifted.ico6.shape.gii
    wb_command -metric-reduce $metrics/$group.merge.curv.affine.dedrifted.ico6.shape.gii MEAN $metrics/$group.mean.curv.affine.dedrifted.ico6.shape.gii
    wb_command -metric-reduce $metrics/$group.merge.curv.affine.dedrifted.ico6.shape.gii STDEV $metrics/$group.stdev.curv.affine.dedrifted.ico6.shape.gii

    wb_command -set-structure $metrics/$group.merge.sulc.affine.dedrifted.ico6.shape.gii CORTEX_LEFT
    wb_command -set-structure $metrics/$group.mean.sulc.affine.dedrifted.ico6.shape.gii CORTEX_LEFT
    wb_command -set-structure $metrics/$group.stdev.sulc.affine.dedrifted.ico6.shape.gii CORTEX_LEFT
    wb_command -set-structure $metrics/$group.merge.curv.affine.dedrifted.ico6.shape.gii CORTEX_LEFT
    wb_command -set-structure $metrics/$group.mean.curv.affine.dedrifted.ico6.shape.gii CORTEX_LEFT
    wb_command -set-structure $metrics/$group.stdev.curv.affine.dedrifted.ico6.shape.gii CORTEX_LEFT
done
