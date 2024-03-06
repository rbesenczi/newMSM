#!/bin/bash -l

groups=(NODE1888 NODE1915 NODE1994 NODE1998)

for group in "${groups[@]}"
do
    input_folder=/scratch/prj/cortical_imaging/Renato/HCP_to_template
    subject_list=$HOME/newMSM_HCP_to_tmp/group_lists/${group}_HPC.csv
    template=/users/k2258483/newMSM_HCP_to_tmp/template/sunet.ico-6.template.surf.gii

    subjects=( $(cat $subject_list | cut -d ',' -f1) )

    datamerge="wb_command -metric-merge $input_folder/metric_merge/$group.merge.sulc.affine.dedrifted.ico6.shape.gii "
    arealmerge="wb_command -metric-merge $input_folder/distortions/$group.areal.distortion.merge.sulc.affine.dedrifted.ico6.shape.gii "
    shapemerge="wb_command -metric-merge $input_folder/distortions/$group.shape.distortion.merge.sulc.affine.dedrifted.ico6.shape.gii "

    for subject in "${subjects[@]}"
    do
        echo "Calculating distortion and metric_merge for subject $subject in group $group."
        datamerge+="-metric $input_folder/1/$subject.MSMSulc.ico6.transformed_and_reprojected.func.gii "
        wb_command -surface-distortion $template $input_folder/1/$subject.MSMSulc.ico6.sphere.reg.surf.gii $input_folder/1/$subject.MSMSulc.ico6.sphere.distortion.func.gii -local-affine-method -log2
        arealmerge+="-metric $input_folder/1/$subject.MSMSulc.ico6.sphere.distortion.func.gii -column 1 "
        shapemerge+="-metric $input_folder/1/$subject.MSMSulc.ico6.sphere.distortion.func.gii -column 2 "
    done

    $datamerge
    $arealmerge
    $shapemerge

    wb_command -metric-reduce $input_folder/metric_merge/$group.merge.sulc.affine.dedrifted.ico6.shape.gii MEAN $input_folder/metric_merge/$group.mean.sulc.affine.dedrifted.ico6.shape.gii
    wb_command -metric-reduce $input_folder/metric_merge/$group.merge.sulc.affine.dedrifted.ico6.shape.gii STDEV $input_folder/metric_merge/$group.stdev.sulc.affine.dedrifted.ico6.shape.gii

    wb_command -set-structure $input_folder/metric_merge/$group.merge.sulc.affine.dedrifted.ico6.shape.gii CORTEX_LEFT
    wb_command -set-structure $input_folder/metric_merge/$group.mean.sulc.affine.dedrifted.ico6.shape.gii CORTEX_LEFT
    wb_command -set-structure $input_folder/metric_merge/$group.stdev.sulc.affine.dedrifted.ico6.shape.gii CORTEX_LEFT
done

