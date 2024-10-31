#!/bin/bash -l

###########################################################
##################### SET BY THE USER #####################
input_folder=$HOME/affined_HCP_sulc_curv
workdir=$HOME/HCP_to_template

clustering=$workdir/frontal_subject_clusters_HCP.csv
template=$workdir/templates/MSMStrain.L.sulc.curv.ico6.shape.gii
tmp_mesh=$workdir/templates/sunet.ico-6.sphere.surf.gii

config_file=$workdir/configs/HCP_to_template_config.txt
outdir=$workdir/output
resultdir=$workdir/results

group_id=NODE2078
###########################################################

mkdir $outdir
mkdir $resultdir

datasulcmerge="wb_command -metric-merge $resultdir/$group_id.merge.sulc.affine.dedrifted.ico6.shape.gii "
datacurvmerge="wb_command -metric-merge $resultdir/$group_id.merge.curv.affine.dedrifted.ico6.shape.gii "
arealmerge="wb_command -metric-merge $resultdir/$group_id.areal.distortion.merge.curv.affine.dedrifted.ico6.shape.gii "
shapemerge="wb_command -metric-merge $resultdir/$group_id.shape.distortion.merge.curv.affine.dedrifted.ico6.shape.gii "
numsub=0

while IFS=',' read -r linenum subject group
do
  if [ $group_id = $group ]; then

    echo "Registering subject $subject to the template..."

    time $HOME/fsldev/bin/newmsm \
      --indata=$input_folder/$subject.sulc.curv.affine.ico6.shape.gii \
      --inmesh=$input_folder/$subject.sunet.ico-6.surf.gii \
      --refdata=$template \
      --refmesh=$tmp_mesh \
      --conf=$config_file \
      --out=$outdir/$subject.MSMSulc.ico6. \
      --verbose

    retval=$?
    if [ $retval -ne 0 ]; then
        echo "Error at subject $subject with code $retval."
        exit $retval
    fi

    wb_command -surface-distortion $input_folder/$subject.sunet.ico-6.surf.gii $outdir/$subject.MSMSulc.ico6.sphere.reg.surf.gii $outdir/$subject.MSMSulc.ico6.sphere.distortion.func.gii -local-affine-method -log2

    datasulcmerge+="-metric $outdir/$subject.MSMSulc.ico6.transformed_and_reprojected.func.gii -column 1 "
    datacurvmerge+="-metric $outdir/$subject.MSMSulc.ico6.transformed_and_reprojected.func.gii -column 2 "
    arealmerge+="-metric $outdir/$subject.MSMSulc.ico6.sphere.distortion.func.gii -column 1 "
    shapemerge+="-metric $outdir/$subject.MSMSulc.ico6.sphere.distortion.func.gii -column 2 "
    ((numsub++))
  fi

done < $clustering

echo "Calculating metrics of group $group_id with $numsub subjects."

$datasulcmerge
$datacurvmerge
$arealmerge
$shapemerge

wb_command -metric-reduce $resultdir/$group_id.merge.sulc.affine.dedrifted.ico6.shape.gii MEAN $resultdir/$group_id.mean.sulc.affine.dedrifted.ico6.shape.gii
wb_command -metric-reduce $resultdir/$group_id.merge.sulc.affine.dedrifted.ico6.shape.gii STDEV $resultdir/$group_id.stdev.sulc.affine.dedrifted.ico6.shape.gii
wb_command -metric-reduce $resultdir/$group_id.merge.curv.affine.dedrifted.ico6.shape.gii MEAN $resultdir/$group_id.mean.curv.affine.dedrifted.ico6.shape.gii
wb_command -metric-reduce $resultdir/$group_id.merge.curv.affine.dedrifted.ico6.shape.gii STDEV $resultdir/$group_id.stdev.curv.affine.dedrifted.ico6.shape.gii

wb_command -set-structure $resultdir/$group_id.merge.sulc.affine.dedrifted.ico6.shape.gii CORTEX_LEFT
wb_command -set-structure $resultdir/$group_id.mean.sulc.affine.dedrifted.ico6.shape.gii CORTEX_LEFT
wb_command -set-structure $resultdir/$group_id.stdev.sulc.affine.dedrifted.ico6.shape.gii CORTEX_LEFT
wb_command -set-structure $resultdir/$group_id.merge.curv.affine.dedrifted.ico6.shape.gii CORTEX_LEFT
wb_command -set-structure $resultdir/$group_id.mean.curv.affine.dedrifted.ico6.shape.gii CORTEX_LEFT
wb_command -set-structure $resultdir/$group_id.stdev.curv.affine.dedrifted.ico6.shape.gii CORTEX_LEFT
