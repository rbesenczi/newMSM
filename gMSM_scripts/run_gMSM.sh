#! /bin/bash -l

###########################################################
dataset=HCP

workdir=$HOME/groupwise/$dataset
input_folder=$HOME/affined_${dataset}_sulc_curv
template=$workdir/templates/sunet.ico-6.template.surf.gii
config_file=$workdir/configs/gMSM_config_8.txt
outdir=$workdir/output
resultdir=$workdir/results
clustering=$workdir/frontal_subject_clusters_${dataset}.csv

if [ $ntasks -gt 15 ]
then
  config_file=$workdir/configs/gMSM_config_16.txt
fi
if [ $ntasks -gt 31 ]
then
  config_file=$workdir/configs/gMSM_config_32.txt
fi
if [ $ntasks -gt 63 ]
then
  config_file=$workdir/configs/gMSM_config_64.txt
fi

###########################################################

mkdir $workdir/file_lists
rm $workdir/file_lists/input_data_$group_id.txt && touch $workdir/file_lists/input_data_$group_id.txt
rm $workdir/file_lists/input_meshes_$group_id.txt && touch $workdir/file_lists/input_meshes_$group_id.txt
subjects=()

while IFS="," read -r linenum subject group
do
  if [ $group_id = $group ]; then
    echo "$input_folder/$subject.sulc.curv.affine.ico6.shape.gii" >> $workdir/file_lists/input_data_$group_id.txt
    echo "$workdir/templates/sunet.ico-6.template.surf.gii" >> $workdir/file_lists/input_meshes_$group_id.txt
    subjects+=($subject)
  fi
done < $clustering

mkdir $outdir
mkdir $outdir/$group_id
mkdir $resultdir
mkdir $resultdir/$group_id

## gMSM registration phase ##
echo "Running gMSM for group $group_id"
time $HOME/fsldev/bin/newmsm \
  --data=$workdir/file_lists/input_data_$group_id.txt \
  --meshes=$workdir/file_lists/input_meshes_$group_id.txt \
  --template=$template \
  --mask=$workdir/NODE2218_frontal_mask.shape.gii \
  --conf=$config_file \
  --out=$outdir/$group_id/groupwise.$group_id. \
  --verbose --groupwise

index=0
## some renaming and setting structures ##
for subject in "${subjects[@]}"
do
  wb_command -set-structure $outdir/$group_id/groupwise.$group_id.sphere-$index.reg.surf.gii CORTEX_LEFT
  mv $outdir/$group_id/groupwise.$group_id.sphere-$index.reg.surf.gii $outdir/$group_id/groupwise.$group_id.sphere-$subject.reg.surf.gii
  wb_command -set-structure $outdir/$group_id/groupwise.$group_id.transformed_and_reprojected-$index.func.gii CORTEX_LEFT
  mv $outdir/$group_id/groupwise.$group_id.transformed_and_reprojected-$index.func.gii $outdir/$group_id/groupwise.$group_id.transformed_and_reprojected-$subject.func.gii
  ((index++))
done

## dedrifting phase ##
## calculating surface average of inverse of registration ##
surf_avg="wb_command -surface-average $outdir/$group_id/groupwise.$group_id.dedriftwarp.ico-6.sphere.surf.gii "

for subject in "${subjects[@]}"
do
  wb_command -surface-sphere-project-unproject \
  $template \
  $outdir/$group_id/groupwise.$group_id.sphere-$subject.reg.surf.gii \
  $template \
  $outdir/$group_id/groupwise.$group_id.sphere-$subject.reg.inv.surf.gii
  surf_avg+="-surf $outdir/$group_id/groupwise.$group_id.sphere-$subject.reg.inv.surf.gii "
done

$surf_avg

## recentre and normalise ##
wb_command -surface-modify-sphere $outdir/$group_id/groupwise.$group_id.dedriftwarp.ico-6.sphere.surf.gii 100 $outdir/$group_id/groupwise.$group_id.dedriftwarp.ico-6.sphere.surf.gii -recenter

for subject in "${subjects[@]}"
do
  ## applying inverse ##
  wb_command -surface-sphere-project-unproject \
  $outdir/$group_id/groupwise.$group_id.sphere-$subject.reg.surf.gii \
  $template \
  $outdir/$group_id/groupwise.$group_id.dedriftwarp.ico-6.sphere.surf.gii \
  $outdir/$group_id/groupwise.$group_id.sphere-$subject.reg.corrected.surf.gii

  ## resampling data for the new dedrifted mesh ##
  wb_command -metric-resample \
  $input_folder/$subject.sulc.curv.affine.ico6.shape.gii \
  $outdir/$group_id/groupwise.$group_id.sphere-$subject.reg.corrected.surf.gii \
  $template \
  ADAP_BARY_AREA \
  $outdir/$group_id/groupwise.$group_id.transformed_and_reprojected.dedrift-$subject.func.gii \
  -area-surfs \
  $outdir/$group_id/groupwise.$group_id.sphere-$subject.reg.corrected.surf.gii \
  $template
done

## calculating mean and stdev, areal and shape distortion ##
merge_sulc="wb_command -metric-merge $resultdir/$group_id/groupwise.$group_id.merge.sulc.affine.dedrifted.ico6.shape.gii "
merge_curv="wb_command -metric-merge $resultdir/$group_id/groupwise.$group_id.merge.curv.affine.dedrifted.ico6.shape.gii "
arealmerge="wb_command -metric-merge $resultdir/$group_id/groupwise.$group_id.areal.distortion.merge.sulc.affine.dedrifted.ico6.shape.gii "
shapemerge="wb_command -metric-merge $resultdir/$group_id/groupwise.$group_id.shape.distortion.merge.sulc.affine.dedrifted.ico6.shape.gii "

for subject in "${subjects[@]}"
do
  merge_sulc+="-metric $outdir/$group_id/groupwise.$group_id.transformed_and_reprojected.dedrift-$subject.func.gii -column 1 "
  merge_curv+="-metric $outdir/$group_id/groupwise.$group_id.transformed_and_reprojected.dedrift-$subject.func.gii -column 2 "
  wb_command -surface-distortion $template $outdir/$group_id/groupwise.$group_id.sphere-$subject.reg.corrected.surf.gii $outdir/$group_id/groupwise.$group_id.sphere-$subject.distortion.func.gii -local-affine-method -log2
  arealmerge+="-metric $outdir/$group_id/groupwise.$group_id.sphere-$subject.distortion.func.gii -column 1 "
  shapemerge+="-metric $outdir/$group_id/groupwise.$group_id.sphere-$subject.distortion.func.gii -column 2 "
done

$merge_sulc
$merge_curv
$arealmerge
$shapemerge

wb_command -metric-reduce $resultdir/$group_id/groupwise.$group_id.merge.sulc.affine.dedrifted.ico6.shape.gii MEAN $resultdir/$group_id/groupwise.$group_id.mean.sulc.affine.dedrifted.ico6.shape.gii
wb_command -metric-reduce $resultdir/$group_id/groupwise.$group_id.merge.curv.affine.dedrifted.ico6.shape.gii MEAN $resultdir/$group_id/groupwise.$group_id.mean.curv.affine.dedrifted.ico6.shape.gii
wb_command -metric-reduce $resultdir/$group_id/groupwise.$group_id.merge.sulc.affine.dedrifted.ico6.shape.gii STDEV $resultdir/$group_id/groupwise.$group_id.stdev.sulc.affine.dedrifted.ico6.shape.gii
wb_command -metric-reduce $resultdir/$group_id/groupwise.$group_id.merge.curv.affine.dedrifted.ico6.shape.gii STDEV $resultdir/$group_id/groupwise.$group_id.stdev.curv.affine.dedrifted.ico6.shape.gii

wb_command -metric-merge $resultdir/$group_id/groupwise.$group_id.mean.sulc.curv.affine.dedrifted.ico6.shape.gii -metric $resultdir/$group_id/groupwise.$group_id.mean.sulc.affine.dedrifted.ico6.shape.gii -metric $resultdir/$group_id/groupwise.$group_id.mean.curv.affine.dedrifted.ico6.shape.gii

wb_command -set-structure $resultdir/$group_id/groupwise.$group_id.merge.sulc.affine.dedrifted.ico6.shape.gii CORTEX_LEFT
wb_command -set-structure $resultdir/$group_id/groupwise.$group_id.merge.curv.affine.dedrifted.ico6.shape.gii CORTEX_LEFT
wb_command -set-structure $resultdir/$group_id/groupwise.$group_id.mean.sulc.affine.dedrifted.ico6.shape.gii CORTEX_LEFT
wb_command -set-structure $resultdir/$group_id/groupwise.$group_id.mean.curv.affine.dedrifted.ico6.shape.gii CORTEX_LEFT
wb_command -set-structure $resultdir/$group_id/groupwise.$group_id.stdev.sulc.affine.dedrifted.ico6.shape.gii CORTEX_LEFT
wb_command -set-structure $resultdir/$group_id/groupwise.$group_id.stdev.curv.affine.dedrifted.ico6.shape.gii CORTEX_LEFT
wb_command -set-structure $resultdir/$group_id/groupwise.$group_id.mean.sulc.curv.affine.dedrifted.ico6.shape.gii CORTEX_LEFT
