#! /bin/bash -l

###########################################################
## The following few lines are for running on CREATE cluster. Ignore it otherwise.
#SBATCH --partition=cpu
#SBATCH --time=0-48:00
#SBATCH --nodes=1
#SBATCH --mem=65536
#SBATCH --ntasks=16
#SBATCH --job-name=gMSM
#SBATCH --output=./gMSM_log.txt
## how to run on create: $ sbatch run_gMSM.sh
###########################################################

###########################################################
## Set the following lines according to your settings and cohort.
dataset=HCP
workdir=$HOME/groupwise/$dataset
input_folder=$HOME/affined_${dataset}_sulc
template=$workdir/templates/sunet.ico-6.template.surf.gii
config_file=$workdir/configs/gMSM_${dataset}_config.txt
###########################################################
groups=(NODE1847 NODE1917)

for g in "${groups[@]}"
do
  group_id=${g}_${dataset}

  subject_list=$workdir/group_lists/$group_id.csv

  rm $workdir/file_lists/input_data_$group_id.txt && touch $workdir/file_lists/input_data_$group_id.txt
  rm $workdir/file_lists/input_meshes_$group_id.txt && touch $workdir/file_lists/input_meshes_$group_id.txt

  while IFS="," read -r rec_column1 rec_column2
  do
    echo "$input_folder/$rec_column1.sulc.affine.ico6.shape.gii" >> $workdir/file_lists/input_data_$group_id.txt
    echo "$workdir/templates/sunet.ico-6.template.surf.gii" >> $workdir/file_lists/input_meshes_$group_id.txt
  done < $subject_list

  subjects=( $(cat $subject_list | cut -d ',' -f1) )

  mkdir $workdir/output
  mkdir $workdir/output/$group_id
  mkdir $workdir/results
  mkdir $workdir/results/$group_id

  ## gMSM registration phase ##
  echo "Running gMSM for group $group_id"
  time $HOME/fsldev/bin/newmsm \
    --data=$workdir/file_lists/input_data_$group_id.txt \
    --meshes=$workdir/file_lists/input_meshes_$group_id.txt \
    --template=$template \
    --conf=$config_file \
    --out=$workdir/output/$group_id/groupwise.$group_id. \
    --verbose --groupwise

  index=0
  ## some renaming and setting structures ##
  for subject in "${subjects[@]}"
  do
    wb_command -set-structure $workdir/output/$group_id/groupwise.$group_id.sphere-$index.reg.surf.gii CORTEX_LEFT
    mv $workdir/output/$group_id/groupwise.$group_id.sphere-$index.reg.surf.gii $workdir/output/$group_id/groupwise.$group_id.sphere-$subject.reg.surf.gii
    wb_command -set-structure $workdir/output/$group_id/groupwise.$group_id.transformed_and_reprojected-$index.func.gii CORTEX_LEFT
    mv $workdir/output/$group_id/groupwise.$group_id.transformed_and_reprojected-$index.func.gii $workdir/output/$group_id/groupwise.$group_id.transformed_and_reprojected-$subject.func.gii
    ((index++))
  done

  ## dedrifting phase ##
  echo "Dedrifting for group $group_id..."

  ## calculating surface average of inverse of registration ##
  surf_avg="wb_command -surface-average $workdir/output/$group_id/groupwise.$group_id.dedriftwarp.ico-6.sphere.surf.gii "

  for subject in "${subjects[@]}"
  do
    wb_command -surface-sphere-project-unproject \
    $template \
    $workdir/output/$group_id/groupwise.$group_id.sphere-$subject.reg.surf.gii \
    $template \
    $workdir/output/$group_id/groupwise.$group_id.sphere-$subject.reg.inv.surf.gii
    surf_avg+="-surf $workdir/output/$group_id/groupwise.$group_id.sphere-$subject.reg.inv.surf.gii "
  done

  $surf_avg

  ## recentre and normalise ##
  wb_command -surface-modify-sphere $workdir/output/$group_id/groupwise.$group_id.dedriftwarp.ico-6.sphere.surf.gii 100 $workdir/output/$group_id/groupwise.$group_id.dedriftwarp.ico-6.sphere.surf.gii -recenter

  for subject in "${subjects[@]}"
  do
    ## applying inverse ##
    wb_command -surface-sphere-project-unproject \
    $workdir/output/$group_id/groupwise.$group_id.sphere-$subject.reg.surf.gii \
    $template \
    $workdir/output/$group_id/groupwise.$group_id.dedriftwarp.ico-6.sphere.surf.gii \
    $workdir/output/$group_id/groupwise.$group_id.sphere-$subject.reg.corrected.surf.gii

    ## resampling data for the new dedrifted mesh ##
    wb_command -metric-resample \
    $input_folder/$subject.sulc.affine.ico6.shape.gii \
    $workdir/output/$group_id/groupwise.$group_id.sphere-$subject.reg.corrected.surf.gii \
    $template \
    ADAP_BARY_AREA \
    $workdir/output/$group_id/groupwise.$group_id.transformed_and_reprojected.dedrift-$subject.func.gii \
    -area-surfs \
    $workdir/output/$group_id/groupwise.$group_id.sphere-$subject.reg.corrected.surf.gii \
    $template
  done

  echo "Dedrifting for group $group_id done."
  echo "Calculating mean and stdev, areal and shape distortion for group $group_id..."

  ## calculating mean and stdev, areal and shape distortion ##
  merge="wb_command -metric-merge $workdir/results/$group_id/groupwise.$group_id.merge.sulc.affine.dedrifted.ico6.shape.gii "
  arealmerge="wb_command -metric-merge $workdir/results/$group_id/groupwise.$group_id.areal.distortion.merge.sulc.affine.dedrifted.ico6.shape.gii "
  shapemerge="wb_command -metric-merge $workdir/results/$group_id/groupwise.$group_id.shape.distortion.merge.sulc.affine.dedrifted.ico6.shape.gii "

  for subject in "${subjects[@]}"
  do
    merge+="-metric $workdir/output/$group_id/groupwise.$group_id.transformed_and_reprojected.dedrift-$subject.func.gii "
    wb_command -surface-distortion $template $workdir/output/$group_id/groupwise.$group_id.sphere-$subject.reg.corrected.surf.gii $workdir/output/$group_id/groupwise.$group_id.sphere-$subject.distortion.func.gii -local-affine-method -log2
    arealmerge+="-metric $workdir/output/$group_id/groupwise.$group_id.sphere-$subject.distortion.func.gii -column 1 "
    shapemerge+="-metric $workdir/output/$group_id/groupwise.$group_id.sphere-$subject.distortion.func.gii -column 2 "
  done

  $merge
  $arealmerge
  $shapemerge

  wb_command -metric-reduce $workdir/results/$group_id/groupwise.$group_id.merge.sulc.affine.dedrifted.ico6.shape.gii MEAN $workdir/results/$group_id/groupwise.$group_id.mean.sulc.affine.dedrifted.ico6.shape.gii
  wb_command -metric-reduce $workdir/results/$group_id/groupwise.$group_id.merge.sulc.affine.dedrifted.ico6.shape.gii STDEV $workdir/results/$group_id/groupwise.$group_id.stdev.sulc.affine.dedrifted.ico6.shape.gii

  wb_command -set-structure $workdir/results/$group_id/groupwise.$group_id.merge.sulc.affine.dedrifted.ico6.shape.gii CORTEX_LEFT
  wb_command -set-structure $workdir/results/$group_id/groupwise.$group_id.mean.sulc.affine.dedrifted.ico6.shape.gii CORTEX_LEFT
  wb_command -set-structure $workdir/results/$group_id/groupwise.$group_id.stdev.sulc.affine.dedrifted.ico6.shape.gii CORTEX_LEFT

  rm $workdir/results/$group_id/groupwise.$group_id.merge.sulc.affine.dedrifted.ico6.shape.gii

  echo "Calculating mean and stdev, areal and shape distortion for group $group_id done."

done
