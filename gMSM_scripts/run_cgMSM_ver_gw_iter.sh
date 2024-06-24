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
path=$workdir/frontal_hierarchical_path_study_exp.csv
###########################################################

while IFS=',' read -r group_A group_B group_id
do

  groups=($group_A $group_B)

  for group in "${groups[@]}"
  do
    echo "$group into $group_id"
  done

  mkdir $workdir/file_lists
  rm $workdir/file_lists/input_data_$group_id.txt && touch $workdir/file_lists/input_data_$group_id.txt
  rm $workdir/file_lists/input_meshes_$group_id.txt && touch $workdir/file_lists/input_meshes_$group_id.txt

  for group in "${groups[@]}"
  do
      echo "$resultdir/$group/groupwise.$group.mean.sulc.curv.affine.dedrifted.ico6.shape.gii" >> $workdir/file_lists/input_data_$group_id.txt
      echo "$workdir/templates/sunet.ico-6.template.surf.gii" >> $workdir/file_lists/input_meshes_$group_id.txt
  done

  mkdir $outdir
  mkdir $outdir/$group_id
  mkdir $resultdir
  mkdir $resultdir/$group_id

  ## gMSM registration phase ##
  echo "Running cgMSM for root group $group_id"
  time $HOME/fsldev/bin/newmsm \
    --data=$workdir/file_lists/input_data_$group_id.txt \
    --meshes=$workdir/file_lists/input_meshes_$group_id.txt \
    --template=$template \
    --conf=$config_file \
    --out=$outdir/$group_id/groupwise.$group_id. \
    --verbose --groupwise

  retval=$?
  if [ $retval -ne 0 ]; then
     echo "Error"
     exit $retval
  fi

  index=0
  ## some renaming and setting structures ##
  for group in "${groups[@]}"
  do
    wb_command -set-structure $outdir/$group_id/groupwise.$group_id.sphere-$index.reg.surf.gii CORTEX_LEFT
    mv $outdir/$group_id/groupwise.$group_id.sphere-$index.reg.surf.gii $outdir/$group_id/groupwise.$group_id.sphere-$group.reg.surf.gii
    wb_command -set-structure $outdir/$group_id/groupwise.$group_id.transformed_and_reprojected-$index.func.gii CORTEX_LEFT
    mv $outdir/$group_id/groupwise.$group_id.transformed_and_reprojected-$index.func.gii $outdir/$group_id/groupwise.$group_id.transformed_and_reprojected-$group.func.gii
    ((index++))
  done

  ## dedrifting phase ##
  ## calculating surface average of inverse of registration ##
  surf_avg="wb_command -surface-average $outdir/$group_id/groupwise.$group_id.dedriftwarp.ico-6.sphere.surf.gii "

  for group in "${groups[@]}"
  do
    wb_command -surface-sphere-project-unproject \
    $template \
    $outdir/$group_id/groupwise.$group_id.sphere-$group.reg.surf.gii \
    $template \
    $outdir/$group_id/groupwise.$group_id.sphere-$group.reg.inv.surf.gii
    surf_avg+="-surf $outdir/$group_id/groupwise.$group_id.sphere-$group.reg.inv.surf.gii "
  done

  $surf_avg

  ## recentre and normalise ##
  wb_command -surface-modify-sphere $outdir/$group_id/groupwise.$group_id.dedriftwarp.ico-6.sphere.surf.gii 100 $outdir/$group_id/groupwise.$group_id.dedriftwarp.ico-6.sphere.surf.gii -recenter

  for group in "${groups[@]}"
  do
    ## applying inverse ##
    wb_command -surface-sphere-project-unproject \
    $outdir/$group_id/groupwise.$group_id.sphere-$group.reg.surf.gii \
    $template \
    $outdir/$group_id/groupwise.$group_id.dedriftwarp.ico-6.sphere.surf.gii \
    $outdir/$group_id/groupwise.$group_id.sphere-$group.reg.corrected.surf.gii

    ## resampling data for the new dedrifted mesh ##
    wb_command -metric-resample \
    $resultdir/$group/groupwise.$group.mean.sulc.curv.affine.dedrifted.ico6.shape.gii \
    $outdir/$group_id/groupwise.$group_id.sphere-$group.reg.corrected.surf.gii \
    $template \
    ADAP_BARY_AREA \
    $outdir/$group_id/groupwise.$group_id.transformed_and_reprojected.dedrift-$group.func.gii \
    -area-surfs \
    $outdir/$group_id/groupwise.$group_id.sphere-$group.reg.corrected.surf.gii \
    $template
  done

  ## calculating mean and stdev, areal and shape distortion ##
  merge_sulc="wb_command -metric-merge $resultdir/$group_id/groupwise.$group_id.merge.sulc.affine.dedrifted.ico6.shape.gii "
  merge_curv="wb_command -metric-merge $resultdir/$group_id/groupwise.$group_id.merge.curv.affine.dedrifted.ico6.shape.gii "
  arealmerge="wb_command -metric-merge $resultdir/$group_id/groupwise.$group_id.areal.distortion.merge.sulc.affine.dedrifted.ico6.shape.gii "
  shapemerge="wb_command -metric-merge $resultdir/$group_id/groupwise.$group_id.shape.distortion.merge.sulc.affine.dedrifted.ico6.shape.gii "

  for group in "${groups[@]}"
  do
    merge_sulc+="-metric $outdir/$group_id/groupwise.$group_id.transformed_and_reprojected.dedrift-$group.func.gii -column 1 "
    merge_curv+="-metric $outdir/$group_id/groupwise.$group_id.transformed_and_reprojected.dedrift-$group.func.gii -column 2 "
    wb_command -surface-distortion $template $outdir/$group_id/groupwise.$group_id.sphere-$group.reg.corrected.surf.gii $outdir/$group_id/groupwise.$group_id.sphere-$group.distortion.func.gii -local-affine-method -log2
    arealmerge+="-metric $outdir/$group_id/groupwise.$group_id.sphere-$group.distortion.func.gii -column 1 "
    shapemerge+="-metric $outdir/$group_id/groupwise.$group_id.sphere-$group.distortion.func.gii -column 2 "
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

  ## calculating mean and stdev, areal and shape distortion ##
  merge_sulc="wb_command -metric-merge $resultdir/$group_id/groupwise.$group_id.merge.sulc.affine.dedrifted.ico6.shape.gii "
  merge_curv="wb_command -metric-merge $resultdir/$group_id/groupwise.$group_id.merge.curv.affine.dedrifted.ico6.shape.gii "
  arealmerge="wb_command -metric-merge $resultdir/$group_id/groupwise.$group_id.areal.distortion.merge.sulc.affine.dedrifted.ico6.shape.gii "
  shapemerge="wb_command -metric-merge $resultdir/$group_id/groupwise.$group_id.shape.distortion.merge.sulc.affine.dedrifted.ico6.shape.gii "

  while IFS="," read -r linenum subject group
  do
    if [[ ${groups[@]} =~ $group ]]; then

      echo "$subject"

      #wb_command -surface-sphere-project-unproject \
      #$outdir/$group/groupwise.$group.sphere-$subject.reg.corrected.surf.gii \
      #$outdir/$group/groupwise.$group.sphere-$subject.reg.corrected.surf.gii \
      #$outdir/$group_id/groupwise.$group_id.sphere-$group.reg.corrected.surf.gii \
      #$outdir/$group_id/groupwise.$group_id.sphere-$subject.reg.corrected.surf.gii

      #wb_command -metric-resample \
      #$outdir/$group/groupwise.$group.transformed_and_reprojected.dedrift-$subject.func.gii \
      #$outdir/$group_id/groupwise.$group_id.sphere-$subject.reg.corrected.surf.gii \
      #$template \
      #ADAP_BARY_AREA \
      #$outdir/$group_id/groupwise.$group_id.transformed_and_reprojected.dedrift-$subject.func.gii \
      #-area-surfs \
      #$outdir/$group_id/groupwise.$group_id.sphere-$subject.reg.corrected.surf.gii \
      #$template

      wb_command -surface-sphere-project-unproject \
	      $outdir/$group/groupwise.$group.sphere-$subject.reg.corrected.surf.gii \
	      $template \
	      $outdir/$group_id/groupwise.$group_id.sphere-$group.reg.corrected.surf.gii \
	      $outdir/$group_id/groupwise.$group_id.sphere-$subject.reg.corrected.surf.gii

      wb_command -metric-resample \
	      $outdir/$group/groupwise.$group.transformed_and_reprojected.dedrift-$subject.func.gii \
	      $outdir/$group/groupwise.$group.sphere-$subject.reg.corrected.surf.gii \
  	      $outdir/$group_id/groupwise.$group_id.sphere-$subject.reg.corrected.surf.gii \
 	      ADAP_BARY_AREA \
  	      $outdir/$group_id/groupwise.$group_id.transformed_and_reprojected.dedrift-$subject.func.gii \
  	      -area-surfs \
    	      $outdir/$group/groupwise.$group.sphere-$subject.reg.corrected.surf.gii \
 	      $outdir/$group_id/groupwise.$group_id.sphere-$subject.reg.corrected.surf.gii


      merge_sulc+="-metric $outdir/$group_id/groupwise.$group_id.transformed_and_reprojected.dedrift-$subject.func.gii -column 1 "
      merge_curv+="-metric $outdir/$group_id/groupwise.$group_id.transformed_and_reprojected.dedrift-$subject.func.gii -column 2 "
      wb_command -surface-distortion $template $outdir/$group_id/groupwise.$group_id.sphere-$subject.reg.corrected.surf.gii $outdir/$group_id/groupwise.$group_id.sphere-$subject.distortion.func.gii -local-affine-method -log2
      arealmerge+="-metric $outdir/$group_id/groupwise.$group_id.sphere-$subject.distortion.func.gii -column 1 "
      shapemerge+="-metric $outdir/$group_id/groupwise.$group_id.sphere-$subject.distortion.func.gii -column 2 "

      #echo "0,$subject,$group_id" >> $clustering
      fi
  done < $clustering

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

done < $path
