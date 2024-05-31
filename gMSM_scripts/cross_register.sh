#!/bin/bash -l

#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --mem=2048
#SBATCH --job-name=cgMSM

dataset=HCP
workdir=/scratch/users/k2258483/groupwise/${dataset}
outdir=$workdir/output
resultdir=$workdir/results

batch=$workdir/blocks/block_${FILENUM}.txt
clustering=$workdir/frontal_subject_clusters_${dataset}.csv
clustering_out=$workdir/new_clusters/frontal_subject_clusters_HCP_${SLURM_ARRAY_TASK_ID}.csv

while IFS=',' read -r linenum subject group_A_id group_B_id root_node_id
do

  mkdir $outdir/$root_node_id
  mkdir $resultdir/$root_node_id

  if [ $linenum -eq 0 ]
  then
    echo "Registering subject $subject from $group_A_id to $group_B_id template into $root_node_id"

    time $HOME/fsldev/bin/newmsm \
      --indata=$outdir/$group_A_id/groupwise.$group_A_id.transformed_and_reprojected.dedrift-${subject}.func.gii \
      --refdata=$resultdir/$group_B_id/groupwise.$group_B_id.mean.sulc.curv.affine.dedrifted.ico6.shape.gii \
      --inmesh=$workdir/templates/sunet.ico-6.template.surf.gii \
      --refmesh=$workdir/templates/sunet.ico-6.template.surf.gii \
      --conf=$workdir/configs/cgMSM_v2_config.txt \
      --out=$outdir/$root_node_id/groupwise.$root_node_id.${SLURM_ARRAY_TASK_ID}. \
      --verbose
      #--inweight=$workdir/NODE2218_frontal_mask.shape.gii \
      #--refweight=$workdir/NODE2218_frontal_mask.shape.gii \

    mv $outdir/$root_node_id/groupwise.$root_node_id.${SLURM_ARRAY_TASK_ID}.sphere.reg.surf.gii $outdir/$root_node_id/groupwise.$root_node_id.sphere-${subject}.reg.surf.gii
    mv $outdir/$root_node_id/groupwise.$root_node_id.${SLURM_ARRAY_TASK_ID}.transformed_and_reprojected.func.gii $outdir/$root_node_id/groupwise.$root_node_id.transformed_and_reprojected.dedrift-${subject}.func.gii
    wb_command -set-structure $outdir/$root_node_id/groupwise.$root_node_id.sphere-${subject}.reg.surf.gii CORTEX_LEFT
    wb_command -set-structure $outdir/$root_node_id/groupwise.$root_node_id.transformed_and_reprojected.dedrift-${subject}.func.gii CORTEX_LEFT
    rm $outdir/$root_node_id/groupwise.$root_node_id.${SLURM_ARRAY_TASK_ID}.sphere.LR.reg.surf.gii
  fi

  if [ $linenum -eq 1 ]
  then
    subjects_A=()
    subjects_B=()

    while IFS="," read -r linenum subject group
    do
      if [ $group_A_id = $group ]; then
        subjects_A+=($subject)
      fi
      if [ $group_B_id = $group ]; then
        subjects_B+=($subject)
      fi
    done < $clustering

    echo "merge $group_A_id $group_B_id $root_node_id"

    all_subjects=( "${subjects_A[@]}" "${subjects_B[@]}" )

    echo "Calculating merge, mean, stdev, areal and shape distortion for $root_node_id..."

    sulcmerge="wb_command -metric-merge $resultdir/$root_node_id/groupwise.$root_node_id.merge.sulc.affine.dedrifted.ico6.shape.gii "
    curvmerge="wb_command -metric-merge $resultdir/$root_node_id/groupwise.$root_node_id.merge.curv.affine.dedrifted.ico6.shape.gii "
    arealmerge="wb_command -metric-merge $resultdir/$root_node_id/groupwise.$root_node_id.areal.distortion.merge.sulc.affine.dedrifted.ico6.shape.gii "
    shapemerge="wb_command -metric-merge $resultdir/$root_node_id/groupwise.$root_node_id.shape.distortion.merge.sulc.affine.dedrifted.ico6.shape.gii "

    index=0
    for subject in "${all_subjects[@]}"
    do
      echo $subject
      sulcmerge+="-metric $outdir/$root_node_id/groupwise.$root_node_id.transformed_and_reprojected.dedrift-$subject.func.gii -column 1 "
      curvmerge+="-metric $outdir/$root_node_id/groupwise.$root_node_id.transformed_and_reprojected.dedrift-$subject.func.gii -column 2 "
      wb_command -surface-distortion $workdir/templates/sunet.ico-6.template.surf.gii $outdir/$root_node_id/groupwise.$root_node_id.sphere-$subject.reg.surf.gii $outdir/$root_node_id/groupwise.$root_node_id.sphere-$subject.distortion.func.gii -local-affine-method -log2
      arealmerge+="-metric $outdir/$root_node_id/groupwise.$root_node_id.sphere-$subject.distortion.func.gii -column 1 "
      shapemerge+="-metric $outdir/$root_node_id/groupwise.$root_node_id.sphere-$subject.distortion.func.gii -column 2 "
      echo "$index,$subject,$root_node_id" >> $clustering_out
      ((index++))
    done

    $sulcmerge
    $curvmerge
    $arealmerge
    $shapemerge

    wb_command -metric-reduce $resultdir/$root_node_id/groupwise.$root_node_id.merge.sulc.affine.dedrifted.ico6.shape.gii MEAN $resultdir/$root_node_id/groupwise.$root_node_id.mean.sulc.affine.dedrifted.ico6.shape.gii
    wb_command -metric-reduce $resultdir/$root_node_id/groupwise.$root_node_id.merge.curv.affine.dedrifted.ico6.shape.gii MEAN $resultdir/$root_node_id/groupwise.$root_node_id.mean.curv.affine.dedrifted.ico6.shape.gii
    wb_command -metric-reduce $resultdir/$root_node_id/groupwise.$root_node_id.merge.sulc.affine.dedrifted.ico6.shape.gii STDEV $resultdir/$root_node_id/groupwise.$root_node_id.stdev.sulc.affine.dedrifted.ico6.shape.gii
    wb_command -metric-reduce $resultdir/$root_node_id/groupwise.$root_node_id.merge.curv.affine.dedrifted.ico6.shape.gii STDEV $resultdir/$root_node_id/groupwise.$root_node_id.stdev.curv.affine.dedrifted.ico6.shape.gii

    wb_command -set-structure $resultdir/$root_node_id/groupwise.$root_node_id.merge.sulc.affine.dedrifted.ico6.shape.gii CORTEX_LEFT
    wb_command -set-structure $resultdir/$root_node_id/groupwise.$root_node_id.mean.sulc.affine.dedrifted.ico6.shape.gii CORTEX_LEFT
    wb_command -set-structure $resultdir/$root_node_id/groupwise.$root_node_id.stdev.sulc.affine.dedrifted.ico6.shape.gii CORTEX_LEFT
    wb_command -set-structure $resultdir/$root_node_id/groupwise.$root_node_id.merge.curv.affine.dedrifted.ico6.shape.gii CORTEX_LEFT
    wb_command -set-structure $resultdir/$root_node_id/groupwise.$root_node_id.mean.curv.affine.dedrifted.ico6.shape.gii CORTEX_LEFT
    wb_command -set-structure $resultdir/$root_node_id/groupwise.$root_node_id.stdev.curv.affine.dedrifted.ico6.shape.gii CORTEX_LEFT

    wb_command -metric-merge $resultdir/$root_node_id/groupwise.$root_node_id.mean.sulc.curv.affine.dedrifted.ico6.shape.gii -metric $resultdir/$root_node_id/groupwise.$root_node_id.mean.sulc.affine.dedrifted.ico6.shape.gii -metric $resultdir/$root_node_id/groupwise.$root_node_id.mean.curv.affine.dedrifted.ico6.shape.gii
  fi

done < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $batch)
