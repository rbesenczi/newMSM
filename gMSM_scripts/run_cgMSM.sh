#! /bin/bash -l

###########################################################
## The following few lines are for running on CREATE cluster. Ignore it otherwise.
#SBATCH --partition=cpu
#SBATCH --time=0-48:00
#SBATCH --nodes=1
#SBATCH --mem=65536
#SBATCH --ntasks=16
#SBATCH --job-name=cgMSM
#SBATCH --output=./cgMSM_log.txt
## how to run on create: $ sbatch run_cgMSM.sh
###########################################################

###########################################################
## Set the following lines according to your settings and cohort.
dataset=HCP
workdir=$HOME/groupwise/$dataset
outdir=$workdir/output
resultdir=$workdir/results

clustering=$workdir/frontal_subject_clusters_hcp.csv
hierarchy=$workdir/frontal_hierarchical_path_study.csv
###########################################################

mkdir $outdir
mkdir $resultdir

while IFS="," read -r group_A_id group_B_id root_node_id
do
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

	rm $workdir/file_lists/mesh_list_$group_A_id.txt
	rm $workdir/file_lists/mesh_list_$group_B_id.txt

	for subject in "${subjects_A[@]}"
	do
		echo "$outdir/$group_A_id/groupwise.$group_A_id.sphere-$subject.reg.corrected.surf.gii" >> $workdir/file_lists/mesh_list_$group_A_id.txt
	done

	for subject in "${subjects_B[@]}"
	do
		echo "$outdir/$group_B_id/groupwise.$group_B_id.sphere-$subject.reg.corrected.surf.gii" >> $workdir/file_lists/mesh_list_$group_B_id.txt
	done

	mkdir $outdir/$root_node_id
	mkdir $resultdir/$root_node_id

	time $HOME/fsldev/bin/newmsm \
		--meanA=$resultdir/$group_A_id/groupwise.$group_A_id.mean.sulc.affine.dedrifted.ico6.shape.gii \
		--meanB=$resultdir/$group_B_id/groupwise.$group_B_id.mean.sulc.affine.dedrifted.ico6.shape.gii \
		--meshA=$workdir/templates/sunet.ico-6.template.surf.gii \
		--meshB=$workdir/templates/sunet.ico-6.template.surf.gii \
		--meshesA=$workdir/file_lists/mesh_list_$group_A_id.txt \
		--meshesB=$workdir/file_lists/mesh_list_$group_B_id.txt \
		--template=$workdir/templates/sunet.ico-6.template.surf.gii \
		--conf=$workdir/configs/cgMSM_config.txt \
		--out=$outdir/$root_node_id/groupwise.$root_node_id. \
		--verbose --cogroup

	if [ $? -ne 0 ]; then
		exit 1
	fi

	index=0
	for subject in "${subjects_A[@]}"
	do
		wb_command -set-structure $outdir/$root_node_id/groupwise.$root_node_id.sphere-0.$index.reg.surf.gii CORTEX_LEFT
		mv $outdir/$root_node_id/groupwise.$root_node_id.sphere-0.$index.reg.surf.gii $outdir/$root_node_id/groupwise.$root_node_id.sphere.$subject.deformation.surf.gii
		((index++))
	done

	index=0
	for subject in "${subjects_B[@]}"
	do
		wb_command -set-structure $outdir/$root_node_id/groupwise.$root_node_id.sphere-1.$index.reg.surf.gii CORTEX_LEFT
		mv $outdir/$root_node_id/groupwise.$root_node_id.sphere-1.$index.reg.surf.gii $outdir/$root_node_id/groupwise.$root_node_id.sphere.$subject.deformation.surf.gii
		((index++))
	done

	mv $outdir/$root_node_id/groupwise.$root_node_id.sphere-0.reg.surf.gii $outdir/$root_node_id/groupwise.$root_node_id.sphere-$group_A_id.reg.surf.gii
	mv $outdir/$root_node_id/groupwise.$root_node_id.sphere-1.reg.surf.gii $outdir/$root_node_id/groupwise.$root_node_id.sphere-$group_B_id.reg.surf.gii
	wb_command -set-structure $outdir/$root_node_id/groupwise.$root_node_id.sphere-$group_A_id.reg.surf.gii CORTEX_LEFT
	wb_command -set-structure $outdir/$root_node_id/groupwise.$root_node_id.sphere-$group_B_id.reg.surf.gii CORTEX_LEFT

	mv $outdir/$root_node_id/groupwise.$root_node_id.transformed_and_reprojected-0.func.gii $outdir/$root_node_id/groupwise.$root_node_id.transformed_and_reprojected-$group_A_id.func.gii
	mv $outdir/$root_node_id/groupwise.$root_node_id.transformed_and_reprojected-1.func.gii $outdir/$root_node_id/groupwise.$root_node_id.transformed_and_reprojected-$group_B_id.func.gii
	wb_command -set-structure $outdir/$root_node_id/groupwise.$root_node_id.transformed_and_reprojected-$group_A_id.func.gii CORTEX_LEFT
	wb_command -set-structure $outdir/$root_node_id/groupwise.$root_node_id.transformed_and_reprojected-$group_B_id.func.gii CORTEX_LEFT

	for subject in "${subjects_A[@]}"
	do
		echo "Warping and resampling $group_A_id $subject."

		$HOME/fsldev/bin/applywarp \
		--to_be_deformed=$outdir/$group_A_id/groupwise.$group_A_id.sphere-$subject.reg.corrected.surf.gii \
		--warp=$outdir/$root_node_id/groupwise.$root_node_id.sphere.$subject.deformation.surf.gii \
		--output=$outdir/$root_node_id/groupwise.$root_node_id.sphere-$subject.

		mv $outdir/$root_node_id/groupwise.$root_node_id.sphere-$subject.warped.surf.gii $outdir/$root_node_id/groupwise.$root_node_id.sphere-$subject.reg.corrected.surf.gii

		wb_command -metric-resample \
		$outdir/$group_A_id/groupwise.$group_A_id.transformed_and_reprojected.dedrift-$subject.func.gii \
		$outdir/$root_node_id/groupwise.$root_node_id.sphere-$subject.reg.corrected.surf.gii \
		$outdir/$group_A_id/groupwise.$group_A_id.sphere-$subject.reg.corrected.surf.gii \
		ADAP_BARY_AREA \
		$outdir/$root_node_id/groupwise.$root_node_id.transformed_and_reprojected.dedrift-$subject.func.gii \
		-area-surfs \
		$outdir/$root_node_id/groupwise.$root_node_id.sphere-$subject.reg.corrected.surf.gii \
		$outdir/$group_A_id/groupwise.$group_A_id.sphere-$subject.reg.corrected.surf.gii
	done

	for subject in "${subjects_B[@]}"
	do
		echo "Warping and resampling $group_B_id $subject."

		$HOME/fsldev/bin/applywarp \
		--to_be_deformed=$outdir/$group_B_id/groupwise.$group_B_id.sphere-$subject.reg.corrected.surf.gii \
		--warp=$outdir/$root_node_id/groupwise.$root_node_id.sphere.$subject.deformation.surf.gii \
		--output=$outdir/$root_node_id/groupwise.$root_node_id.sphere-$subject.

		mv $outdir/$root_node_id/groupwise.$root_node_id.sphere-$subject.warped.surf.gii $outdir/$root_node_id/groupwise.$root_node_id.sphere-$subject.reg.corrected.surf.gii

		wb_command -metric-resample \
		$outdir/$group_B_id/groupwise.$group_B_id.transformed_and_reprojected.dedrift-$subject.func.gii \
		$outdir/$root_node_id/groupwise.$root_node_id.sphere-$subject.reg.corrected.surf.gii \
		$outdir/$group_B_id/groupwise.$group_B_id.sphere-$subject.reg.corrected.surf.gii \
		ADAP_BARY_AREA \
		$outdir/$root_node_id/groupwise.$root_node_id.transformed_and_reprojected.dedrift-$subject.func.gii \
		-area-surfs \
		$outdir/$root_node_id/groupwise.$root_node_id.sphere-$subject.reg.corrected.surf.gii \
		$outdir/$group_B_id/groupwise.$group_B_id.sphere-$subject.reg.corrected.surf.gii
	done

	all_subjects=( "${subjects_A[@]}" "${subjects_B[@]}" )

	echo "Calculating merge, mean, stdev, areal and shape distortion for $root_node_id..."

	merge="wb_command -metric-merge $resultdir/$root_node_id/groupwise.$root_node_id.merge.sulc.affine.dedrifted.ico6.shape.gii "
	arealmerge="wb_command -metric-merge $resultdir/$root_node_id/groupwise.$root_node_id.areal.distortion.merge.sulc.affine.dedrifted.ico6.shape.gii "
	shapemerge="wb_command -metric-merge $resultdir/$root_node_id/groupwise.$root_node_id.shape.distortion.merge.sulc.affine.dedrifted.ico6.shape.gii "

	index=0
	for subject in "${all_subjects[@]}"
	do
		merge+="-metric $outdir/$root_node_id/groupwise.$root_node_id.transformed_and_reprojected.dedrift-$subject.func.gii "
		wb_command -surface-distortion $workdir/templates/sunet.ico-6.template.surf.gii $outdir/$root_node_id/groupwise.$root_node_id.sphere-$subject.reg.corrected.surf.gii $outdir/$root_node_id/groupwise.$root_node_id.sphere-$subject.distortion.func.gii -local-affine-method -log2
		arealmerge+="-metric $outdir/$root_node_id/groupwise.$root_node_id.sphere-$subject.distortion.func.gii -column 1 "
		shapemerge+="-metric $outdir/$root_node_id/groupwise.$root_node_id.sphere-$subject.distortion.func.gii -column 2 "
		echo "$index,$subject,$root_node_id" >> $clustering
		((index++))
	done

	$merge
	$arealmerge
	$shapemerge

	wb_command -metric-reduce $resultdir/$root_node_id/groupwise.$root_node_id.merge.sulc.affine.dedrifted.ico6.shape.gii MEAN $resultdir/$root_node_id/groupwise.$root_node_id.mean.sulc.affine.dedrifted.ico6.shape.gii
	wb_command -metric-reduce $resultdir/$root_node_id/groupwise.$root_node_id.merge.sulc.affine.dedrifted.ico6.shape.gii STDEV $resultdir/$root_node_id/groupwise.$root_node_id.stdev.sulc.affine.dedrifted.ico6.shape.gii

	wb_command -set-structure $resultdir/$root_node_id/groupwise.$root_node_id.merge.sulc.affine.dedrifted.ico6.shape.gii CORTEX_LEFT
	wb_command -set-structure $resultdir/$root_node_id/groupwise.$root_node_id.mean.sulc.affine.dedrifted.ico6.shape.gii CORTEX_LEFT
	wb_command -set-structure $resultdir/$root_node_id/groupwise.$root_node_id.stdev.sulc.affine.dedrifted.ico6.shape.gii CORTEX_LEFT

done < $hierarchy
