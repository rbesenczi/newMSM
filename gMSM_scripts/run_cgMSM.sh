#/bin/bash

dataset=HCP
workdir=$HOME/groupwise/$dataset

mkdir $workdir/lists

while IFS="," read -r group_A_id group_B_id root_id
do
	group_A_id=${group_A_id}_HPC
	group_B_id=${group_B_id}_HPC
	root_node_id=${root_id}_HPC
	group_A_list=$workdir/group_lists/$group_A_id.csv
	group_B_list=$workdir/group_lists/$group_B_id.csv

	subjects_A=( $(cat $group_A_list | cut -d ',' -f1) )
	subjects_B=( $(cat $group_B_list | cut -d ',' -f1) )

	rm $workdir/lists/mesh_list_$group_A_id.txt
	rm $workdir/lists/mesh_list_$group_B_id.txt

	for subject in "${subjects_A[@]}"
	do
		echo "$workdir/output/$group_A_id/groupwise.$group_A_id.sphere-$subject.reg.corrected.surf.gii" >> $workdir/lists/mesh_list_$group_A_id.txt
	done

	for subject in "${subjects_B[@]}"
	do
		echo "$workdir/output/$group_B_id/groupwise.$group_B_id.sphere-$subject.reg.corrected.surf.gii" >> $workdir/lists/mesh_list_$group_B_id.txt
	done

	mkdir $workdir/output
	mkdir $workdir/output/$root_node_id
	mkdir $workdir/results
	mkdir $workdir/results/$root_node_id

	time $HOME/fsldev/bin/newmsm \
		--meanA=$workdir/results/$group_A_id/groupwise.$group_A_id.mean.sulc.affine.dedrifted.ico6.shape.gii \
		--meanB=$workdir/results/$group_B_id/groupwise.$group_B_id.mean.sulc.affine.dedrifted.ico6.shape.gii \
		--meshA=$workdir/sunet.ico-6.template.surf.gii \
		--meshB=$workdir/sunet.ico-6.template.surf.gii \
		--meshesA=$workdir/lists/mesh_list_$group_A_id.txt \
		--meshesB=$workdir/lists/mesh_list_$group_B_id.txt \
		--template=$workdir/sunet.ico-6.template.surf.gii \
		--conf=$workdir/cgMSM_HCP_config.txt \
		--out=$workdir/output/$root_node_id/groupwise.$root_node_id. \
		--verbose --cogroup

	if [ $? -ne 0 ]; then
		exit 1
	fi

	index=0
	for subject in "${subjects_A[@]}"
	do
		wb_command -set-structure $workdir/output/$root_node_id/groupwise.$root_node_id.sphere-0.$index.reg.surf.gii CORTEX_LEFT
		mv $workdir/output/$root_node_id/groupwise.$root_node_id.sphere-0.$index.reg.surf.gii $workdir/output/$root_node_id/groupwise.$root_node_id.sphere.$subject.deformation.surf.gii
		((index++))
	done

	index=0
	for subject in "${subjects_B[@]}"
	do
		wb_command -set-structure $workdir/output/$root_node_id/groupwise.$root_node_id.sphere-1.$index.reg.surf.gii CORTEX_LEFT
		mv $workdir/output/$root_node_id/groupwise.$root_node_id.sphere-1.$index.reg.surf.gii $workdir/output/$root_node_id/groupwise.$root_node_id.sphere.$subject.deformation.surf.gii
		((index++))
	done

	mv $workdir/output/$root_node_id/groupwise.$root_node_id.sphere-0.reg.surf.gii $workdir/output/$root_node_id/groupwise.$root_node_id.sphere-$group_A_id.reg.surf.gii
	mv $workdir/output/$root_node_id/groupwise.$root_node_id.sphere-1.reg.surf.gii $workdir/output/$root_node_id/groupwise.$root_node_id.sphere-$group_B_id.reg.surf.gii
	wb_command -set-structure $workdir/output/$root_node_id/groupwise.$root_node_id.sphere-$group_A_id.reg.surf.gii CORTEX_LEFT
	wb_command -set-structure $workdir/output/$root_node_id/groupwise.$root_node_id.sphere-$group_B_id.reg.surf.gii CORTEX_LEFT

	mv $workdir/output/$root_node_id/groupwise.$root_node_id.transformed_and_reprojected-0.func.gii $workdir/output/$root_node_id/groupwise.$root_node_id.transformed_and_reprojected-$group_A_id.func.gii
	mv $workdir/output/$root_node_id/groupwise.$root_node_id.transformed_and_reprojected-1.func.gii $workdir/output/$root_node_id/groupwise.$root_node_id.transformed_and_reprojected-$group_B_id.func.gii
	wb_command -set-structure $workdir/output/$root_node_id/groupwise.$root_node_id.transformed_and_reprojected-$group_A_id.func.gii CORTEX_LEFT
	wb_command -set-structure $workdir/output/$root_node_id/groupwise.$root_node_id.transformed_and_reprojected-$group_B_id.func.gii CORTEX_LEFT

	for subject in "${subjects_A[@]}"
	do
		echo "Warping and resampling $group_A_id $subject."

		$HOME/fsldev/bin/applywarp \
		--to_be_deformed=$workdir/output/$group_A_id/groupwise.$group_A_id.sphere-$subject.reg.corrected.surf.gii \
		--warp=$workdir/output/$root_node_id/groupwise.$root_node_id.sphere.$subject.deformation.surf.gii \
		--output=$workdir/output/$root_node_id/groupwise.$root_node_id.sphere-$subject.

		mv $workdir/output/$root_node_id/groupwise.$root_node_id.sphere-$subject.warped.surf.gii $workdir/output/$root_node_id/groupwise.$root_node_id.sphere-$subject.reg.corrected.surf.gii

		wb_command -metric-resample \
		$workdir/output/$group_A_id/groupwise.$group_A_id.transformed_and_reprojected.dedrift-$subject.func.gii \
		$workdir/output/$root_node_id/groupwise.$root_node_id.sphere-$subject.reg.corrected.surf.gii \
		$workdir/output/$group_A_id/groupwise.$group_A_id.sphere-$subject.reg.corrected.surf.gii \
		ADAP_BARY_AREA \
		$workdir/output/$root_node_id/groupwise.$root_node_id.transformed_and_reprojected.dedrift-$subject.func.gii \
		-area-surfs \
		$workdir/output/$root_node_id/groupwise.$root_node_id.sphere-$subject.reg.corrected.surf.gii \
		$workdir/output/$group_A_id/groupwise.$group_A_id.sphere-$subject.reg.corrected.surf.gii
	done

	for subject in "${subjects_B[@]}"
	do
		echo "Warping and resampling $group_B_id $subject."

		$HOME/fsldev/bin/applywarp \
		--to_be_deformed=$workdir/output/$group_B_id/groupwise.$group_B_id.sphere-$subject.reg.corrected.surf.gii \
		--warp=$workdir/output/$root_node_id/groupwise.$root_node_id.sphere.$subject.deformation.surf.gii \
		--output=$workdir/output/$root_node_id/groupwise.$root_node_id.sphere-$subject.

		mv $workdir/output/$root_node_id/groupwise.$root_node_id.sphere-$subject.warped.surf.gii $workdir/output/$root_node_id/groupwise.$root_node_id.sphere-$subject.reg.corrected.surf.gii

		wb_command -metric-resample \
		$workdir/output/$group_B_id/groupwise.$group_B_id.transformed_and_reprojected.dedrift-$subject.func.gii \
		$workdir/output/$root_node_id/groupwise.$root_node_id.sphere-$subject.reg.corrected.surf.gii \
		$workdir/output/$group_B_id/groupwise.$group_B_id.sphere-$subject.reg.corrected.surf.gii \
		ADAP_BARY_AREA \
		$workdir/output/$root_node_id/groupwise.$root_node_id.transformed_and_reprojected.dedrift-$subject.func.gii \
		-area-surfs \
		$workdir/output/$root_node_id/groupwise.$root_node_id.sphere-$subject.reg.corrected.surf.gii \
		$workdir/output/$group_B_id/groupwise.$group_B_id.sphere-$subject.reg.corrected.surf.gii
	done

	all_subjects=( "${subjects_A[@]}" "${subjects_B[@]}" )

	echo "calculating mean and stdev, areal and shape distortion"

	merge="wb_command -metric-merge $workdir/results/$root_node_id/groupwise.$root_node_id.merge.sulc.affine.dedrifted.ico6.shape.gii "
	arealmerge="wb_command -metric-merge $workdir/results/$root_node_id/groupwise.$root_node_id.areal.distortion.merge.sulc.affine.dedrifted.ico6.shape.gii "
	shapemerge="wb_command -metric-merge $workdir/results/$root_node_id/groupwise.$root_node_id.shape.distortion.merge.sulc.affine.dedrifted.ico6.shape.gii "

	for subject in "${all_subjects[@]}"
	do
		merge+="-metric $workdir/output/$root_node_id/groupwise.$root_node_id.transformed_and_reprojected.dedrift-$subject.func.gii "
		wb_command -surface-distortion $workdir/sunet.ico-6.template.surf.gii $workdir/output/$root_node_id/groupwise.$root_node_id.sphere-$subject.reg.corrected.surf.gii $workdir/output/$root_node_id/groupwise.$root_node_id.sphere-$subject.distortion.func.gii -local-affine-method -log2
		arealmerge+="-metric $workdir/output/$root_node_id/groupwise.$root_node_id.sphere-$subject.distortion.func.gii -column 1 "
		shapemerge+="-metric $workdir/output/$root_node_id/groupwise.$root_node_id.sphere-$subject.distortion.func.gii -column 2 "
	done

	$merge
	$arealmerge
	$shapemerge

	wb_command -metric-reduce $workdir/results/$root_node_id/groupwise.$root_node_id.merge.sulc.affine.dedrifted.ico6.shape.gii MEAN $workdir/results/$root_node_id/groupwise.$root_node_id.mean.sulc.affine.dedrifted.ico6.shape.gii
	wb_command -metric-reduce $workdir/results/$root_node_id/groupwise.$root_node_id.merge.sulc.affine.dedrifted.ico6.shape.gii STDEV $workdir/results/$root_node_id/groupwise.$root_node_id.stdev.sulc.affine.dedrifted.ico6.shape.gii

	wb_command -set-structure $workdir/results/$root_node_id/groupwise.$root_node_id.merge.sulc.affine.dedrifted.ico6.shape.gii CORTEX_LEFT
	wb_command -set-structure $workdir/results/$root_node_id/groupwise.$root_node_id.mean.sulc.affine.dedrifted.ico6.shape.gii CORTEX_LEFT
	wb_command -set-structure $workdir/results/$root_node_id/groupwise.$root_node_id.stdev.sulc.affine.dedrifted.ico6.shape.gii CORTEX_LEFT

done < $HOME/groupwise/data/frontal_hierarchical_path_rem_S.csv
