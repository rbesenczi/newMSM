#/bin/bash

input_folder=$HOME/gpMSM

group_A_id=2126
group_A_size=8
subjects_A=(104820 116221 127226 194645 207628 284646 394956 647858)

group_B_id=2133
group_B_size=7
subjects_B=(117122 117123 128026 139637 561444 660951 677968)

find $input_folder/output/$group_A_id/*reg.corrected* | sort -V > $input_folder/output/$group_A_id/mesh_list_$group_A_id.txt
find $input_folder/output/$group_B_id/*reg.corrected* | sort -V > $input_folder/output/$group_B_id/mesh_list_$group_B_id.txt
mkdir $input_folder/results/cogroup/

time $HOME/fsldev/bin/newmsm \
	--meanA=$input_folder/results/$group_A_id/groupwise.$group_A_id.mean.L.sulc.affine.dedrifted.ico6.shape.gii \
	--meanB=$input_folder/results/$group_B_id/groupwise.$group_B_id.mean.L.sulc.affine.dedrifted.ico6.shape.gii \
	--meshA=$input_folder/input/$group_A_id/sunet.ico-6.template.surf.gii \
	--meshB=$input_folder/input/$group_B_id/sunet.ico-6.template.surf.gii \
	--meshesA=$input_folder/output/$group_A_id/mesh_list_$group_A_id.txt \
	--meshesB=$input_folder/output/$group_B_id/mesh_list_$group_B_id.txt \
	--template=$input_folder/input/$group_A_id/sunet.ico-6.template.surf.gii \
	--conf=$input_folder/input/cogroup_conf.txt \
	--out=$input_folder/results/cogroup/cogroup.$group_A_id.$group_B_id. \
	--verbose --cogroup

for (( i=0; i<$group_A_size; i++ ))
do
	wb_command -set-structure $input_folder/results/cogroup/cogroup.$group_A_id.$group_B_id.sphere-0.$i.reg..surf.gii CORTEX_LEFT
	mv $input_folder/results/cogroup/cogroup.$group_A_id.$group_B_id.sphere-0.$i.reg..surf.gii $input_folder/results/cogroup/cogroup.$group_A_id.$group_B_id.sphere-$group_A_id.${subjects_A[$i]}.reg..surf.gii
done

for (( i=0; i<$group_B_size; i++ ))
do
	wb_command -set-structure $input_folder/results/cogroup/cogroup.$group_A_id.$group_B_id.sphere-1.$i.reg..surf.gii CORTEX_LEFT
	mv $input_folder/results/cogroup/cogroup.$group_A_id.$group_B_id.sphere-1.$i.reg..surf.gii $input_folder/results/cogroup/cogroup.$group_A_id.$group_B_id.sphere-$group_B_id.${subjects_B[$i]}.reg..surf.gii
done

mv $input_folder/results/cogroup/cogroup.$group_A_id.$group_B_id.sphere-0.reg.surf.gii $input_folder/results/cogroup/cogroup.$group_A_id.$group_B_id.sphere-$group_A_id.reg.surf.gii
mv $input_folder/results/cogroup/cogroup.$group_A_id.$group_B_id.sphere-1.reg.surf.gii $input_folder/results/cogroup/cogroup.$group_A_id.$group_B_id.sphere-$group_B_id.reg.surf.gii
wb_command -set-structure $input_folder/results/cogroup/cogroup.$group_A_id.$group_B_id.sphere-$group_A_id.reg.surf.gii CORTEX_LEFT
wb_command -set-structure $input_folder/results/cogroup/cogroup.$group_A_id.$group_B_id.sphere-$group_B_id.reg.surf.gii CORTEX_LEFT

mv $input_folder/results/cogroup/cogroup.$group_A_id.$group_B_id.transformed_and_reprojected-0.func.gii $input_folder/results/cogroup/cogroup.$group_A_id.$group_B_id.transformed_and_reprojected-$group_A_id.func.gii
mv $input_folder/results/cogroup/cogroup.$group_A_id.$group_B_id.transformed_and_reprojected-1.func.gii $input_folder/results/cogroup/cogroup.$group_A_id.$group_B_id.transformed_and_reprojected-$group_B_id.func.gii

wb_command -set-structure $input_folder/results/cogroup/cogroup.$group_A_id.$group_B_id.transformed_and_reprojected-$group_A_id.func.gii CORTEX_LEFT
wb_command -set-structure $input_folder/results/cogroup/cogroup.$group_A_id.$group_B_id.transformed_and_reprojected-$group_B_id.func.gii CORTEX_LEFT

for (( i=0; i<$group_A_size; i++ ))
do
	$HOME/fsldev/bin/applywarp \
	--to_be_deformed=$input_folder/output/$group_A_id/groupwise.$group_A_id.sphere-${subjects_A[$i]}.reg.corrected.surf.gii \
	--warp=$input_folder/results/cogroup/cogroup.$group_A_id.$group_B_id.sphere-$group_A_id.${subjects_A[$i]}.reg..surf.gii \
	--output=$input_folder/results/cogroup/cogroup.$group_A_id.${subjects_A[$i]}.

	wb_command -metric-resample \
	$input_folder/input/$group_A_id/data/${subjects_A[$i]}.L.sulc.affine.ico6.shape.gii \
	$input_folder/results/cogroup/cogroup.$group_A_id.${subjects_A[$i]}.-warped.surf.gii \
	$input_folder/input/$group_A_id/meshes/${subjects_A[$i]}.surf.gii \
	ADAP_BARY_AREA \
	$input_folder/output/$group_A_id/cogroup.$group_A_id.transformed_and_reprojected.dedrift.warped-${subjects_A[$i]}.func.gii \
	-area-surfs \
	$input_folder/results/cogroup/cogroup.$group_A_id.${subjects_A[$i]}.-warped.surf.gii \
	$input_folder/input/$group_A_id/meshes/${subjects_A[$i]}.surf.gii
done

for (( i=0; i<$group_B_size; i++ ))
do
	$HOME/fsldev/bin/applywarp \
	--to_be_deformed=$input_folder/output/$group_B_id/groupwise.$group_B_id.sphere-${subjects_B[$i]}.reg.corrected.surf.gii \
	--warp=$input_folder/results/cogroup/cogroup.$group_A_id.$group_B_id.sphere-$group_B_id.${subjects_B[$i]}.reg..surf.gii \
	--output=$input_folder/results/cogroup/cogroup.$group_B_id.${subjects_B[$i]}.
	
	wb_command -metric-resample \
	$input_folder/input/$group_B_id/data/${subjects_B[$i]}.L.sulc.affine.ico6.shape.gii \
	$input_folder/results/cogroup/cogroup.$group_B_id.${subjects_B[$i]}.-warped.surf.gii \
	$input_folder/input/$group_B_id/meshes/${subjects_B[$i]}.surf.gii \
	ADAP_BARY_AREA \
	$input_folder/output/$group_B_id/cogroup.$group_B_id.transformed_and_reprojected.dedrift.warped-${subjects_B[$i]}.func.gii \
	-area-surfs \
	$input_folder/results/cogroup/cogroup.$group_B_id.${subjects_B[$i]}.-warped.surf.gii \
	$input_folder/input/$group_B_id/meshes/${subjects_B[$i]}.surf.gii

done
