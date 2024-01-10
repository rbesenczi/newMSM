#!/bin/bash -l

###########################################################
## The following few lines are for running on CREATE cluster. Ignore it otherwise.
#SBATCH --partition=cpu
#SBATCH --time=0-48:00
#SBATCH --nodes=1
#SBATCH --output=./study_2133.txt
#SBATCH --mem=32768
#SBATCH --ntasks=7
#SBATCH --job-name=gMSM
## how to run on create: $ sbatch run_gMSM.sh
###########################################################

###########################################################
## Set the following lines according to your settings and cohort.
workdir=$HOME/groupwise

group_id=2133
group_size=7
subjects=(117122 117123 128026 139637 561444 660951 677968)

input_folder=$workdir/input/$group_id
template=$input_folder/sunet.ico-6.template.surf.gii
config_file=$input_folder/gMSM_config.txt
###########################################################

## Make required files and folders
find $input_folder/data | sort -V > $input_folder/input_data_$group_id.txt
echo "$(tail -n +2 $input_folder/input_data_$group_id.txt)" > $input_folder/input_data_$group_id.txt
find $input_folder/meshes | sort -V > $input_folder/input_meshes_$group_id.txt
echo "$(tail -n +2 $input_folder/input_meshes_$group_id.txt)" > $input_folder/input_meshes_$group_id.txt
mkdir $workdir/output
mkdir $workdir/output/$group_id
mkdir $workdir/results
mkdir $workdir/results/$group_id

## gMSM registration phase ##
echo "Running gMSM for group $group_id"
time $HOME/msm-env/bin/newmsm \
	--data=$input_folder/input_data_$group_id.txt \
	--meshes=$input_folder/input_meshes_$group_id.txt \
	--template=$template \
	--conf=$config_file \
	--out=$workdir/output/$group_id/groupwise.$group_id. \
	--verbose --groupwise

## some renaming and setting structures ##
for (( i=0; i<$group_size; i++ ))
do
	wb_command -set-structure $workdir/output/$group_id/groupwise.$group_id.sphere-$i.reg.surf.gii CORTEX_LEFT
	mv $workdir/output/$group_id/groupwise.$group_id.sphere-$i.reg.surf.gii $workdir/output/$group_id/groupwise.$group_id.sphere-${subjects[$i]}.reg.surf.gii
	wb_command -set-structure $workdir/output/$group_id/groupwise.$group_id.transformed_and_reprojected-$i.func.gii CORTEX_LEFT
	mv $workdir/output/$group_id/groupwise.$group_id.transformed_and_reprojected-$i.func.gii $workdir/output/$group_id/groupwise.$group_id.transformed_and_reprojected-${subjects[$i]}.func.gii
done

## dedrifting phase ##
echo "Dedrifting for group $group_id..."

## calculating surface average of inverse of registration ##
surf_avg="wb_command -surface-average $workdir/output/$group_id/groupwise.$group_id.dedriftwarp.ico-6.sphere.surf.gii "

for subject in "${subjects[@]}"
do
	wb_command -surface-sphere-project-unproject \
	$input_folder/meshes/$subject.surf.gii \
	$workdir/output/$group_id/groupwise.$group_id.sphere-$subject.reg.surf.gii \
	$input_folder/meshes/$subject.surf.gii \
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
	$input_folder/meshes/$subject.surf.gii \
	$workdir/output/$group_id/groupwise.$group_id.dedriftwarp.ico-6.sphere.surf.gii \
	$workdir/output/$group_id/groupwise.$group_id.sphere-$subject.reg.corrected.surf.gii

	## resampling data for the new dedrifted mesh ##
	wb_command -metric-resample \
	$input_folder/data/$subject.L.sulc.affine.ico6.shape.gii \
	$workdir/output/$group_id/groupwise.$group_id.sphere-$subject.reg.corrected.surf.gii \
	$input_folder/meshes/$subject.surf.gii \
	ADAP_BARY_AREA \
	$workdir/output/$group_id/groupwise.$group_id.transformed_and_reprojected.dedrift-$subject.func.gii \
	-area-surfs \
	$workdir/output/$group_id/groupwise.$group_id.sphere-$subject.reg.corrected.surf.gii \
	$input_folder/meshes/$subject.surf.gii
done

echo "Dedrifting for group $group_id done."
echo "Calculating mean and stdev, areal and shape distortion for group $group_id..."

## calculating mean and stdev, areal and shape distortion ##
merge="wb_command -metric-merge $workdir/results/$group_id/groupwise.$group_id.merge.L.sulc.affine.dedrifted.ico6.shape.gii "
arealmerge="wb_command -metric-merge $workdir/results/$group_id/groupwise.$group_id.areal.distortion.merge.L.sulc.affine.dedrifted.ico6.shape.gii "
shapemerge="wb_command -metric-merge $workdir/results/$group_id/groupwise.$group_id.shape.distortion.merge.L.sulc.affine.dedrifted.ico6.shape.gii "

for subject in "${subjects[@]}"
do
	merge+="-metric $workdir/output/$group_id/groupwise.$group_id.transformed_and_reprojected.dedrift-$subject.func.gii "
    wb_command -surface-distortion $input_folder/meshes/$subject.surf.gii $workdir/output/$group_id/groupwise.$group_id.sphere-$subject.reg.corrected.surf.gii $workdir/output/$group_id/groupwise.$group_id.sphere-$subject.distortion.func.gii -local-affine-method -log2
    arealmerge+="-metric $workdir/output/$group_id/groupwise.$group_id.sphere-$subject.distortion.func.gii -column 1 "
    shapemerge+="-metric $workdir/output/$group_id/groupwise.$group_id.sphere-$subject.distortion.func.gii -column 2 "
done

$merge
$arealmerge
$shapemerge

wb_command -metric-reduce $workdir/results/$group_id/groupwise.$group_id.merge.L.sulc.affine.dedrifted.ico6.shape.gii MEAN $workdir/results/$group_id/groupwise.$group_id.mean.L.sulc.affine.dedrifted.ico6.shape.gii
wb_command -metric-reduce $workdir/results/$group_id/groupwise.$group_id.merge.L.sulc.affine.dedrifted.ico6.shape.gii STDEV $workdir/results/$group_id/groupwise.$group_id.stdev.L.sulc.affine.dedrifted.ico6.shape.gii

wb_command -set-structure $workdir/results/$group_id/groupwise.$group_id.merge.L.sulc.affine.dedrifted.ico6.shape.gii CORTEX_LEFT
wb_command -set-structure $workdir/results/$group_id/groupwise.$group_id.mean.L.sulc.affine.dedrifted.ico6.shape.gii CORTEX_LEFT
wb_command -set-structure $workdir/results/$group_id/groupwise.$group_id.stdev.L.sulc.affine.dedrifted.ico6.shape.gii CORTEX_LEFT

rm $workdir/results/$group_id/groupwise.$group_id.merge.L.sulc.affine.dedrifted.ico6.shape.gii	

echo "Calculating mean and stdev, areal and shape distortion for group $group_id done."
