#!/bin/bash 

###########################################################
## Set the following lines according to your settings and cohort.
nthreads=8
parallel_tasks=$(( 350/$nthreads ))

dataset="HCP"
workdir=$HOME/${dataset}_to_template
input_folder=$HOME/affined_${dataset}_sulc
template=$workdir/MSMSulc.L.sulc.ico6.shape.gii
mesh=$workdir/sunet.ico-6.sphere.surf.gii
config_file=$workdir/HCP_to_template_config.txt
sublist=$workdir/subjects_in_study.txt
outdir=$workdir/output
###########################################################

num_subjects=$(cat $sublist | wc -l)
bunch_size=$(( $num_subjects / $parallel_tasks + 1 ))

touch $workdir/newMSM_HCP_to_template.sh
mkdir $outdir
mkdir $workdir/logs

for (( bunch=0; bunch<parallel_tasks; bunch++ ))
do
  echo "#!/bin/bash -l

###########################################################
#SBATCH --partition=cpu
#SBATCH --time=0-12:00
#SBATCH --nodes=1
#SBATCH --mem=8192
#SBATCH --ntasks=$nthreads
#SBATCH --job-name=MSM
#SBATCH --output=${workdir}/logs/HCP_to_template_bunch_${bunch}.txt
###########################################################" > $workdir/newMSM_HCP_to_template.sh

  from=$(( bunch * bunch_size + 1 ))
  to=$(( from + bunch_size - 1 ))

  echo "
while IFS=',' read -r subject group
do
  echo \"Registering subject \$subject to the template...\"

  time $HOME/fsldev/bin/newmsm \\
  --inmesh=$mesh \\
  --refmesh=$mesh \\
  --indata="$input_folder/\$subject.sulc.affine.ico6.shape.gii" \\
  --refdata=$template \\
  --conf=$config_file \\
  --out=$outdir/\$subject.MSMSulc.ico6. \\
  --verbose

done < <(sed -n \"${from},${to}p\" $sublist)
  " >> $workdir/newMSM_HCP_to_template.sh

  echo "Submitting job #$(( $bunch + 1 )) to slurm."
  sbatch $workdir/newMSM_HCP_to_template.sh

  if [ $? -ne 0 ];
  then
    exit 1
  fi

done
