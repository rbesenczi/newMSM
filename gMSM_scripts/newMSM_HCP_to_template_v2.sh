#!/bin/bash -l

#SBATCH --partition=cpu
#SBATCH --time=0-1:00
#SBATCH --nodes=1
#SBATCH --mem=2048
#SBATCH --ntasks=4
#SBATCH --job-name=MSM

while IFS=',' read -r subject
do
  echo "Registering subject $subject to the template..."

  time $HOME/fsldev/bin/newmsm \
    --inmesh=$mesh \
    --refmesh=$mesh \
    --indata="$input_folder/$subject.sulc.affine.ico6.shape.gii" \
    --refdata=$template \
    --conf=$config_file \
    --out=$outdir/$subject.MSMSulc.ico6. \
    --verbose

  wb_command -surface-distortion $mesh $outdir/$subject.MSMSulc.ico6.sphere.reg.surf.gii $outdir/$subject.MSMSulc.ico6.sphere.distortion.func.gii -local-affine-method -log2

done < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $sublist)
