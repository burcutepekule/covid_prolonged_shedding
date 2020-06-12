#!/usr/bin/env bash
#SBATCH --time=23:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --output=OUT_SWITZERLAND_%A_%a.out
#SBATCH --array=0-152

module load anaconda3 
eval "$(conda shell.bash hook)"
conda activate COVID_IC_A7

rcmin=0
rcmax=100
steprc=2

mrelaxinmin=0
mrelaxinmax=10
stepmrelax=5

rcArr=()
mrelaxArr=()

for rc in $(seq ${rcmin} ${steprc} ${rcmax}); do
	for mrelax in $(seq ${mrelaxinmin} ${stepmrelax} ${mrelaxinmax}); do
		rcArr+=($rc)
		mrelaxArr+=($mrelax)
	done
done

rc=${rcArr[$SLURM_ARRAY_TASK_ID]}
mrelax=${mrelaxArr[$SLURM_ARRAY_TASK_ID]}

pathOut=result_SWITZERLAND_${rc}_${mrelax}.out

srun --nodes=1 --ntasks=1 --cpus-per-task=1 --exclusive \
  Rscript --vanilla analysis_SWITZERLAND.R ${rc} ${mrelax} > $pathOut