#!/usr/bin/env bash
#SBATCH --time=23:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --output=OUT_ITALY_%A_%a.out
#SBATCH --array=0-305

mkdir /net/cephfs/data/btepek/COVID_IC/IC_OUT/OUT_ITALY/
mkdir /net/cephfs/data/btepek/COVID_IC/IC_OUT/OUT_ITALY/CSVS/
mkdir /net/cephfs/data/btepek/COVID_IC/IC_OUT/OUT_ITALY/FIGS/
mkdir /net/cephfs/data/btepek/COVID_IC/IC_OUT/OUT_ITALY/RDATA/

module load anaconda3 
eval "$(conda shell.bash hook)"
conda activate COVID_IC_A1

rcmin=0
rcmax=100
steprc=2


endindexinmin=1
endindexinmax=6
endindexinstep=1

mrelax=0

rcArr=()
endIndexArr=()

for rc in $(seq ${rcmin} ${steprc} ${rcmax}); do
	for endindex in $(seq ${endindexinmin} ${endindexinstep} ${endindexinmax}); do
		rcArr+=($rc)
		endIndexArr+=($endindex)
	done
done

rc=${rcArr[$SLURM_ARRAY_TASK_ID]}
endindex=${endIndexArr[$SLURM_ARRAY_TASK_ID]}
pathOut=result_ITALY_${rc}_${mrelax}_${endindex}.out

srun --nodes=1 --ntasks=1 --cpus-per-task=1 --exclusive \
  Rscript --vanilla analysis_ITALY.R ${rc} ${mrelax} ${endindex} > $pathOut