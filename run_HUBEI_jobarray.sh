#!/usr/bin/env bash
#SBATCH --time=23:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --output=OUT_HUBEI_%A_%a.out
#SBATCH --array=0-713

mkdir /net/cephfs/data/btepek/COVID_IC/IC_OUT/OUT_HUBEI/
mkdir /net/cephfs/data/btepek/COVID_IC/IC_OUT/OUT_HUBEI/CSVS/
mkdir /net/cephfs/data/btepek/COVID_IC/IC_OUT/OUT_HUBEI/FIGS/
mkdir /net/cephfs/data/btepek/COVID_IC/IC_OUT/OUT_HUBEI/RDATA/

module load anaconda3 
eval "$(conda shell.bash hook)"
conda activate COVID_IC_A6

rcmin=0
rcmax=100
steprc=2


endindexinmin=1
endindexinmax=14
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
pathOut=result_HUBEI_${rc}_${mrelax}_${endindex}.out

srun --nodes=1 --ntasks=1 --cpus-per-task=1 --exclusive \
  Rscript --vanilla analysis_HUBEI.R ${rc} ${mrelax} ${endindex} > $pathOut