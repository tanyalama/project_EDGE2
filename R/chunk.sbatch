#!/bin/bash
#SBATCH -A pi_tlama_smith_edu
#SBATCH -p cpu
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -t 04:00:00
#SBATCH --array=1-20
#SBATCH -o logs/%x_%A_%a.out
export EDGE2_WS=/work/pi_tlama_smith_edu/project_EDGE2_v2 CHUNK_SIZE=50 EXT_PROB=Isaac BASE_SEED=20240601
/work/pi_tlama_smith_edu/envs/edge2/bin/Rscript $EDGE2_WS/run/run_chunk.R $SLURM_ARRAY_TASK_ID
