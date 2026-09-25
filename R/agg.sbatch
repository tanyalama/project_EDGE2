#!/bin/bash
#SBATCH -A pi_tlama_smith_edu
#SBATCH -p cpu
#SBATCH -c 1
#SBATCH --mem=24G
#SBATCH -t 01:00:00
#SBATCH -o logs/%x_%j.out
export EDGE2_WS=/work/pi_tlama_smith_edu/project_EDGE2_v2
/work/pi_tlama_smith_edu/envs/edge2/bin/Rscript $EDGE2_WS/run/aggregate.R
