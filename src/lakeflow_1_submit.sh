#!/bin/bash

###########################################################################
## environment & variable setup
####### job customization
#SBATCH --job-name dwnld
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5GB
#SBATCH --time 5:00:00
#SBATCH -p normal_q
#SBATCH -A swot
#SBATCH --array=1-50
####### end of job customization
# end of environment & variable setup

module load apptainer/1.4.0
apptainer exec \
    --pwd /projects/swot/hana/LakeFlow_Confluence \
    --bind /projects/swot/hana/LakeFlow_Confluence \
    --cleanenv \
    --env SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID \
    /projects/swot/hana/LakeFlow_Confluence/lakeflow_input.sif Rscript src/lakeflow_1.R -c "in/lakeids/lakeid${SLURM_ARRAY_TASK_ID}.csv" -w 1 -i in/ --index -256