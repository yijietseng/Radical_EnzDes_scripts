#!/bin/bash

#SBATCH --time=120:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=5120M   # memory per CPU core
#SBATCH -J "4U0P_hp_design"   # job name
#SBATCH --mail-user=yijietseng@gmail.com   # email address
#SBATCH --mail-type=END
#SBATCH --array=1-8

#Make sure to change the numbers of array according to the passed number of hbDesign

# Compatibility variables for PBS. Delete if not needed.
export PBS_NODEFILE=`/fslapps/fslutils/generate_pbs_nodefile`
export PBS_JOBID=$SLURM_JOB_ID
export PBS_O_WORKDIR="$SLURM_SUBMIT_DIR"
export PBS_QUEUE=batch

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
source /fslhome/tseng3/.bashrc

cd ../scripts

ID=$1

paramsFile=params4marylou/${ID}_params.json

HBfile=../${ID}/passedHB/${ID}_CLUSTERX_01_passedHB.txt

#struc=`cat 4M7T_CLUSTERX_02_passedCLS.txt | awk -v "line_start=${SLURM_ARRAY_TASK_ID}" -v "line_end=${SLURM_ARRAY_TASK_ID}" 'NR==line_start,NR==line_end {print $4}'`

python3.7 Design_driver_marylou.py $paramsFile $HBfile hp ${SLURM_ARRAY_TASK_ID}
