#!/bin/bash

#SBATCH --time=120:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=5120M   # memory per CPU core
#SBATCH -J "4U0P_hb_design"   # job name
#SBATCH --mail-user=yijietseng@gmail.com   # email address
#SBATCH --mail-type=END
#SBATCH --array=1-375

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

ID='4U0P'

paramsFile=params4marylou/${ID}_params.json

CLSfile=../${ID}/passedCLS/${ID}_CLUSTERX_01_passedCLS.txt

python3.7 Design_driver_marylou.py $paramsFile $CLSfile hb ${SLURM_ARRAY_TASK_ID}