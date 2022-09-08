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

ID='4U0P'

paramsFile=../scripts/params4marylou/${ID}_params.json

#Second, to read in a passedCLS file to set structures to be designed
cd ../scripts
for i in $(seq 0 $((PDBID_No-1)))
do
    #Get PDBID
    ID=${PDBIDs[$i]%/*}
    ID=${ID##*/}

    #Open passed CLS check master file
    #CLSfile=../${ID}/passedCLS/${ID}_CLUSTERX_01_passedCLS.txt
    paramsFile=../scripts/params4marylou/${ID}_params.json
    python3.7 SCF_Check_Driver.py $paramsFile ${SLURM_ARRAY_TASK_ID}

done


python3.7 SCF_Check_analyzer.py 
