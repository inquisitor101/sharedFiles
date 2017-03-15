#!/bin/bash -l

#The name of the script 
#SBATCH -J steady180

# specify he wall-clock time to be given to this job
#SBATCH -t 24:00:00 -A m.2016-1-440

# Number of nodes
#SBATCH --nodes=1

# Number of MPI processes per node (default 32)
#SBATCH --ntasks-per-node=32

# Number of cores to be allocated
NCORES=32

# Swap to PrgEnv-gnu 
module swap PrgEnv-cray PrgEnv-gnu

# Load the openfoam module 
module load openfoam/3.0.1

# Set the openfoam environment variables
. $FOAM_BASHRC

# Start the OpenFOAM job
aprun -n $NCORES pimpleFoam -parallel > logFiles/runAvg.log
aprun reconstructPar -latestTime 
