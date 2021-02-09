#!/bin/bash -l
#SBATCH -n 700                  # Number of cores
#SBATCH -t 2-00:00:00           # Runtime in D-HH:MM:SS
#SBATCH -J CVG
#SBATCH -o ./jobid_%j.out       # File to which STDOUT will be written
#SBATCH -e ./jobid_%j.err       # File to which STDERR will be written
#SBATCH -p cosma7               # Partition to submit to
#SBATCH -A dp004
#SBATCH --exclusive

# perpare
module purge
module load intel_comp/2018-update2 intel_mpi/2018

# run
mpirun -np $SLURM_NTASKS ./src/bin/ramses3d ./src/namelist/cosmos.nml

# report
echo ""
sacct -j $SLURM_JOBID --format=JobID,JobName,Partition,AllocCPUS,Elapsed,ExitCode

exit
