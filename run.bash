#!/bin/bash -l
#SBATCH -n 700                  # Number of cores
#SBATCH -t 2-00:00:00           # Runtime in D-HH:MM:SS
#SBATCH -J CVG
#SBATCH -o ./jobid_%j.out      # File to which STDOUT will be written
#SBATCH -e ./jobid_%j.err      # File to which STDERR will be written
#SBATCH -p cosma7               # Partition to submit to
#SBATCH -A dp004
#SBATCH --exclusive
#SBATCH --mail-user=christoph.becker@durham.ac.uk
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# perpare
module purge
module load intel_comp/2018-update2
module load intel_mpi/2018
# module load hdfview
unset I_MPI_HYDRA_BOOTSTRAP

# run
mpirun -np $SLURM_NTASKS ./bin/ramses3d ./namelist/cosmo_box1.nml

# report
echo ""
sacct -j $SLURM_JOBID --format=JobID,JobName,Partition,AllocCPUS,Elapsed,ExitCode

exit
