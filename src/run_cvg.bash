#!/bin/bash -l
#SBATCH -n 392                   # Number of cores
#SBATCH -t 0-13:00:00            # Runtime in D-HH:MM:SS
#SBATCH -J cvGv5
#SBATCH -o ./logs/jobid_%j.out  # File to which STDOUT will be written
#SBATCH -e ./logs/jobid_%j.err  # File to which STDERR will be written
#SBATCH -p cosma7               # Partition to submit to
#SBATCH -A dp004
#SBATCH --exclusive

# perpare
module purge
module load intel_comp/2018-update2 intel_mpi/2018
unset I_MPI_HYDRA_BOOTSTRAP

# run
mpirun -np $SLURM_NTASKS ./bin/ramses3d ./namelist/cosmo_FULL.nml

#test="tophat"  #[sine,gauss,tophat]
## store output regarding cv-Galileon in seperate file
#awk '/cvg/' ./logs/jobid_${SLURM_JOBID}.out >> ./logs/jobid_${SLURM_JOBID}_cvg_${test}.out &&
#awk '{$1=""}1' ./logs/jobid_${SLURM_JOBID}_cvg_${test}.out >> ./logs/jobid_${SLURM_JOBID}_cvg_1_${test}.out &&
#mv ./logs/jobid_${SLURM_JOBID}_cvg_1_${test}.out ./logs/jobid_${SLURM_JOBID}_cvg_${test}.out

#awk '/alpha_sf/' ./logs/jobid_${SLURM_JOBID}.out >> ./logs/jobid_${SLURM_JOBID}_alphasf.out &&
#awk '{$1=""}1' ./logs/jobid_${SLURM_JOBID}_alphasf.out >> ./logs/jobid_${SLURM_JOBID}_alphasf_1.out &&
#mv ./logs/jobid_${SLURM_JOBID}_alphasf_1.out ./logs/jobid_${SLURM_JOBID}_alphasf.out

# report
echo ""
sacct -j $SLURM_JOBID --format=JobID,JobName,Partition,AllocCPUS,Elapsed,ExitCode

exit
