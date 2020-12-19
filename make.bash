# perpare
module purge
module load intel_comp/2018-update2 intel_mpi/2018
# module load hdfview
unset I_MPI_HYDRA_BOOTSTRAP

cd ./src/bin
# compile
make clean && make
