#!/bin/bash

# run a pmx_mdrun on one node using 1 GPU
#SBATCH --job-name=JOBNAME
#SBATCH --get-user-env
#SBATCH --gres=gpu:1             # number of GPUs requested
#SBATCH -p p20                   # partitions to use, each partition has one GPU
#SBATCH -t 5:30:00               # hours:min:sec
#SBATCH --signal TERM@120        # send a TERM signal 2 minutes before the queue time ends
#===============================================================================


# Do not ignore signals:
#$ -notify
#
# We install a signal handler for the case that this job in some emergency needs to
# be terminated (e.g. node or cluster shutdown). In such a case, both mdrun as well
# as this bash jobscript will get a TERM signal. We tell this jobscript to ignore the
# TERM signal (using the null "" action string below) so that mdrun has a few seconds
# to properly write a checkpoint file before it exits.
trap "" TERM

source /etc/profile.d/modules.sh             # enables us to use the 'module load' command
module load owl/intel-mpi-default            # for multi-node parallel runs using MPI


echo "Env. variable NHOSTS not set, assuming we are using SLURM."
NHOSTS=$SLURM_JOB_NUM_NODES
NSLOTS=$(( $SLURM_CPUS_ON_NODE * $NHOSTS ))
echo Nhosts is $NHOSTS, Nslots is $NSLOTS, cores per node is $SLURM_CPUS_ON_NODE

source /usr/local/gromacs/GMXRC2022

source /home/chui/Software-2023-04-11/miniforge3/bin/activate pmx_MC_dev


# Report the exit code of any pipeline as the exit code of the last program to return a non-zero exit code:
set -o pipefail

# Determine the optimal mdrun SIMD acceleration for GROMACS 2018 and later
func.getAccel ( )
{
    local CPUFLAGS=$( cat /proc/cpuinfo | grep flags     | head -1 )
    local VENDORID=$( cat /proc/cpuinfo | grep vendor_id | head -1 )
    local ACCEL

    if     [ $( echo $VENDORID | grep -c "AuthenticAMD" ) -eq "1" ]; then
        ACCEL="_AVX2_128"
    elif   [ $( echo $CPUFLAGS | grep -c "avx2"         ) -eq "1" ]; then
        ACCEL="_AVX2_256"
    elif   [ $( echo $CPUFLAGS | grep -c "avx"          ) -eq "1" ]; then
        ACCEL="_AVX_256"
    elif   [ $( echo $CPUFLAGS | grep -c "sse4_1"       ) -eq "1" ]; then
        ACCEL="_SSE4_1"
    elif   [ $( echo $CPUFLAGS | grep -c "sse2"         ) -eq "1" ]; then
        ACCEL="_SSE2"
    else
        ACCEL=""
    fi

    echo "$ACCEL"
}


#===============================================================================


NMPI=2  # multi-simulation with one MPI rank per replica
NOMP=$(( $NSLOTS / $NMPI ))

SUF=$( func.getAccel )
export MDRUN="$(which gmx_mpi$SUF) mdrun -ntomp $NOMP" # set openmp accoring to your partition
MPIMDRUN="$MPIRUN -n $NMPI --cpu-bind no $MDRUN"
export GROMPP="$(which gmx_threads$SUF) grompp -maxwarn 1"

echo "Check software paths:"
which nvcc
which pmx_mdrun
echo $MPIMDRUN
echo $GROMPP

pmx_mdrun \
-mdp_folder mdp/ -p ../../../topol.top -cycle 10 \
-MDRUN "$MPIMDRUN" \
-GROMPP "$GROMPP" \
-maxh 5 
