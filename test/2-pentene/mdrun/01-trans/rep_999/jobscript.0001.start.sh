#!/bin/bash

#===============================================================================
# Definitions for SLURM must come before any executables are called!
#===============================================================================
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

if [ -n "$NHOSTS" ]; then
    echo "Env. variable NHOSTS found, assuming we are using SGE."
    module load shared                           # access to modules in /cm/shared
else
    echo "Env. variable NHOSTS not set, assuming we are using SLURM."
    NHOSTS=$SLURM_JOB_NUM_NODES
    NSLOTS=$(( $SLURM_CPUS_ON_NODE * $NHOSTS ))
    echo Nhosts is $NHOSTS, Nslots is $NSLOTS, cores per node is $SLURM_CPUS_ON_NODE
fi

source /usr/local/gromacs/GMXRC2022

# use your own conda environment
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/home/chui/Software-2023-04-11/miniforge3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/chui/Software-2023-04-11/miniforge3/etc/profile.d/conda.sh" ]; then
        . "/home/chui/Software-2023-04-11/miniforge3/etc/profile.d/conda.sh"
    else
        export PATH="/home/chui/Software-2023-04-11/miniforge3/bin:$PATH"
    fi
fi
unset __conda_setup

if [ -f "/home/chui/Software-2023-04-11/miniforge3/etc/profile.d/mamba.sh" ]; then
    . "/home/chui/Software-2023-04-11/miniforge3/etc/profile.d/mamba.sh"
fi
# <<< conda initialize <<<
mamba activate pmx_MC_dev



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


# Get the number of GPUs we may use on this node
func.getNumberOfGpus ( )
{
    # CUDA_VISIBLE_DEVICES contains a comma-separated list of GPU IDs.
    # This list tells us how many and which GPUs we can use for the job.
    # (Other GPUs might be reserved for other jobs which we share the node with).
    # Find out how many entries (=GPUs) are in the list with the following line:
    echo $( echo $CUDA_VISIBLE_DEVICES | awk -F, '{print NF}' )

    return 0
}


# Write a host- or machinefile to be used by IntelMPI on Owl
# Thereby we control how many MPI processes are started on each of the nodes
# (which can be less than the number of slots if we use OpenMP threads
# Needs the desired number of MPI processes per node as only argument
func.writeHostfile ( )
{
    if [ $# -ne 1 ]; then
        echo "ERROR: g_submit's ${FUNCNAME[0]} needs the number of MPI processes per node as argument!" >&2
        echo "       it got: '$@'" >&2
        return 1
    fi

    # Local name for the hostfile:
    local PE_HOSTFILE_LOCAL="./pehostfile$JOB_ID"

    # Patch the PE_HOSTFILE to get the required number of MPI processes per node:
    awk -v var=$1 '{ print $1, var, $3, $4 }' $PE_HOSTFILE > $PE_HOSTFILE_LOCAL
    export PE_HOSTFILE=$PE_HOSTFILE_LOCAL
    return 0
}



# Guesses whether the node is exclusively used by this job
func.guessHaveExclusiveNode ( )
{
    local NCORES=$( cat /proc/cpuinfo | grep processor | wc -l )

    local bEXCLUSIVE=0

    if [ "$NHOSTS" -eq "1" ]; then
        # If we get slots on a single node, we might have to
        # share that node, unless we get all the available cores
        if [ "$NCORES" -eq "$NSLOTS" ]; then
            bEXCLUSIVE=1
        fi
    else
        # If we get more than a single node on Owl,
        # it is always for exclusive usage
        bEXCLUSIVE=1
    fi

    echo "$bEXCLUSIVE"
    return 0
}



# Sets various environment variables in a function
# so that this does not clutter up the rest of the job script
# Takes the number of threads per process as input argument
# The environment variables NSLOTS and NHOSTS are created by SGE
# NHOSTS is the number of physical nodes the job is running on
# NSLOTS is the number of available cores for the job.
func.setNeededVariables ( )
{
    if [ $# -ne 1 ]; then
        echo "ERROR: g_submit's ${FUNCNAME[0]} needs the number of threads per MPI rank as 1st argument!" >&2
        echo "       it got: '$@'" >&2
        return 1
    fi

    if [ -z $NHOSTS ]; then
        echo "ERROR: g_submit needs the environment variable NHOSTS to be set to the number of nodes assigned to the job." >&2
        return 1
    fi
    if [ -z $NSLOTS ]; then
        echo "ERROR: g_submit needs the environment variable NSLOTS to be set to the number of cores assigned to the job." >&2
        return 1
    fi

    # Set temporary directory:
    export MPD_CON_EXT="sge_$JOB_ID.$SGE_TASK_ID"
    export MPD_TMPDIR=$TMPDIR

    # The following variables are non-local and can be used later:
    if [ -z $NMPI ]; then  # set NMPI if not set already
        NMPI=$( echo "$NSLOTS / $1" | bc )
    fi
    if [ $( echo "$NMPI % $NHOSTS" | bc ) -ne "0" ] ; then
        echo "ERROR: g_submit refuses to distribute $NMPI MPI ranks on $NHOSTS hosts." >&2
        return 1
    fi

    if [ $NMPI -lt 1 ] ; then
        NMPI=1
    fi
    NMPI_PER_HOST=$( echo "$NMPI / $NHOSTS" | bc )

    echo "g_submit's job script was executed on host $( hostname )."
    echo "         In total there are $NSLOTS slots on $NHOSTS host(s) for this job."
    echo "         Will start the simulation using $NMPI_PER_HOST MPI rank(s) per node ($NMPI rank(s) in total)."
    echo "         Each MPI rank will have $1 OpenMP thread(s)."

    NGPU_PER_HOST=$( func.getNumberOfGpus )
    
    if [ -z $CUDA_VISIBLE_DEVICES ] ; then
        echo "         WARNING: Environment variable CUDA_VISIBLE_DEVICES is not set, although it should be!"
    else
        echo "         Will use $NGPU_PER_HOST GPU(s) per node (CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES)."
    fi

    return 0
}


#===============================================================================

# Do we have to share the node(s) with other users?
bEXCLUSIVE=$( func.guessHaveExclusiveNode )

NMPI=2  # multi-simulation with one MPI rank per replica
func.setNeededVariables 0

if (( $? )); then
    echo "func.setNeededVariables returned an exit code. Exiting." >&2
    exit 1
fi











SUF=$( func.getAccel )
func.writeHostfile $NMPI_PER_HOST
export MDRUN="$(which gmx_mpi$SUF) mdrun -ntomp 10" # set openmp accoring to your partition
GMXRUN="$MPIRUN -n $NMPI --cpu-bind no $MDRUN"
export GROMPP="$(which gmx_threads$SUF) grompp -maxwarn 1"

echo "Check software paths:"
which python
which mpirun
which gmx_mpi
which gmx
which nvcc
which pmx_mdrun
echo $GMXRUN
echo $GROMPP

pmx_mdrun \
-mdp_folder mdp/ -p ../../../topol.top -folder_start 000000 -cycle 200 \
-MDRUN "$GMXRUN" \
-GROMPP "$GROMPP" \
-maxh 5 
