
export OMP_NUM_THREADS=$(($(nproc) / 2))
echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"

base=$PWD
# This box has +1 charge, use -maxwarn 1
cd $base/nvt/0
gmx grompp -f eq.mdp -c ../../../confout.gro -p $base/../topol.top -maxwarn 1 -o nvt > grompp.log 2>&1
cd $base/nvt/1
gmx grompp -f eq.mdp -c ../../../confout.gro -p $base/../topol.top -maxwarn 1 -o nvt > grompp.log 2>&1
cd $base/nvt
mpirun -np 2 --bind-to none gmx_mpi mdrun -s nvt.tpr -deffnm nvt -multidir 0 1 -maxh 0.5 > mdrun.log 2>&1

# prepare the eq.tpr for the first cycle
cd $base
for l in 0 1
do
    mkdir $base/000000/$l -p
    cd    $base/000000/$l
    gmx grompp -f ../../mdp/eq$l.mdp \
        -c ../../nvt/$l/nvt.gro \
        -t ../../nvt/$l/nvt.cpt \
        -p $base/../topol.top \
        -o eq -maxwarn 1 > grompp.log 2>&1
done

cd $base
pmx_mdrun -MDRUN "mpirun -np 2 --bind-to none gmx_mpi mdrun -ntomp $OMP_NUM_THREADS" -GROMPP "gmx grompp -maxwarn 1" -mdp_folder ./mdp -p $base/../topol.top -cycle 3     # run 5 cycles, 0-4
pmx_mdrun -MDRUN "mpirun -np 2 --bind-to none gmx_mpi mdrun -ntomp $OMP_NUM_THREADS" -GROMPP "gmx grompp -maxwarn 1" -mdp_folder ./mdp -p $base/../topol.top -cyc_until 6 # run until 6 cycles, 0-5
pmx_mdrun -MDRUN "mpirun -np 2 --bind-to none gmx_mpi mdrun -ntomp $OMP_NUM_THREADS" -GROMPP "gmx grompp -maxwarn 1" -mdp_folder ./mdp -p $base/../topol.top -cycle 4     # append 4 more cycles, 0-9
ls 000009/*
analysis_bar
